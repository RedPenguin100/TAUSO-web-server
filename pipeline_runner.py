import logging
import os
import tempfile
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from typing import Optional

import pandas as pd

from email_service import send_processing_completed, send_processing_started
from tauso.aso_generation import design_asos, design_details, summarize_design
from tauso.data.consts import ASO_SEQUENCE
from tauso.off_target.search import find_all_gene_off_targets_BULK

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# design_asos runs inside the daemonic ProcessPool worker below, which cannot spawn its own child
# processes, so feature computation stays single-process. Raising TAUSO_DESIGN_JOBS above 1 requires
# a non-daemonic executor. TAUSO_FIRST_N tiles the whole transcript when unset; set it to featurize
# only the first N candidates (faster, partial).
DESIGN_JOBS = int(os.environ.get("TAUSO_DESIGN_JOBS", "1"))
FIRST_N = int(os.environ["TAUSO_FIRST_N"]) if os.environ.get("TAUSO_FIRST_N") else None
OFFTARGET_THREADS = int(os.environ.get("TAUSO_BOWTIE_THREADS", "4"))

# Single isolated background process queue: submit returns instantly, keeping the UI responsive.
executor = ProcessPoolExecutor(max_workers=1)


@dataclass
class JobConfig:
    """All mandatory and optional parameters for one design job."""

    target_data: str
    target_mrna_name: str
    source_info: str
    user_email: str
    cell_line: Optional[str] = None


def generate_off_target_table(ranked, target_gene_name, genome="GRCh38", threads=OFFTARGET_THREADS):
    """Long-format table of the genes each candidate ASO aligns to (bulk bowtie), excluding the
    intended target. Complements the aggregate off-target burden in the safety detail."""
    if ranked.empty or ASO_SEQUENCE not in ranked.columns:
        return pd.DataFrame(columns=[ASO_SEQUENCE, "off_target_gene"])

    unique_seqs = ranked[ASO_SEQUENCE].dropna().unique()

    fd, temp_fasta = tempfile.mkstemp(suffix=".fasta", prefix="tauso_web_bulk_")
    try:
        with os.fdopen(fd, "w") as f:
            for seq in unique_seqs:
                dna = seq.upper().replace("U", "T")
                f.write(f">{dna}\n{dna}\n")
        seq_to_genes_map = find_all_gene_off_targets_BULK(temp_fasta, genome, threads)
    finally:
        if os.path.exists(temp_fasta):
            os.remove(temp_fasta)

    records = []
    for seq in unique_seqs:
        dna = seq.upper().replace("U", "T")
        for gene in seq_to_genes_map.get(dna, []):
            if gene != target_gene_name:
                records.append({ASO_SEQUENCE: seq, "off_target_gene": gene})
    return pd.DataFrame(records)


def execute_tauso_pipeline(config: JobConfig):
    """Design ASOs for the target end-to-end and email the ranked results, safety detail, and
    off-target gene hits. Runs in an isolated background process."""
    logger.info(
        f"Design job for {config.user_email} | gene={config.target_mrna_name} | cell_line={config.cell_line}"
    )
    send_processing_started(config.user_email, config.source_info)

    try:
        # Tile candidate ASOs across the target, featurize them, and score with the bundled model.
        # A DB-gene selection leaves target_data empty -> the target is looked up from the genome cache.
        ranked = design_asos(
            config.target_mrna_name,
            gene_sequence=(config.target_data or None),
            cell_line=config.cell_line,
            aso_sizes=[20],
            first_n=FIRST_N,
            n_jobs=DESIGN_JOBS,
        )
        logger.info(f"Ranked {len(ranked)} candidate ASOs; building result tables...")

        designed = summarize_design(ranked)
        safety = design_details(ranked)
        off_targets = generate_off_target_table(ranked, target_gene_name=config.target_mrna_name)

        name = config.target_mrna_name
        results_files = [
            (f"{name}_designed_asos.csv", designed.to_csv(index=False).encode("utf-8")),
            (f"{name}_safety_detail.csv", safety.to_csv(index=False).encode("utf-8")),
            (f"{name}_off_target_genes.csv", off_targets.to_csv(index=False).encode("utf-8")),
        ]

        send_processing_completed(config.user_email, config.source_info, results_files)
        logger.info(f"Design job complete for {config.user_email}.")

    except Exception as e:
        logger.exception(f"Design job failed for {config.user_email}: {e}")


def trigger_background_job(config: JobConfig):
    """Submit the config to the background pool and return immediately."""
    executor.submit(execute_tauso_pipeline, config)
