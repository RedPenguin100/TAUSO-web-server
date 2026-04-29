import logging
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from typing import Optional
from email_service import send_processing_started, send_processing_completed
from tauso.aso_generation import generate_aso_features, generate_stub_data
from tauso.populate.calculators.cache import AssetCache
from tauso.data.consts import CANONICAL_GENE, SEQUENCE
import os
import uuid
import pandas as pd
from tauso.off_target.search import find_all_gene_off_targets_BULK
import tempfile


# Force Python to show INFO logs in the Docker terminal
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# Initialize a single isolated background process queue
executor = ProcessPoolExecutor(max_workers=1)

@dataclass
class JobConfig:
    """A clean struct to hold all mandatory and optional pipeline parameters."""
    target_data: str
    target_mrna_name: str
    source_info: str
    user_email: str
    cell_line: Optional[str] = None
    include_feature_breakdown: bool = False


def generate_off_target_table(df, target_gene_name, genome="GRCh38", threads=8):
    """
    Runs bulk off-target alignment and returns a completely separate
    long-format DataFrame of off-target hits. Uses OS temp dir to avoid permission errors.
    """
    if df.empty or SEQUENCE not in df.columns:
        return pd.DataFrame(columns=[SEQUENCE, "off_target_gene"])

    unique_seqs = df[SEQUENCE].dropna().unique()

    # 1. Create a secure temporary file in the OS temp directory (e.g., /tmp/)
    fd, temp_fasta = tempfile.mkstemp(suffix=".fasta", prefix="tauso_web_bulk_")

    try:
        # Write to the temp file
        with os.fdopen(fd, "w") as f:
            for seq in unique_seqs:
                dna = seq.upper().replace("U", "T")
                f.write(f">{dna}\n{dna}\n")

        # 2. Run alignment
        seq_to_genes_map = find_all_gene_off_targets_BULK(temp_fasta, genome, threads)

    finally:
        # 3. Cleanup: guaranteed to run even if the alignment crashes
        if os.path.exists(temp_fasta):
            os.remove(temp_fasta)

    # 4. Build the separate relational table
    records = []
    for seq in unique_seqs:
        dna = seq.upper().replace("U", "T")
        hits = seq_to_genes_map.get(dna, [])

        # Isolate true off-targets
        off_targets = [g for g in hits if g != target_gene_name]
        for ot in off_targets:
            records.append({
                SEQUENCE: seq,
                "off_target_gene": ot
            })

    return pd.DataFrame(records)

def execute_tauso_pipeline(config: JobConfig):
    """
    The main pipeline execution function that runs in an isolated process.
    """
    logger.info(f"Process started for {config.user_email} | Cell Line: {config.cell_line}")
    send_processing_started(config.user_email, config.source_info)

    try:
        # 1. INITIALIZE CACHE
        cache = AssetCache()

        # Check if a custom sequence was provided via File Upload or Paste
        # (If they selected a DB gene, target_data will be empty)
        if config.target_data:
            cache.set_custom_gene(name=config.target_mrna_name, sequence=config.target_data)
            sequence_to_pass = config.target_data
        else:
            sequence_to_pass = None

        logger.info(f"Generating stub data for {config.target_mrna_name}...")

        # 2. GENERATE STUB DATA
        # Note: I removed 'first_n=10' so it processes the whole sequence,
        # but you can add it back if you strictly want to limit the pipeline.
        stub_data = generate_stub_data(
            target_gene=config.target_mrna_name,
            gene_sequence=sequence_to_pass,
            version=f"{config.target_mrna_name}_v1",
            first_n=10
        )

        logger.info(f"Generating ASO features...")

        # 3. GENERATE FEATURES & PREDICTIONS
        df, aso_features = generate_aso_features(stub_data, cache)

        # ---> NEW INTEGRATION HERE <---
        logger.info(f"Generating off-targets table...")

        df_off_targets = generate_off_target_table(df, target_gene_name=config.target_mrna_name)

        # 4. PACKAGE RESULTS FOR EMAIL
        csv_bytes_main = df.to_csv(index=False).encode('utf-8')
        csv_bytes_ot = df_off_targets.to_csv(index=False).encode('utf-8')

        results_files = [
            (f"{config.target_mrna_name}_aso_predictions.csv", csv_bytes_main),
            (f"{config.target_mrna_name}_off_targets.csv", csv_bytes_ot)
        ]

        # Send the success email with the attached CSVs
        send_processing_completed(config.user_email, config.source_info, results_files)
        logger.info(f"Process completed successfully for {config.user_email}")

    except Exception as e:
        logger.error(f"Pipeline crashed for {config.user_email}: {str(e)}")

def trigger_background_job(config: JobConfig):
    """
    Submits the configuration object to the ProcessPool.
    Returns instantly, leaving the UI perfectly responsive.
    """
    executor.submit(execute_tauso_pipeline, config)