import logging
import os
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from typing import Optional

from email_service import send_processing_completed, send_processing_started
from tauso.aso_generation import design_asos, summarize_design, tox_details

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# design_asos runs inside the daemonic ProcessPool worker below, which cannot spawn its own child
# processes, so feature computation stays single-process. Raising TAUSO_DESIGN_JOBS above 1 requires
# a non-daemonic executor.
DESIGN_JOBS = int(os.environ.get("TAUSO_DESIGN_JOBS", "1"))
# Featurizing tiled candidates is the dominant cost and scales with transcript length, so
# TAUSO_FIRST_N bounds it. design_asos tiles 5'->3' and takes the first N windows, so this covers
# only the 5' end of the target (N=50 twenty-mers spans positions 1-69), not a sample across it.
# Set TAUSO_FIRST_N=0 to tile the whole transcript, which takes hours for a full-length mRNA.
FIRST_N = int(os.environ.get("TAUSO_FIRST_N", "50")) or None
# Report only the top-ranked shortlist so the per-candidate bowtie off-target search stays bounded.
TOP_N = int(os.environ.get("TAUSO_TOP_N", "100"))
# Mismatch tolerance for the sequence off-target search (0 = perfect matches only).
OFFTARGET_MAX_DISTANCE = int(os.environ.get("TAUSO_OFFTARGET_MAX_DISTANCE", "2"))

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


def execute_tauso_pipeline(config: JobConfig):
    """Design ASOs for the target end-to-end and email the ranked results, safety detail, and
    per-candidate sequence off-target hits. Runs in an isolated background process."""
    logger.info(
        f"Design job for {config.user_email} | gene={config.target_mrna_name} | cell_line={config.cell_line}"
    )
    send_processing_started(config.user_email, config.source_info)

    try:
        # Tile candidate ASOs across the target, featurize them, and score with the bundled model.
        # A DB-gene selection leaves target_data empty -> the target is looked up from the genome cache.
        ranked, off_targets = design_asos(
            config.target_mrna_name,
            gene_sequence=(config.target_data or None),
            cell_line=config.cell_line,
            aso_sizes=[20],
            first_n=FIRST_N,
            top_n=TOP_N,
            n_jobs=DESIGN_JOBS,
            off_targets=True,
            off_target_max_distance=OFFTARGET_MAX_DISTANCE,
        )
        logger.info(f"Ranked {len(ranked)} candidate ASOs; building result tables...")

        designed = summarize_design(ranked)
        safety = tox_details(ranked)

        name = config.target_mrna_name
        results_files = [
            (f"{name}_designed_asos.csv", designed.to_csv(index=False).encode("utf-8")),
            (f"{name}_safety_detail.csv", safety.to_csv(index=False).encode("utf-8")),
            (f"{name}_off_targets.csv", off_targets.to_csv(index=False).encode("utf-8")),
        ]

        send_processing_completed(config.user_email, config.source_info, results_files)
        logger.info(f"Design job complete for {config.user_email}.")

    except Exception as e:
        logger.exception(f"Design job failed for {config.user_email}: {e}")


def trigger_background_job(config: JobConfig):
    """Submit the config to the background pool and return immediately."""
    executor.submit(execute_tauso_pipeline, config)
