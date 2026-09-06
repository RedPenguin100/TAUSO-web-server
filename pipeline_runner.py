import logging
import os
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool
from dataclasses import dataclass
from typing import Optional

import jobs
from email_service import send_processing_completed, send_processing_failed, send_processing_started
from tauso.aso_generation import (
    _sequence_offtarget_table,
    default_config,
    design_asos,
    summarize_design,
    tox_details,
)

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

# design_asos runs inside the daemonic ProcessPool worker below, which cannot spawn its own child
# processes, so feature computation stays single-process. Raising TAUSO_DESIGN_JOBS above 1 requires
# a non-daemonic executor.
DESIGN_JOBS = int(os.environ.get("TAUSO_DESIGN_JOBS", "1"))
# TAUSO_FIRST_N bounds how many tiled candidates are featurised. design_asos tiles 5'->3' and
# takes the first N windows, so this covers the 5' end of the target rather than sampling across
# it. Most of a run is fixed cost -- 50 candidates take 87 seconds and 500 take 121 -- so the
# marginal candidate is cheap and the bound can be generous. Set TAUSO_FIRST_N=0 to tile the whole
# transcript, which for a 8.8 kb target is around a quarter of an hour.
FIRST_N = int(os.environ.get("TAUSO_FIRST_N", "500")) or None
# How many candidates make the shortlist: the table, and the per-candidate bowtie off-target
# search, which is what this bounds. Every candidate scored is still charted and downloadable --
# a plot of the winners alone reads as a landscape while hiding exactly the stretches a reader
# most needs to see are bad.
TOP_N = int(os.environ.get("TAUSO_TOP_N", "100"))
# Mismatch tolerance for the sequence off-target search (0 = perfect matches only).
OFFTARGET_MAX_DISTANCE = int(os.environ.get("TAUSO_OFFTARGET_MAX_DISTANCE", "2"))

# Single isolated background process queue: submit returns instantly, keeping the UI responsive.
# It is created on demand rather than held for the life of the process: a worker that dies abruptly
# -- an OOM kill, or a segfault in one of the native libraries -- leaves the pool permanently broken,
# while Streamlit itself stays up and healthy, so nothing else would notice.
_executor = None


def _get_executor():
    global _executor
    if _executor is None:
        _executor = ProcessPoolExecutor(max_workers=1)
    return _executor


# Gapmer chemistries offered to the user. `pattern` is TAUSO's per-sugar code -- 'M' 2'-MOE,
# 'C' cEt, 'd' deoxy -- and its length is the ASO length, so each chemistry designs the oligo
# length its wing geometry implies: 5-10-5 for 2'-MOE, 3-10-3 for cEt. `ps_pattern` is one
# character per inter-nucleotide linkage, so one shorter than the oligo; '*' is phosphorothioate.
# `modification` gates the MOE-specific hybridization features, which look for "MOE" in it,
# while the cEt features key off 'C' in the pattern.
CHEMISTRIES = {
    "2'-MOE": {
        "pattern": "MMMMMddddddddddMMMMM",
        "modification": "MOE/5-methylcytosines/deoxy",
    },
    "cEt": {
        "pattern": "CCCddddddddddCCC",
        "modification": "cEt/5-methylcytosines/deoxy",
    },
}
for _spec in CHEMISTRIES.values():
    _spec["ps_pattern"] = "*" * (len(_spec["pattern"]) - 1)

DEFAULT_CHEMISTRY = "2'-MOE"

# Sugar codes TAUSO understands, and the two linkage codes: '*' phosphorothioate, 'o' phosphodiester.
SUGAR_CODES = "MCLd"
SUGAR_NAMES = {"M": "2'-MOE", "C": "cEt", "L": "LNA", "d": "DNA"}
LINKAGE_CODES = "*o"
# The model is not calibrated outside the ASO lengths seen in training.
ASO_LENGTH_RANGE = (12, 28)
# The lengths offered for design. Narrower than what the model accepts: a gapmer keeps a 10-nt
# DNA gap, so these are the lengths whose wings come out at 3-7 nt either side.
DESIGN_LENGTH_RANGE = (16, 24)
DNA_GAP = 10

# Delivery methods offered. TAUSO also accepts "Other", which is the catch-all the training data
# uses for cohorts whose method was never recorded -- not something a user can meaningfully pick.
TRANSFECTION_METHODS = ["Gymnosis", "Electroporation", "Lipofection"]

# Assay conditions are model inputs, so they are bounded by what the training experiments covered
# and default to the median of that distribution rather than to the low edge.
DOSAGE_RANGE_NM = (2, 20000)
DEFAULT_DOSAGE_NM = 4000
CELL_DENSITY_RANGE = (85, 300000)
DEFAULT_CELL_DENSITY = 20000

# Two features carried into the results so a reader can see why a site scored as it did.
# Accessibility is how unpaired the target site is over a 60-nt window; the hybridization term is
# the DNA:RNA duplex free energy, which is computed for every chemistry rather than only for MOE.
# A 20-nt window opened over the target site, folded with 60 nt of flank on each side.
ACCESSIBILITY_FEATURE = "access_f60_sinf_u20_a5"
HYBRIDIZATION_FEATURE = "hybr_dna_rna_dg"
RNASE_FEATURE = "rnase_score_dinucleotide_R4a_dinuc_dynamic"
# Folding energy of the site itself. The narrow window is the one that varies candidate to
# candidate, and the one the model leans on; the wide windows barely move along a transcript.
MFE_FEATURE = "fold_mfe_win25_flank30_step4"

# design_asos tiles a window at every position of the target before any of them are scored, so the
# target length sets how much memory the job needs up front. The longest human mRNA is around
# 110 kb, so this accepts any real transcript while keeping a pasted mistake from taking the worker
# down -- which would break the pool for every job after it.
MAX_TARGET_LENGTH = 200000


def exon_layout(locus) -> dict:
    """Exons as offsets into the pre-mRNA, which is the coordinate a candidate's start is given in.

    The annotation holds genomic coordinates, and the pre-mRNA runs from gene_start to gene_end, so
    a positive-strand exon is measured from the start and a negative-strand one from the end."""
    length = len(locus.full_mrna) if locus.full_mrna else 0
    reverse = str(getattr(locus, "strand", "")).endswith("NEG") or getattr(locus, "strand", 1) == -1
    exons = []
    for start, end in getattr(locus, "_exon_indices", []):
        if reverse:
            exons.append([locus.gene_end - end, locus.gene_end - start])
        else:
            exons.append([start - locus.gene_start, end - locus.gene_start])
    return {"length": length, "exons": sorted(exons)}


def describe_chemistry(chemical_pattern: str, ps_pattern: str) -> str:
    """The chemistry in the terms it is normally written: wing-gap-wing, the modified sugar, and
    how much of the backbone is phosphorothioate. A length on its own is not a chemistry."""
    runs = [len(r) for r in chemical_pattern.replace("d", " ").split()]
    deoxy = len(chemical_pattern) - sum(runs)
    geometry = f"{runs[0]}-{deoxy}-{runs[-1]}" if len(runs) >= 2 else f"{len(chemical_pattern)}-mer"
    sugar = next((SUGAR_NAMES[c] for c in chemical_pattern if c != "d"), "DNA")
    thio = ps_pattern.count("*")
    backbone = "full PS" if thio == len(ps_pattern) else f"{thio}/{len(ps_pattern)} PS"
    return f"{geometry} {sugar}, {backbone}"


def to_idt_notation(sequence: str, chemical_pattern: str, ps_pattern: str) -> Optional[str]:
    """The IDT order string for one designed ASO, or None when the chemistry has no IDT equivalent.

    TAUSO renders these; cEt is not an IDT catalogue product, so a cEt oligo has no order string and
    raises there rather than returning something unorderable."""
    from tauso.common.modifications import to_idt_notation as render

    modification = f"{'MOE' if 'M' in chemical_pattern else 'LNA'}/5-methylcytosines/deoxy"
    try:
        return render(sequence, chemical_pattern, ps_pattern, modification)
    except ValueError:
        return None


def describe_pattern_problem(chemical_pattern: str, ps_pattern: str) -> Optional[str]:
    """Explain why this sugar/backbone pair cannot be designed, or None if it can. Several features
    return NaN rather than failing when the pattern does not line up with the oligo, so the pair is
    checked here instead of being discovered as a blank column in the results."""
    from tauso.common.modifications import is_gapmer

    low, high = ASO_LENGTH_RANGE
    if not low <= len(chemical_pattern) <= high:
        return f"The sugar pattern is {len(chemical_pattern)} long; the model covers {low}–{high}."
    unknown = sorted(set(chemical_pattern) - set(SUGAR_CODES))
    if unknown:
        return f"The sugar pattern may only use {', '.join(SUGAR_CODES)} — found {', '.join(unknown)}."
    if not is_gapmer(chemical_pattern):
        return "The sugar pattern must be a gapmer: a run of d flanked by modified sugars on both sides."
    if len(ps_pattern) != len(chemical_pattern) - 1:
        return (
            f"The backbone describes the bonds between sugars, so it must be "
            f"{len(chemical_pattern) - 1} long, not {len(ps_pattern)}."
        )
    unknown = sorted(set(ps_pattern) - set(LINKAGE_CODES))
    if unknown:
        return f"The backbone may only use {' or '.join(LINKAGE_CODES)} — found {', '.join(unknown)}."
    return None


@dataclass
class JobConfig:
    """All mandatory and optional parameters for one design job."""

    target_data: str
    target_mrna_name: str
    source_info: str
    user_email: str
    job_id: Optional[str] = None
    cell_line: Optional[str] = None
    chemical_pattern: str = CHEMISTRIES[DEFAULT_CHEMISTRY]["pattern"]
    ps_pattern: str = CHEMISTRIES[DEFAULT_CHEMISTRY]["ps_pattern"]
    # Left unset, these reach the model as missing. Around an eighth of the training experiments
    # record none of them, so the booster has a branch for each.
    transfection: Optional[str] = None
    dosage_nm: Optional[int] = None
    cell_density: Optional[int] = None

    @property
    def modification(self) -> str:
        """The MOE hybridization features look for "MOE" in this string while the cEt features key
        off 'C' in the sugar pattern, so it is derived from the pattern rather than chosen apart."""
        wings = set(self.chemical_pattern) - {"d"}
        return f"{'MOE' if 'M' in wings else 'cEt' if 'C' in wings else 'LNA'}/5-methylcytosines/deoxy"


def _layout_for(config: JobConfig):
    """The gene model for this job's target, or None for a sequence the user supplied."""
    if config.target_data:
        return None
    try:
        from tauso.populate.calculators.cache import AssetCache

        locus = AssetCache(genome="GRCh38").get_full_gene_data().get(config.target_mrna_name)
        return exon_layout(locus) if locus else None
    except Exception:
        logger.warning("Could not read the gene model for %s", config.target_mrna_name)
        return None


def execute_tauso_pipeline(config: JobConfig):
    """Design ASOs for the target end-to-end and email the ranked results, safety detail, and
    per-candidate sequence off-target hits. Runs in an isolated background process."""
    logger.info(
        f"Design job for {config.user_email} | gene={config.target_mrna_name} | "
        f"cell_line={config.cell_line} | sugars={config.chemical_pattern} | "
        f"transfection={config.transfection} | {config.dosage_nm} nM | {config.cell_density} cells/well"
    )
    jobs.mark(config.job_id, jobs.RUNNING)
    send_processing_started(config.user_email, config.source_info)

    try:
        design_config = default_config()
        design_config.standard_chemical_pattern = config.chemical_pattern
        design_config.standard_ps_pattern = config.ps_pattern
        design_config.standard_modification = config.modification
        # An unrecognised transfection label one-hot encodes to NaN, which is how "not recorded"
        # is spelled for all three of these.
        design_config.transfection_method = config.transfection
        design_config.volume = float("nan") if config.dosage_nm is None else config.dosage_nm
        design_config.cell_per_well = float("nan") if config.cell_density is None else config.cell_density

        # Tile candidate ASOs across the target, featurize them, and score with the bundled model.
        # A DB-gene selection leaves target_data empty -> the target is looked up from the genome cache.
        # The oligo length comes from the sugar pattern: several features return NaN unless the
        # pattern is exactly as long as the ASO.
        layout = None
        # Scored without a cutoff, so the whole scan is available to chart; the shortlist below
        # is what the off-target search and the table are bounded to.
        ranked = design_asos(
            config.target_mrna_name,
            gene_sequence=(config.target_data or None),
            cell_line=config.cell_line,
            aso_sizes=[len(config.chemical_pattern)],
            config=design_config,
            first_n=FIRST_N,
            top_n=None,
            n_jobs=DESIGN_JOBS,
            off_targets=False,
        )
        logger.info(f"Scored {len(ranked)} candidate ASOs; building result tables...")

        shortlist = ranked.head(TOP_N)
        off_targets = _sequence_offtarget_table(
            shortlist,
            genome="GRCm39" if design_config.organism_name == "mouse" else "GRCh38",
            max_distance=OFFTARGET_MAX_DISTANCE,
            exclude_genes=None,
        )
        logger.info(f"{len(off_targets)} off-target hits across the top {len(shortlist)}.")

        designed = summarize_design(ranked)
        for column in (ACCESSIBILITY_FEATURE, MFE_FEATURE, HYBRIDIZATION_FEATURE, RNASE_FEATURE):
            if column in ranked.columns:
                designed[column] = ranked[column].to_numpy()
        safety = tox_details(shortlist)
        jobs.save_results(
            config.job_id,
            {
                "designed_asos.csv": designed,
                "safety_detail.csv": safety,
                "off_targets.csv": off_targets,
            },
        )
        jobs.save_features(config.job_id, ranked)
        jobs.save_layout(config.job_id, layout or _layout_for(config))
        jobs.mark(config.job_id, jobs.DONE)

        send_processing_completed(config.user_email, config.source_info, jobs.public_url(config.job_id))
        logger.info(f"Design job {config.job_id} complete for {config.user_email}.")

    except Exception as e:
        logger.exception(f"Design job failed for {config.user_email}: {e}")
        reason = f"{type(e).__name__}: {e}"
        jobs.mark(config.job_id, jobs.FAILED, error=reason)
        # The submitter has already had the "started" mail, so without this they would wait on a
        # result that is never coming.
        send_processing_failed(config.user_email, config.source_info, reason)


def _report_lost_job(config: JobConfig, future):
    """A job that raises inside execute_tauso_pipeline mails the user from there. This covers the
    other case: the worker process dying, which takes that handler down with it."""
    try:
        future.result()
    except Exception as e:
        logger.exception(f"Design job for {config.user_email} was lost with the worker: {e}")
        reason = f"{type(e).__name__}: {e}"
        jobs.mark(config.job_id, jobs.FAILED, error=reason)
        send_processing_failed(config.user_email, config.source_info, reason)


def trigger_background_job(config: JobConfig):
    """Submit the config to the background pool and return immediately."""
    global _executor
    try:
        future = _get_executor().submit(execute_tauso_pipeline, config)
    except BrokenProcessPool:
        logger.warning("The background pool was broken by an earlier job; starting a new one.")
        _executor.shutdown(wait=False)
        _executor = None
        future = _get_executor().submit(execute_tauso_pipeline, config)
    future.add_done_callback(lambda finished: _report_lost_job(config, finished))
