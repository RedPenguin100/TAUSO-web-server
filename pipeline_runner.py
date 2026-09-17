import logging
import os
import sys
import threading
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool
from dataclasses import dataclass
from typing import Optional

import pandas as pd

import jobs
from email_service import send_processing_completed, send_processing_failed, send_processing_started
from tauso.aso_generation import (
    default_config,
    design_asos,
    summarize_design,
    tox_details,
)
from tauso.data.consts import ASO_SEQUENCE, CANONICAL_GENE_NAME
from tauso.off_target.search import annotate_hits, run_bowtie_search_many

LOG_FORMAT = "%(asctime)s %(levelname)s: %(message)s"


class ColourFormatter(logging.Formatter):
    """Dim the timestamp and colour the level, so an ERROR is findable in a wall of INFO.

    Only installed when the stream is a terminal -- `tty: true` in the compose file is what makes
    that true inside the container. Redirect the logs to a file and this falls back to plain text
    rather than filling it with escape codes.
    """

    DIM = "\033[38;5;245m"
    RESET = "\033[0m"
    LEVEL = {
        logging.DEBUG: "\033[38;5;245m",
        logging.INFO: "\033[36m",
        logging.WARNING: "\033[33m",
        logging.ERROR: "\033[31m",
        logging.CRITICAL: "\033[1;31m",
    }

    def format(self, record):
        colour = self.LEVEL.get(record.levelno, "")
        stamp = f"{self.DIM}{self.formatTime(record, self.datefmt)}{self.RESET}"
        level = f"{colour}{record.levelname}{self.RESET}"
        message = record.getMessage()
        if record.exc_info:
            message = f"{message}\n{self.formatException(record.exc_info)}"
        return f"{stamp} {level}: {message}"


_handler = logging.StreamHandler()
_handler.setFormatter(
    ColourFormatter(datefmt="%H:%M:%S")
    if getattr(_handler.stream, "isatty", lambda: False)()
    else logging.Formatter(LOG_FORMAT, datefmt="%H:%M:%S")
)
logging.basicConfig(level=logging.INFO, handlers=[_handler])
logger = logging.getLogger(__name__)

# The machine budget, declared once in the compose file and derived from here by everything that
# can spend it -- feature computation and the bowtie off-target search alike. Stating it beats
# hardcoding a number per call site: the two were already inconsistent, with featurisation on eight
# workers while bowtie ran single-threaded under a twelve-core limit.
CORES = int(os.environ.get("TAUSO_CORES", "0")) or (os.cpu_count() or 1)
MEMORY_MB = int(os.environ.get("TAUSO_MEMORY_MB", "0"))

# Measured on this image rather than guessed: the parent peaked at ~4.4 GB while parsing the SAM
# from a whole-transcript run, and each joblib/loky worker sat near 130 MB RSS, so 300 MB a worker
# leaves roughly two-fold headroom.
PARENT_MB = 4400
PER_WORKER_MB = 300


def _affordable_cores() -> int:
    """Cores we can actually pay for, which is not always the cores we were given.

    Every worker is a process with its own copy of the working set, so past a point more of them
    only buys an OOM kill. Where the memory budget is stated, it caps the core count.
    """
    cores = max(1, min(CORES, os.cpu_count() or CORES))
    if MEMORY_MB > 0:
        affordable = max(1, (MEMORY_MB - PARENT_MB) // PER_WORKER_MB)
        if affordable < cores:
            logger.info(
                f"Memory budget {MEMORY_MB} MB allows {affordable} workers, not {cores}; using "
                f"{affordable}."
            )
            return affordable
    return cores


WORKERS = _affordable_cores()

# The off-target stage runs after featurisation, on top of everything featurisation still holds --
# the cached genome and transcriptomes, the feature frame, the worker processes -- and it is where
# the container has been pushed past its memory limit. It therefore gets its own, smaller budget
# rather than the machine-wide one. Featurisation is unaffected and keeps WORKERS.
OFFTARGET_WORKERS = max(1, min(int(os.environ.get("TAUSO_OFFTARGET_CORES", "4")), WORKERS))

# An explicit TAUSO_DESIGN_JOBS still wins, for pinning one job's featurisation without touching
# the rest. An earlier note here claimed the ProcessPool worker below was daemonic and so could
# not spawn children; that is no longer true -- on the image's Python 3.12 the executor's workers
# are non-daemonic and nested pools work. Verified in-container 2026-09-07.
DESIGN_JOBS = int(os.environ.get("TAUSO_DESIGN_JOBS", str(WORKERS)))
# TAUSO_FIRST_N bounds how many tiled candidates are featurised. design_asos tiles 5'->3' and
# takes the first N windows, so this covers the 5' end of the target rather than sampling across
# it. Featurising is around 33 seconds of fixed cost plus 0.04 seconds a candidate, so the marginal
# candidate is cheap and the bound can be generous. Set TAUSO_FIRST_N=0 to tile the whole transcript,
# which for a 8.8 kb target is around seven minutes.
FIRST_N = int(os.environ.get("TAUSO_FIRST_N", "500")) or None
# How many candidates make the shortlist: the table, and the bowtie off-target search, which is
# what this bounds. Every candidate scored is still charted and downloadable -- a plot of the
# winners alone reads as a landscape while hiding exactly the stretches a reader most needs to
# see are bad.
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

# Backbones offered as one-click presets for the 20-mer MOE gapmer, which is what most of the
# training data is. A 20-mer has 19 linkages, and the alphabet is LINKAGE_CODES -- '*' for
# phosphorothioate, 'o' for phosphodiester -- so the PO substitutions in the wings are written 'o'.
# The names are positional because these patterns have no accepted names.
BACKBONE_PRESETS = {
    "Full PS": "*" * 19,
    "Var 1": "*ooo***********oo**",
    "Var 2": "*oooo**********oo**",
}

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
# GC across the ASO itself. It sets how tightly the duplex binds, so it moves several of the
# hybridization features with it.
GC_FEATURE = "seq_gc_content"
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

    def offsets(intervals):
        out = []
        for start, end in intervals or []:
            if reverse:
                out.append([locus.gene_end - end, locus.gene_end - start])
            else:
                out.append([start - locus.gene_start, end - locus.gene_start])
        return sorted(out)

    # The UTRs are exonic, so they are subsets of the exon spans rather than a separate track --
    # the chart draws them over the exons to show which part of an exon is untranslated.
    return {
        "length": length,
        # Carried so the results page can sanity-check the target against the cell line without
        # reopening a 1.4 GB GTF: a chrY gene in a female-derived line cannot be expressed.
        "chrom": getattr(locus, "chrom", None),
        "exons": offsets(getattr(locus, "_exon_indices", [])),
        "utr5": offsets(getattr(locus, "_5utr_indices", [])),
        "utr3": offsets(getattr(locus, "_3utr_indices", [])),
    }


SITE_REGION = "site_region"
REPEAT_NOTE = "repeat_note"
COMPLEMENT = str.maketrans("ACGTU", "TGCAA")


def _region_at(position: int, length: int, exons, span) -> str:
    """Where one site sits in the gene model: exonic if it overlaps an exon, intronic if it falls
    between them inside the transcript, unannotated if it is outside the transcript altogether.

    Coarser than tauso's own `target_region`, which separates the UTRs -- the layout carries exon
    bounds only. It is computed per site, so the copies of a repeated candidate are each labelled
    where they actually are rather than inheriting the first occurrence's label.
    """
    if not exons:
        return "unannotated"
    start, end = position, position + length
    if end <= span[0] or start >= span[1]:
        return "unannotated"
    return "exon" if any(start < b and end > a for a, b in exons) else "intron"


def annotate_repeat_sites(designed: pd.DataFrame, pre_mrna: Optional[str], layout: Optional[dict]):
    """Give every copy of a repeated candidate its own position, and say which copy was scored.

    design_asos keys candidates by sequence, so a k-mer occurring several times in the target --
    which a tandem repeat guarantees -- collapses onto its first occurrence: each copy is emitted
    as its own row, but all of them carry that one start. The chart then draws a hole over the
    repeat while the surplus rows stack on a single point. The true positions are recoverable by
    re-scanning the target, so they are handed back here, and every copy after the first says so
    explicitly: it was not scored where it is drawn.
    """
    if designed.empty or not pre_mrna or ASO_SEQUENCE not in designed or "target_start" not in designed:
        return designed

    exons = [tuple(e) for e in (layout or {}).get("exons", [])]
    span = (min(a for a, _ in exons), max(b for _, b in exons)) if exons else (0, 0)

    # An ASO is antisense to its target, so a candidate's sequence is the reverse complement of the
    # window it binds. The index is keyed that way round so a candidate looks itself up directly.
    sequences = designed[ASO_SEQUENCE].astype(str).str.upper()
    target = pre_mrna.upper()
    sites: dict[str, list[int]] = {}
    for width in sorted({len(s) for s in sequences}):
        for i in range(len(target) - width + 1):
            antisense = target[i:i + width].translate(COMPLEMENT)[::-1]
            sites.setdefault(antisense, []).append(i)

    starts, regions, notes = [], [], []
    seen: dict[str, int] = {}
    for sequence, original in zip(sequences, designed["target_start"]):
        found = sites.get(sequence, [])
        index = seen.get(sequence, 0)
        seen[sequence] = index + 1
        position = found[index] if index < len(found) else (found[0] if found else int(original))
        starts.append(position)
        regions.append(_region_at(position, len(sequence), exons, span))
        if len(found) < 2:
            notes.append("")
        elif index == 0:
            others = ", ".join(f"{p:,}" for p in found[1:])
            notes.append(f"first of {len(found)} identical sites (also at {others})")
        else:
            notes.append(
                f"identical to the site at {found[0]:,}, which is where it was scored -- "
                f"this copy carries that score, not one computed here"
            )

    annotated = designed.copy()
    annotated["target_start"] = starts
    annotated[SITE_REGION] = regions
    annotated[REPEAT_NOTE] = notes
    repeats = sum(1 for n in notes if n)
    if repeats:
        logger.info(f"{repeats} candidates share a sequence with another site; positions restored.")
    return annotated


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


def _locus_for(config: JobConfig):
    """This job's target locus, or None for a sequence the user supplied. Fetched once: the layout
    and the pre-mRNA the repeat scan needs both come off it."""
    if config.target_data:
        return None
    try:
        from tauso.populate.calculators.cache import AssetCache

        return AssetCache(genome="GRCh38").get_full_gene_data().get(config.target_mrna_name)
    except Exception:
        logger.warning("Could not read the gene model for %s", config.target_mrna_name)
        return None


def _layout_for(config: JobConfig):
    """The gene model for this job's target, or None for a sequence the user supplied."""
    locus = _locus_for(config)
    return exon_layout(locus) if locus else None


OFFTARGET_COLS = ["rank", "aso_sequence", "off_target_gene", "distance", "region", "chrom", "start", "strand"]


def _offtarget_table(ranked, *, genome, max_distance, exclude_genes=None):
    """One row per hit the shortlist makes to a gene other than its target: the mismatch count,
    the region, and the locus.

    The whole shortlist goes through one bowtie run, so the genome index is loaded once rather
    than once per oligo -- 0.6 s of index against 25 ms of alignment. That batching used to live
    here; tauso's `run_bowtie_search_many` now does it, and streams the alignments off the pipe
    as columns besides, so this calls that instead of keeping a second copy of the parser.
    """
    from tauso.data.data import get_paths

    paths = get_paths(genome)
    sentinel = os.path.join(os.path.dirname(paths["fasta"]), f"{genome}_bowtie_index", "SUCCESS")
    if not os.path.exists(sentinel):
        raise FileNotFoundError(
            f"Bowtie index for {genome} not found; run `tauso setup-bowtie --genome {genome}`."
        )

    exclude = set(exclude_genes or []) | set(pd.Series(ranked[CANONICAL_GENE_NAME]).dropna().unique())
    rank_of = {seq: i + 1 for i, seq in enumerate(ranked[ASO_SEQUENCE].tolist())}
    unique_seqs = list(dict.fromkeys(ranked[ASO_SEQUENCE].tolist()))
    if not unique_seqs:
        return pd.DataFrame(columns=OFFTARGET_COLS)

    hits, _counts = run_bowtie_search_many(
        unique_seqs, genome=genome, max_mismatches=max_distance, threads=OFFTARGET_WORKERS
    )

    annotated = annotate_hits(hits, genome=genome)
    if annotated.empty:
        return pd.DataFrame(columns=OFFTARGET_COLS)
    tbl = annotated[annotated["gene_name"].notna() & ~annotated["gene_name"].isin(exclude)]
    if tbl.empty:
        return pd.DataFrame(columns=OFFTARGET_COLS)

    out = pd.DataFrame(
        {
            "rank": tbl["sequence"].map(rank_of).to_numpy(),
            "aso_sequence": tbl["sequence"].to_numpy(),
            "off_target_gene": tbl["gene_name"].to_numpy(),
            "distance": tbl["mismatches"].to_numpy(),
            "region": tbl["region_type"].to_numpy(),
            "chrom": tbl["chrom"].to_numpy(),
            "start": tbl["start"].to_numpy(),
            "strand": tbl["strand"].to_numpy(),
        }
    )
    return out.sort_values(["rank", "distance", "off_target_gene", "start"]).reset_index(drop=True)


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

    stop_watching = watch_memory()
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
        locus = _locus_for(config)
        layout = exon_layout(locus) if locus else None
        pre_mrna = getattr(locus, "full_mrna", None) if locus else None
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
        memory_note("featurisation and scoring")
        logger.info(f"Scored {len(ranked)} candidate ASOs; building result tables...")

        shortlist = ranked.head(TOP_N)
        off_targets = _offtarget_table(
            shortlist,
            genome="GRCm39" if design_config.organism_name == "mouse" else "GRCh38",
            max_distance=OFFTARGET_MAX_DISTANCE,
            exclude_genes=None,
        )
        memory_note("the off-target table")
        logger.info(f"{len(off_targets)} off-target hits across the top {len(shortlist)}.")

        designed = summarize_design(ranked)
        for column in (ACCESSIBILITY_FEATURE, MFE_FEATURE, HYBRIDIZATION_FEATURE, RNASE_FEATURE,
                       GC_FEATURE):
            if column in ranked.columns:
                designed[column] = ranked[column].to_numpy()
        designed = annotate_repeat_sites(designed, pre_mrna, layout)
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
        jobs.save_layout(config.job_id, layout)
        jobs.mark(config.job_id, jobs.DONE)

        send_processing_completed(config.user_email, config.source_info, jobs.public_url(config.job_id))
        memory_note("saving results")
        logger.info(f"Design job {config.job_id} complete for {config.user_email}.")

    except Exception as e:
        logger.exception(f"Design job failed for {config.user_email}: {e}")
        reason = f"{type(e).__name__}: {e}"
        jobs.mark(config.job_id, jobs.FAILED, error=reason)
        # The submitter has already had the "started" mail, so without this they would wait on a
        # result that is never coming.
        send_processing_failed(config.user_email, config.source_info, reason)
    finally:
        stop_watching()


# Colour only when the logs go to a terminal, which `tty: true` in the compose file arranges.
# Redirected to a file, the escape codes would be noise rather than emphasis.
MEMORY_COLOUR = getattr(sys.stderr, "isatty", lambda: False)()


def _cgroup_mb(name: str):
    """One cgroup memory figure in MB, or None where the file is absent (cgroup v1, no container)."""
    try:
        with open(f"/sys/fs/cgroup/memory.{name}") as handle:
            raw = handle.read().strip()
        return None if raw == "max" else int(raw) / 2**20
    except (OSError, ValueError):
        return None


def _cgroup_stat_mb(field: str):
    """One field of memory.stat in MB, or None. `anon` is what the kernel kills on; `file` is
    page cache, which it reclaims instead. Bowtie maps its index rather than reading it in, so
    memory.current now counts a couple of gigabytes of index as cache: reporting that as though
    it were the job's footprint says 94% of the limit when the job holds less than half of it."""
    try:
        with open("/sys/fs/cgroup/memory.stat") as handle:
            for line in handle:
                key, _, value = line.partition(" ")
                if key == field:
                    return int(value) / 2**20
    except (OSError, ValueError):
        return None
    return None


def _tint(share):
    """Bold green, amber or red by how close a reading came to the limit."""
    if not MEMORY_COLOUR:
        return "", ""
    if share >= 0.85:
        return "\033[1;31m", "\033[0m"
    if share >= 0.6:
        return "\033[1;33m", "\033[0m"
    return "\033[1;32m", "\033[0m"


def watch_memory(interval: float = 4.0, jump_mb: float = 150.0):
    """Sample memory while the job runs, logging only when it moves.

    tauso logs "Starting step: X" but no memory, and those steps are upstream. Sampling on a thread
    and printing only material changes puts a reading between those lines, so growth can be pinned
    to the step that caused it -- without a line every four seconds saying nothing happened.
    """
    stop = threading.Event()

    def sample():
        # Anonymous memory, not memory.current: the mapped bowtie index is page cache, and
        # following the total would report the job climbing when the kernel is only caching a
        # file it can drop again.
        def held():
            return _cgroup_stat_mb("anon") or _cgroup_mb("current")

        last = held() or 0.0
        limit = _cgroup_mb("max")
        while not stop.wait(interval):
            current = held()
            if current is None:
                return
            if abs(current - last) >= jump_mb:
                arrow = "up" if current > last else "down"
                share = current / limit if limit else 0.0
                tint, reset = _tint(share)
                logger.info(f"{tint}MEMORY {arrow} to {current:.0f} MB held"
                            + (f" ({share * 100:.0f}% of limit)" if limit else "") + reset)
                last = current

    threading.Thread(target=sample, daemon=True).start()
    return stop.set


def memory_note(stage: str) -> None:
    """Log what the container is holding, at a point where it matters.

    Read from the cgroup rather than the process: the peak is the sum across the parent, the pool
    worker and every joblib child, which is what the limit is applied to and what the kernel kills
    on. memory.peak is the high-water mark since the container started, so it survives the moment
    that caused it.
    """
    try:
        readings = {}
        for name in ("current", "peak", "max"):
            with open(f"/sys/fs/cgroup/memory.{name}") as handle:
                readings[name] = handle.read().strip()
        current = int(readings["current"]) / 2**20
        peak = int(readings["peak"]) / 2**20
        limit = readings["max"]
        # The kill decision is made on anonymous memory; page cache is reclaimed first. Both are
        # in memory.current, so it is reported as its two parts rather than as one number that
        # reads like the job is at the limit when half of it is a mapped index.
        anon = _cgroup_stat_mb("anon")
        cache = _cgroup_stat_mb("file")
        if limit == "max":
            share, ceiling = 0.0, "no limit"
        else:
            ceiling_mb = int(limit) / 2**20
            share = (anon if anon is not None else peak) / ceiling_mb if ceiling_mb else 0.0
            ceiling = f"{ceiling_mb:.0f} MB limit ({share * 100:.0f}% of it held)"
        # Coloured by how close the held memory came to the limit, so the line that matters is the
        # one that catches the eye: this is the number that decides whether a job survives.
        tint, reset = _tint(share)
        split = "" if anon is None else f" ({anon:.0f} MB held, {cache or 0:.0f} MB reclaimable)"
        logger.info(
            f"{tint}MEMORY after {stage}: {current:.0f} MB now{split}, {peak:.0f} MB peak, {ceiling}{reset}"
        )
    except (OSError, ValueError):
        # cgroup v1, or not in a container at all. Not worth a warning on every job.
        pass


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
