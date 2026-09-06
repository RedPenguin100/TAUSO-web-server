import hashlib
import io
import json
import os
from pathlib import Path

import altair as alt
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from Bio import SeqIO

import jobs
from pipeline_runner import (
    ACCESSIBILITY_FEATURE,
    HYBRIDIZATION_FEATURE,
    MFE_FEATURE,
    RNASE_FEATURE,
    describe_chemistry,
    DESIGN_LENGTH_RANGE,
    CELL_DENSITY_RANGE,
    CHEMISTRIES,
    DEFAULT_CELL_DENSITY,
    DEFAULT_CHEMISTRY,
    DEFAULT_DOSAGE_NM,
    DOSAGE_RANGE_NM,
    MAX_TARGET_LENGTH,
    TRANSFECTION_METHODS,
    JobConfig,
    describe_pattern_problem,
    trigger_background_job,
)

st.set_page_config(page_title="TAUSO | ASO design", layout="centered")

_pattern_editor = components.declare_component(
    "pattern_editor", path=str(Path(__file__).parent / "components" / "pattern_editor")
)

# TAUSO also has a mouse genome, but every cell line with expression data here is human.
ORGANISM = "human"

GENE_PLACEHOLDER = "Search a gene"

# One colour per modified sugar, matching the circles the editor draws.
CHEMISTRY_COLOURS = {"2'-MOE": "#2A78D6", "cEt": "#EB6834"}


def _sugar_code(pattern: str) -> str:
    """The single modified-sugar code in a gapmer pattern."""
    return next((c for c in pattern if c != "d"), "M")


@st.cache_resource
def _clear_interrupted_jobs():
    """Once per server process: nothing in flight survived the restart that got us here."""
    return jobs.fail_interrupted()


@st.cache_data(ttl=3600)
def fetch_genes():
    db_path = os.path.join(os.environ.get("TAUSO_DATA_DIR", "/home/mambauser/.tauso_data"), "available_genes.json")
    if not os.path.exists(db_path):
        return []

    with open(db_path, "r") as f:
        return json.load(f)


@st.cache_data(ttl=3600)
def fetch_cell_lines(organism: str):
    """Cell lines of `organism` this deployment can actually condition on: the expression files
    present in the data directory, named so that design_asos resolves them. A DepMap id whose
    expression was never downloaded, or whose name TAUSO cannot resolve, is left out rather than
    offered and ignored."""
    from tauso.data.consts import CELL_LINE_TO_DEPMAP, CELL_LINE_TO_DEPMAP_PROXY_DICT, resolve_depmap_id

    if organism != "human":
        return []

    expression_dir = os.path.join(
        os.environ.get("TAUSO_DATA_DIR", "/home/mambauser/.tauso_data"), "processed_expression"
    )
    if not os.path.isdir(expression_dir):
        return []
    available = {f.replace("_expression.csv", "") for f in os.listdir(expression_dir)}

    # Several dataset spellings map to one DepMap id; collect them so each line is offered once.
    names_by_id = {}
    for name in CELL_LINE_TO_DEPMAP_PROXY_DICT:
        depmap_id = resolve_depmap_id(name)
        if depmap_id in available:
            names_by_id.setdefault(depmap_id, []).append(name)

    canonical = {v: k for k, v in CELL_LINE_TO_DEPMAP.items()}
    chosen = []
    for depmap_id, names in names_by_id.items():
        preferred = canonical.get(depmap_id)
        # The canonical spelling is the one to show when it resolves; punctuation differences do
        # not matter to the lookup, but some canonical names have no entry of their own.
        chosen.append(preferred if preferred and resolve_depmap_id(preferred) == depmap_id else min(names, key=len))
    return sorted(chosen)


def parse_fasta_input(raw_text: str):
    """Parse pasted or uploaded FASTA into (name, sequence). The name carries a hash of the
    sequence so two different sequences under one header stay distinguishable downstream."""
    raw_text = raw_text.strip()

    if not raw_text.startswith(">"):
        cleaned_seq = "".join(raw_text.split())
        base_name = "Custom_Sequence"
    else:
        with io.StringIO(raw_text) as string_stream:
            try:
                record = next(SeqIO.parse(string_stream, "fasta"))
                base_name = record.id
                cleaned_seq = str(record.seq)
            except StopIteration:
                raise ValueError("The provided FASTA format is invalid or empty.")

    cleaned_seq = cleaned_seq.upper().replace("T", "U")
    seq_hash = hashlib.md5(cleaned_seq.encode("utf-8")).hexdigest()[:8]
    return f"{base_name}_{seq_hash}", cleaned_seq


def target_section():
    """The target to design against: a gene from the reference, or a sequence the user supplies.
    Returns (name, sequence, description); sequence is empty for a gene, which the worker looks up."""
    source = st.radio(
        "Target", ["Gene", "FASTA"], horizontal=True, label_visibility="collapsed"
    )

    if source == "Gene":
        genes = fetch_genes()
        if not genes:
            st.error("The gene database is not initialised yet.")
            return None, None, None
        gene = st.selectbox(
            "Gene", genes, index=None, placeholder=GENE_PLACEHOLDER, label_visibility="collapsed"
        )
        if gene is None:
            return None, None, None
        return gene, "", f"Selected Gene: {gene}"

    pasted = st.text_area(
        "Paste FASTA or a raw sequence",
        height=140,
        placeholder=">MyTranscript\nAUGCGUACGUUAG…",
    )
    uploaded = st.file_uploader("…or upload a file", type=["fasta", "fa", "txt"])

    # A file is the more deliberate action of the two, so it wins if both are present.
    if uploaded is not None:
        name, sequence = parse_fasta_input(uploaded.getvalue().decode("utf-8"))
        return name, sequence, f"Uploaded File: {uploaded.name}"
    if pasted.strip():
        name, sequence = parse_fasta_input(pasted.strip())
        return name, sequence, "Pasted Sequence"
    return None, None, None


def conditions_section():
    """The oligo, then the assay. Every one of these is a model input, so the defaults are stated
    rather than hidden: the sugar/backbone pair defines the oligo, while transfection, dosage, cell
    density and cell line describe the experiment the prediction is conditioned on."""
    chemistry_column, edit_column = st.columns([3, 2])
    with chemistry_column:
        preset_name = st.segmented_control(
            "Chemistry", list(CHEMISTRIES), default=DEFAULT_CHEMISTRY
        ) or DEFAULT_CHEMISTRY
    with edit_column:
        st.markdown("<div style='height:1.75rem'></div>", unsafe_allow_html=True)
        editing = st.toggle("Edit sugars and backbone")

    preset = CHEMISTRIES[preset_name]
    sugar, backbone = preset["pattern"], preset["ps_pattern"]
    st.caption(f"{len(sugar)} nt · {describe_chemistry(sugar, backbone)}")

    if editing:
        low, high = DESIGN_LENGTH_RANGE
        length = st.slider("Length (nt)", min_value=low, max_value=high, value=len(preset["pattern"]))
        edited = _pattern_editor(
            length=length,
            code=_sugar_code(preset["pattern"]),
            colour=CHEMISTRY_COLOURS[preset_name],
            label=preset_name,
            sugar=preset["pattern"] if length == len(preset["pattern"]) else None,
            backbone=preset["ps_pattern"] if length == len(preset["pattern"]) else None,
            key="pattern_editor",
            default={"sugar": preset["pattern"], "backbone": preset["ps_pattern"]},
        )
        sugar = (edited or {}).get("sugar") or sugar
        backbone = (edited or {}).get("backbone") or backbone
        st.caption(
            "Click a sugar to swap it with DNA, or drag across several. The row above is the "
            "backbone: filled is phosphorothioate, hollow is phosphodiester. "
            "\\* One modified chemistry per oligo — mixmers are not supported yet."
        )

    # These four describe the experiment the prediction is conditioned on rather than the oligo.
    # Each may be left blank, which reaches the model as missing rather than as a made-up value.
    st.markdown("**Experimental conditions** &nbsp;·&nbsp; *optional*", unsafe_allow_html=True)
    transfection_column, cell_line_column = st.columns(2)
    with transfection_column:
        transfection = st.selectbox(
            "Transfection", TRANSFECTION_METHODS, index=None, placeholder="Not specified"
        )
    with cell_line_column:
        cell_line = st.selectbox(
            "Cell line", fetch_cell_lines(ORGANISM), index=None, placeholder="Not specified"
        )

    dosage_column, density_column = st.columns(2)
    with dosage_column:
        dosage = st.number_input(
            "Dosage (nM)", min_value=DOSAGE_RANGE_NM[0], max_value=DOSAGE_RANGE_NM[1],
            value=None, step=100, placeholder="Not specified",
        )
    with density_column:
        density = st.number_input(
            "Cells per well", min_value=CELL_DENSITY_RANGE[0], max_value=CELL_DENSITY_RANGE[1],
            value=None, step=1000, placeholder="Not specified",
        )

    st.caption(
        "The model predicts the best sequences from the information it has. "
        "The more you give it, the better the sequences."
    )

    return sugar, backbone, transfection, dosage, density, cell_line



# One scale for every track: red at the low end of the values, green at the high end.
TRACK_RANGE = ["#cf4c41", "#e9b23c", "#4aa058"]


def _track(data, field, label, title):
    """One row of circles under the plot, a candidate each, coloured by `field`: red at the low end
    of the values, green at the high end."""
    colours = TRACK_RANGE
    return (
        alt.Chart(data.assign(track=label))
        .mark_circle(size=110)
        .encode(
            x=alt.X("target_start:Q", axis=None, scale=alt.Scale(zero=False, nice=False)),
            y=alt.Y("track:N", axis=None),
            color=alt.Color(
                f"{field}:Q",
                scale=alt.Scale(range=colours),
                legend=alt.Legend(
                    title=title, orient="left", direction="horizontal",
                    gradientLength=90, gradientThickness=8, titleLimit=150,
                    labelFontSize=9, titleFontSize=10, format=".3~g",
                ),
            ),
            tooltip=[
                alt.Tooltip("rank:Q"),
                alt.Tooltip("target_start:Q", title="start"),
                alt.Tooltip(f"{field}:Q", title=title, format=".2f"),
            ],
        )
        .properties(height=44)
    )


def _position_chart(designed, score_column):
    """Score against transcript position, with the structure and binding of each candidate on their
    own rows beneath, sharing the x scale so a column of marks is one candidate."""
    data = designed.copy()
    scatter = (
        alt.Chart(data)
        .mark_circle(size=70, opacity=0.85, color="#2A78D6")
        .encode(
            x=alt.X("target_start:Q", title=None, scale=alt.Scale(zero=False, nice=False)),
            # The scores of a shortlist sit close together, so the axis follows them rather than
            # reaching down to zero and flattening the differences.
            y=alt.Y(f"{score_column}:Q", title="score", scale=alt.Scale(zero=False, nice=True)),
            tooltip=[
                alt.Tooltip("rank:Q"),
                alt.Tooltip("aso_sequence:N", title="sequence"),
                alt.Tooltip("target_start:Q", title="start"),
                alt.Tooltip(f"{score_column}:Q", title="score", format=".2f"),
            ],
        )
        .properties(height=230)
    )

    # Accessibility and folding energy describe the same thing, so they are pushed together into
    # one block and the other tracks keep their spacing.
    structure = [
        _track(data, field, label, label)
        for field, label in ((ACCESSIBILITY_FEATURE, "open site"), (MFE_FEATURE, "MFE"))
        if field in data
    ]
    rows = [scatter]
    if structure:
        rows.append(alt.vconcat(*structure, spacing=0) if len(structure) > 1 else structure[0])
    if HYBRIDIZATION_FEATURE in data:
        # Free energy in kcal/mol, as measured: the more negative, the tighter, and red marks the
        # tight end of the scale.
        rows.append(_track(data, HYBRIDIZATION_FEATURE, "binding", "binding dG (kcal/mol)"))
    if RNASE_FEATURE in data:
        rows.append(_track(data, RNASE_FEATURE, "RNase H1", "RNase H1 motif fit"))

    return alt.vconcat(*rows, spacing=10).resolve_scale(x="shared", color="independent")


def _liability_chips(row):
    """The flags worth scrutinising on one candidate, as short labels."""
    chips = []
    if row.get("tox_cpg_count", 0) > 0:
        chips.append(f"CpG x{int(row['tox_cpg_count'])}")
    if abs(row.get("tox_g4hunter_max", 0) or 0) >= 1.5 or (row.get("tox_grun_count", 0) or 0) > 0:
        chips.append(f"G4 {row.get('tox_g4hunter_max', 0):.1f}")
    if (row.get("offtarget_rrna", 0) or 0) > 0:
        chips.append("rRNA")
    return ", ".join(chips) if chips else "-"


def results_page(job_id: str):
    """Everything one finished job produced, opened from its own address."""
    job = jobs.get(job_id)
    if job is None:
        st.error(f"No job called {job_id}.")
        st.page_link("app.py", label="Design something new")
        return

    st.title("TAUSO")
    parameters = job["parameters"]

    if job["status"] in (jobs.QUEUED, jobs.RUNNING):
        st.info("This design is still running. The page will show the results when it finishes.")
        if st.button("Check again"):
            st.rerun()
        return
    if job["status"] == jobs.FAILED:
        st.error("This design did not finish.")
        st.code(job["error"] or "no reason recorded", language="text")
        return
    if not jobs.has_results(job_id):
        st.error("This job is marked finished but its result tables are missing.")
        return

    designed = pd.read_csv(jobs.results_path(job_id, "designed_asos.csv"))
    safety = pd.read_csv(jobs.results_path(job_id, "safety_detail.csv"))
    off_targets = pd.read_csv(jobs.results_path(job_id, "off_targets.csv"))
    # Named, not positional: the explanatory feature columns are appended after it.
    score_column = next(c for c in designed.columns if c.startswith("tauso_score_"))

    chemistry = describe_chemistry(
        parameters.get("chemical_pattern", ""), parameters.get("ps_pattern", "")
    )
    st.caption(chemistry)

    top = st.columns(4)
    top[0].metric("Candidates", len(designed))
    top[1].metric("Length", f"{len(parameters.get('chemical_pattern', ''))} nt")
    top[2].metric("Cell line", parameters.get("cell_line") or "none")
    top[3].metric("Off-target hits", len(off_targets))

    st.subheader("Score along the transcript")
    st.caption(
        "Each point is one candidate, placed where it binds. Higher is better predicted knockdown "
        "relative to the others here — it ranks candidates, it is not a percent. The tracks beneath "
        "carry the same candidates. Red is the low end of each scale and green the high end, so a "
        "red **binding** mark is the most negative free energy, meaning the tightest duplex."
    )
    st.altair_chart(_position_chart(designed, score_column), use_container_width=True)

    starts = designed.head(10)["target_start"].sort_values().tolist()
    if len(starts) > 1 and starts[-1] - starts[0] < 2 * len(parameters.get("chemical_pattern", "x" * 20)):
        st.warning(
            f"The top 10 all start between {starts[0]} and {starts[-1]}. Tiling moves one nucleotide "
            "at a time, so these overlap heavily — they are one site rather than ten choices."
        )

    st.subheader("Candidates")
    merged = designed.merge(safety, on="aso_sequence", how="left")
    hits = off_targets.groupby("aso_sequence")["distance"].value_counts().unstack(fill_value=0)
    accessibility = merged.get(ACCESSIBILITY_FEATURE)
    binding = merged.get(HYBRIDIZATION_FEATURE)
    one_mismatch = merged["aso_sequence"].map(hits.get(1, {})).fillna(0).astype(int)
    two_mismatch = merged["aso_sequence"].map(hits.get(2, {})).fillna(0).astype(int)
    table = pd.DataFrame(
        {
            "#": merged["rank"],
            "sequence (5'->3')": merged["aso_sequence"],
            "start": merged["target_start"],
            "score": merged[score_column].round(2),
            "open": accessibility.round(2) if accessibility is not None else None,
            "binding": binding.round(1) if binding is not None else None,
            "RNase H1": merged[RNASE_FEATURE].round(2) if RNASE_FEATURE in merged else None,
            "MFE": merged[MFE_FEATURE].round(3) if MFE_FEATURE in merged else None,
            "liabilities": merged.apply(_liability_chips, axis=1),
            "1mm": one_mismatch,
            "2mm": two_mismatch,
        }
    ).dropna(axis=1, how="all")
    st.dataframe(
        table,
        hide_index=True,
        use_container_width=True,
    )
    st.caption(
        "**open** is how unpaired the target site is over a 60-nt window; **binding** is the "
        "DNA:RNA duplex free energy in kcal/mol, more negative being a tighter duplex; "
        "**MFE** is the folding energy of the site itself, more negative being more structured; **RNase H1** is how well the local dinucleotide context suits the enzyme "
        "that cuts. "
        "**1mm** and **2mm** count genomic hits to a gene other than the target."
    )

    st.subheader("Downloads")
    for name in jobs.RESULT_FILES:
        path = jobs.results_path(job_id, name)
        st.download_button(name, path.read_bytes(), file_name=f"{job['target']}_{name}", mime="text/csv")

    with st.expander("What this run was"):
        st.json({"job": job_id, "target": job["target"], **parameters})


def main():
    _clear_interrupted_jobs()
    opened = st.query_params.get("job")
    if opened:
        results_page(opened)
        return

    st.title("TAUSO")
    st.caption("Design antisense oligonucleotides against a human transcript.")

    target_name, target_sequence, source_info = target_section()

    with st.container(border=True):
        sugar, backbone, transfection, dosage, density, cell_line = conditions_section()

    email = st.text_input("Email for results", placeholder="you@lab.org")

    if not st.button("Design ASOs", type="primary", use_container_width=True):
        return

    if target_name is None:
        st.error("Choose a gene, or paste or upload a sequence.")
        return
    if not email or "@" not in email:
        st.error("Enter an email address — results are delivered by email.")
        return
    if target_sequence and len(target_sequence) > MAX_TARGET_LENGTH:
        st.error(
            f"That sequence is {len(target_sequence):,} nt. The limit is {MAX_TARGET_LENGTH:,} — "
            "longer than any human transcript, so this is usually a paste that went wrong."
        )
        return
    problem = describe_pattern_problem(sugar, backbone)
    if problem:
        st.error(problem)
        return

    # "None" is TAUSO's no-cell-line sentinel and is passed through as that string: the half-life,
    # codon-usage and off-target features each handle it explicitly, while a Python None would leave
    # design_asos on its own default cell line.
    parameters = {
        "cell_line": cell_line,
        "chemical_pattern": sugar,
        "ps_pattern": backbone,
        "transfection": transfection,
        "dosage_nm": dosage,
        "cell_density": density,
    }
    job_id = jobs.create(target_name, source_info, email, parameters)
    trigger_background_job(
        JobConfig(
            target_mrna_name=target_name,
            target_data=target_sequence,
            source_info=source_info,
            user_email=email,
            job_id=job_id,
            cell_line="None" if cell_line is None else cell_line,
            chemical_pattern=sugar,
            ps_pattern=backbone,
            transfection=transfection,
            dosage_nm=dosage,
            cell_density=density,
        )
    )
    st.success("Queued. This page is where the results will appear.")
    st.markdown(f"**Your results:** {jobs.public_url(job_id)}")
    st.caption(
        f"A run takes a few minutes. Keep the link — it is emailed to {email} when the design "
        "finishes, and again if it fails."
    )


if __name__ == "__main__":
    main()
