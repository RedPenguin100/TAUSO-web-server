import hashlib
import io
import json
import os

import streamlit as st
from Bio import SeqIO

from pipeline_runner import (
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

# TAUSO also has a mouse genome, but every cell line with expression data here is human.
ORGANISM = "human"

# Label for the no-cell-line option, set apart from the real cell lines around it.
NO_CELL_LINE = "— None (no cell metadata) —"

GENE_PLACEHOLDER = "Gene name, e.g. MALAT1"


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


@st.cache_data(ttl=3600)
def fetch_gene_index(name_upper: str):
    """The gene with this name, matched without regard to case, or None."""
    return {gene.upper(): gene for gene in fetch_genes()}.get(name_upper)


@st.cache_data(ttl=3600)
def suggest_genes(query: str, limit: int = 8):
    """Gene names to offer when the typed one is not a name itself: those starting with it first,
    then those merely containing it."""
    query = query.upper()
    genes = fetch_genes()
    starts = [g for g in genes if g.upper().startswith(query)]
    if len(starts) >= limit:
        return starts[:limit]
    contains = [g for g in genes if query in g.upper() and g not in starts]
    return (starts + contains)[:limit]


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

        # A text box rather than a list of every gene: the reference holds tens of thousands, and
        # putting them all in a select control makes the browser render the lot whenever the filter
        # is cleared.
        query = st.text_input("Gene", placeholder=GENE_PLACEHOLDER, label_visibility="collapsed")
        query = query.strip()
        if not query:
            return None, None, None

        gene = fetch_gene_index(query.upper())
        if gene is None:
            matches = suggest_genes(query)
            if not matches:
                st.warning(f"No gene named {query!r}. Names come from the GRCh38 annotation.")
                return None, None, None
            gene = st.radio(f"Did you mean:", matches, horizontal=True)
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
    """Chemistry and assay conditions. Every one of these is a model input, so the defaults are
    stated rather than hidden: the sugar/backbone pair defines the oligo, and transfection, dosage
    and cell density describe the experiment the prediction is conditioned on."""
    preset_column, transfection_column = st.columns([2, 2])
    with preset_column:
        preset_name = st.selectbox("Chemistry", list(CHEMISTRIES), index=list(CHEMISTRIES).index(DEFAULT_CHEMISTRY))
    with transfection_column:
        transfection = st.selectbox("Transfection", TRANSFECTION_METHODS)

    preset = CHEMISTRIES[preset_name]
    sugar_column, backbone_column = st.columns([2, 2])
    with sugar_column:
        sugar = st.text_input("Sugar", value=preset["pattern"], key=f"sugar_{preset_name}")
    with backbone_column:
        backbone = st.text_input("Backbone", value=preset["ps_pattern"], key=f"backbone_{preset_name}")

    dosage_column, density_column = st.columns([2, 2])
    with dosage_column:
        dosage = st.number_input(
            "Dosage (nM)", min_value=DOSAGE_RANGE_NM[0], max_value=DOSAGE_RANGE_NM[1],
            value=DEFAULT_DOSAGE_NM, step=100,
        )
    with density_column:
        density = st.number_input(
            "Cells per well", min_value=CELL_DENSITY_RANGE[0], max_value=CELL_DENSITY_RANGE[1],
            value=DEFAULT_CELL_DENSITY, step=1000,
        )

    st.caption(
        f"{len(sugar)}-mer · M 2'-MOE · C cEt · L LNA · d DNA · ∗ phosphorothioate. "
        "Defaults are the median of the training experiments."
    )
    return sugar, backbone, transfection, int(dosage), int(density)


def main():
    st.title("TAUSO")
    st.caption("Design antisense oligonucleotides against a human transcript.")

    target_name, target_sequence, source_info = target_section()

    with st.container(border=True):
        sugar, backbone, transfection, dosage, density = conditions_section()

        cell_line = st.selectbox(
            "Cell line (optional)", [NO_CELL_LINE] + fetch_cell_lines(ORGANISM)
        )
        if cell_line == NO_CELL_LINE:
            st.markdown(
                ":orange[**No cell metadata:** some features revert to NaN or to a default. "
                "For detail, see TAUSO Supplementary Material S1.]"
            )

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
    trigger_background_job(
        JobConfig(
            target_mrna_name=target_name,
            target_data=target_sequence,
            source_info=source_info,
            user_email=email,
            cell_line="None" if cell_line == NO_CELL_LINE else cell_line,
            chemical_pattern=sugar,
            ps_pattern=backbone,
            transfection=transfection,
            dosage_nm=dosage,
            cell_density=density,
        )
    )
    st.success(f"Queued. Results for {source_info} will be emailed to {email}.")
    st.caption("A run takes a few minutes. You will get an email either way, including if it fails.")


if __name__ == "__main__":
    main()
