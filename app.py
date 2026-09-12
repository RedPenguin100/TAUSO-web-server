import hashlib
import io
import logging
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit as st
import streamlit.components.v1 as components
from Bio import SeqIO

import jobs
from email_service import send_contact_message

logger = logging.getLogger(__name__)
from pipeline_runner import (
    BACKBONE_PRESETS,
    TOP_N,
    ACCESSIBILITY_FEATURE,
    GC_FEATURE,
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

# The design form reads better in a column; the results page is a wide chart and wants the room.
st.set_page_config(
    page_title="TAUSO | ASO design",
    layout="wide" if st.query_params.get("job") else "centered",
)

_gene_select = components.declare_component(
    "gene_select", path=str(Path(__file__).parent / "components" / "gene_select")
)
_pattern_editor = components.declare_component(
    "pattern_editor", path=str(Path(__file__).parent / "components" / "pattern_editor")
)

# TAUSO also has a mouse genome, but every cell line with expression data here is human.
ORGANISM = "human"

GENE_PLACEHOLDER = "Search a gene"

# Parent-document markup so the note can overlap what follows it. Pure CSS hover: no script, which
# Streamlit would strip anyway.
GENE_HELP_HTML = """
<div class="tauso-help">
  <span class="tauso-help-mark">?</span>
  <div class="tauso-help-box">
    GENCODE v38 (Ensembl 104), GRCh38 primary assembly<br>
    Targeting pre-mRNA transcripts, introns included<br>
    Protein coding and lncRNA genes only
  </div>
</div>
<style>
/* Same height as the input beside it, so the marker centres against the box rather than
   sitting near its top. */
.tauso-help { position:relative; display:flex; align-items:center; height:36px; }
.tauso-help-mark { display:inline-block; width:16px; height:16px; border-radius:50%;
  border:1px solid #D7DDE5; color:#5A6473; background:#FBFCFD; font:10px/14px system-ui;
  text-align:center; cursor:default; user-select:none; }
.tauso-help:hover .tauso-help-mark { border-color:#2A78D6; color:#2A78D6; }
.tauso-help-box { display:none; position:absolute; left:24px; top:50%;
  transform:translateY(-50%); width:268px;
  padding:7px 10px; border:1px solid #D7DDE5; border-radius:8px; background:#fff; color:#141A22;
  font:11px/1.45 system-ui,-apple-system,"Segoe UI",sans-serif;
  box-shadow:0 6px 18px rgba(20,26,34,.16); z-index:9999; }
.tauso-help:hover .tauso-help-box { display:block; }
</style>
"""

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


def _squash(name: str) -> str:
    """A cell line name with punctuation and case dropped, so "SK-N-AS" and "SKNAS" compare equal."""
    return re.sub(r"[^A-Za-z0-9]", "", str(name)).upper()


@st.cache_data(ttl=3600)
def fetch_cell_lines(organism: str):
    """Every cell line this deployment can actually condition on.

    The names come from the DepMap model table tauso ships -- all 2,132 of them -- rather than the
    55-name proxy dict that used to be the only source: that dict was a curated subset, and a line
    outside it could not be offered even with its expression on disk. What is still required is the
    expression itself, so the list is whatever DepMap publishes intersected with what has been
    built here. A line with no expression is left out rather than offered and then failing.
    """
    if organism != "human":
        return []

    data_dir = os.environ.get("TAUSO_DATA_DIR", "/home/mambauser/.tauso_data")
    expression_dir = os.path.join(data_dir, "processed_expression")
    if not os.path.isdir(expression_dir):
        return []
    available = {f.replace("_expression.csv", "") for f in os.listdir(expression_dir)}

    names_by_id = {}
    try:
        import tauso
        from tauso.data.consts import CELL_LINE_TO_DEPMAP_PROXY_DICT, resolve_depmap_id

        models = pd.read_csv(
            os.path.join(os.path.dirname(tauso.__file__), "data", "cell_lines", "depmap_models.csv")
        )
        table_name = {}
        for name, depmap_id in zip(models["cell_line"], models["depmap_id"]):
            if depmap_id in available:
                table_name[depmap_id] = str(name)
                names_by_id.setdefault(depmap_id, []).append(str(name))
        # The curated spellings are nicer to read ("HEK-293" over "HEK293"), so they win where a
        # line has one. Several curated names can point at one id, though, and some of those are
        # typos ("A459" for A549), aliases ("Ts-24" for T24) or experiment labels
        # ("Angptl2/Actin"); taking whichever came last would put those on the menu. The one that
        # matches the table's own name once punctuation and case are dropped is the real spelling,
        # and where none does the shortest wins, which prefers "Hep3B" over "HepG2/Hep3B".
        curated_by_id = {}
        for name in CELL_LINE_TO_DEPMAP_PROXY_DICT:
            depmap_id = resolve_depmap_id(name)
            if depmap_id in available:
                curated_by_id.setdefault(depmap_id, []).append(name)
        for depmap_id, names in curated_by_id.items():
            canonical = _squash(table_name.get(depmap_id, ""))
            best = sorted(names, key=lambda n: (_squash(n) != canonical, len(n), n))[0]
            names_by_id.setdefault(depmap_id, []).insert(0, best)
    except Exception:
        logger.exception("Could not read the DepMap model table; falling back to the proxy dict.")
        from tauso.data.consts import CELL_LINE_TO_DEPMAP_PROXY_DICT, resolve_depmap_id

        for name in CELL_LINE_TO_DEPMAP_PROXY_DICT:
            depmap_id = resolve_depmap_id(name)
            if depmap_id in available:
                names_by_id.setdefault(depmap_id, []).append(name)

    return sorted(names[0] for names in names_by_id.values() if names)


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
        # st.selectbox builds a DOM node per option, so clearing a typed query re-rendered all
        # ~59k of them and stalled for seconds. This component keeps the whole list -- every gene
        # is scrollable and the count is shown -- but only the rows inside the viewport are ever
        # in the DOM, so the cost no longer scales with the size of the genome.
        # Drawn as raw HTML in the parent page, not inside the component and not as a native
        # tooltip: an iframe cannot paint outside its own rectangle, and Streamlit's own tooltip
        # picks its own side. Here the marker sits right of the box, the note opens to its right,
        # and z-index puts it over the chemistry section rather than moving it.
        picker_column, help_column = st.columns([4, 6])
        with picker_column:
            gene = _gene_select(
                genes=genes,
                placeholder=GENE_PLACEHOLDER,
                value=st.session_state.get("gene_select"),
                key="gene_select",
                default=None,
            )
        with help_column:
            st.markdown(GENE_HELP_HTML, unsafe_allow_html=True)
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


def conditions_section(target_name=None):
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
        editing = st.toggle("Edit sugars and linkages")

    preset = CHEMISTRIES[preset_name]
    sugar, backbone = preset["pattern"], preset["ps_pattern"]
    st.caption(f"{len(sugar)} nt · {describe_chemistry(sugar, backbone)}")

    if editing:
        low, high = DESIGN_LENGTH_RANGE
        length = st.slider("Length (nt)", min_value=low, max_value=high, value=len(preset["pattern"]))
        # One-click backbones, for the 20-mer MOE gapmer only: the patterns are 19 linkages long
        # and mean nothing at another length or on another sugar.
        if preset_name == DEFAULT_CHEMISTRY and length == len(preset["pattern"]):
            button_columns = st.columns([1, 1, 1, 3])
            for column, (name, pattern) in zip(button_columns, BACKBONE_PRESETS.items()):
                with column:
                    if st.button(name, key=f"backbone_{name}", width="stretch"):
                        st.session_state["backbone_set"] = pattern
                        st.session_state["backbone_nonce"] = (
                            st.session_state.get("backbone_nonce", 0) + 1
                        )
            with button_columns[3]:
                st.markdown(
                    "",
                    help=(
                        "Backbones seen most often in the training data, so the model is best "
                        "equipped to predict them. `*` is phosphorothioate, `o` is "
                        "phosphodiester.  \n"
                        + "  \n".join(f"**{n}** `{p}`" for n, p in BACKBONE_PRESETS.items())
                    ),
                )

        edited = _pattern_editor(
            length=length,
            code=_sugar_code(preset["pattern"]),
            colour=CHEMISTRY_COLOURS[preset_name],
            label=preset_name,
            sugar=preset["pattern"] if length == len(preset["pattern"]) else None,
            backbone=preset["ps_pattern"] if length == len(preset["pattern"]) else None,
            backbone_set=st.session_state.get("backbone_set"),
            backbone_nonce=st.session_state.get("backbone_nonce", 0),
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
        pairing_notice(target_name, cell_line)

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
# How much of the transcript the chart opens on, and the width at which its marks are drawn at
# full weight. A whole scan is tens of thousands of candidates across 468 px, which is a smear
# rather than a landscape.
WINDOW_NT = 2000

# The name every row of the chart refers to its candidates by.
CANDIDATES = "candidates"

TRACK_SCALE = [(0.0, "#cf4c41"), (0.5, "#e9b23c"), (1.0, "#4aa058")]
# Fixed cut points on the raw score, not percentiles of the run: a percentile always finds a
# "top 1%" even on a gene with no good sites, while an empty top band says the true thing. The
# ramp goes dark to pale so better reads as stronger without needing the legend.
# label, lower bound, upper bound, colour, extra marker radius, opacity. The good bands are drawn
# larger and at full strength and the poor ones recede: on a scan of tens of thousands, the few
# marks worth looking at should not be the same size as the thousands that are not.
SCORE_BANDS = [
    ("Great", "above +20", 20.0, float("inf"), "#0B2E5C", 4.0, 1.00),
    ("Good", "+10 to +20", 10.0, 20.0, "#1F6FD0", 3.0, 1.00),
    ("Average", "0 to +10", 0.0, 10.0, "#7FB0E8", 0.5, 0.85),
    ("Poor", "-15 to 0", -15.0, 0.0, "#B9C4D4", 0.0, 0.65),
    ("No knockdown", "below -15", float("-inf"), -15.0, "#DCE2EA", 0.0, 0.50),
]

GENE_COLOURS = {"exon": "#3D4653", "intron": "#8792A2", "utr5": "#2A78D6", "utr3": "#B4762E",
                # Sequence the canonical transcript does not cover: dashed and faded, so it
                # reads as "not part of this transcript" rather than as intron.
                "outside": "#C6CCD6"}

def _position_figure(designed, score_column, layout=None):
    """Score against transcript position, with the structure and binding of each candidate on their
    own rows beneath, sharing the x axis so a column of marks is one candidate.

    The scores are a WebGL trace and each feature row is a single raster, so panning and zooming
    move work the browser has already done rather than redrawing thousands of shapes. A browser
    without WebGL is handled in the page itself, which downgrades the traces to SVG rather than
    letting Plotly fail with "WebGL is not supported by your browser"."""
    tracks = [
        (ACCESSIBILITY_FEATURE, "open site"),
        (MFE_FEATURE, "MFE"),
        (HYBRIDIZATION_FEATURE, "binding dG"),
        (RNASE_FEATURE, "RNase H1"),
        (GC_FEATURE, "GC"),
    ]
    tracks = [(column, label) for column, label in tracks if column in designed.columns]
    rows_at = []

    data = designed.sort_values("target_start")
    x = data["target_start"].to_numpy()
    scores = data[score_column].to_numpy(dtype=float)
    low, high = float(x.min()), float(x.max())

    # Score, then the transcript, then a row a feature, then the whole scan again small enough to
    # navigate by.
    heights = [0.60, 0.05] + [0.345 / len(tracks)] * len(tracks) + [0.005]
    overview_row = 3 + len(tracks)
    figure = make_subplots(rows=overview_row, cols=1, shared_xaxes=True,
                           vertical_spacing=0.006, row_heights=heights)

    # A full window is thousands of candidates over the panel's width, where solid marks read as a
    # block of colour; a short scan is a few dozen, where they read as scattered specks.
    crowding = min(1.0, float((x < x.min() + WINDOW_NT).sum()) / WINDOW_NT)

    if layout:
        for begin, finish in layout.get("exons", []):
            if finish > low and begin < high:
                figure.add_vrect(x0=max(begin, low), x1=min(finish, high), row=1, col=1,
                                 fillcolor="#5A6473", opacity=0.09, line_width=0, layer="below")
        # Solid only between the first and last exon, which is the canonical transcript. Outside
        # that the axis is still gene -- the span belongs to some other isoform -- and drawing it
        # in the intron colour said "intron" about sequence that is not in this transcript at all.
        # HBB is the plain case: 59% of its axis sits upstream of the canonical transcript, and
        # read as an intron before the 5'UTR.
        exon_spans = layout.get("exons") or []
        body_start = min((begin for begin, _ in exon_spans), default=low)
        body_end = max((finish for _, finish in exon_spans), default=high)
        figure.add_shape(type="line", x0=max(body_start, low), x1=min(body_end, high),
                         y0=0.5, y1=0.5, row=2, col=1,
                         layer="below", line=dict(color=GENE_COLOURS["intron"], width=2))
        for outside_start, outside_end in ((low, body_start), (body_end, high)):
            if outside_end > outside_start:
                figure.add_shape(type="line", x0=max(outside_start, low), x1=min(outside_end, high),
                                 y0=0.5, y1=0.5, row=2, col=1, layer="below",
                                 line=dict(color=GENE_COLOURS["outside"], width=1.5, dash="dash"))
        for begin, finish in layout.get("exons", []):
            if finish > low and begin < high:
                figure.add_shape(type="rect", x0=max(begin, low), x1=min(finish, high),
                                 y0=0.12, y1=0.88, row=2, col=1, layer="above",
                                 fillcolor=GENE_COLOURS["exon"], line_width=0)
        # Drawn over the exons, and shorter, so an exon reads as coding where no UTR covers it.
        # Older jobs were saved before the layout carried UTRs, hence the default.
        for key in ("utr5", "utr3"):
            for begin, finish in layout.get(key, []):
                if finish > low and begin < high:
                    figure.add_shape(type="rect", x0=max(begin, low), x1=min(finish, high),
                                     y0=0.26, y1=0.74, row=2, col=1, layer="above",
                                     fillcolor=GENE_COLOURS[key], line_width=0)

    figure.update_yaxes(range=[0, 1], showticklabels=False, ticks="", showgrid=False,
                        zeroline=False, row=2, col=1)
    gene_domain = figure.layout.yaxis2.domain
    # "canonical" heads the key because every span below it comes from the canonical transcript,
    # while the axis spans the whole gene -- sequence outside that transcript is neither exon nor
    # intron here, and saying so stops the unshaded stretches being read as intronic.
    figure.add_annotation(
        text=(f"<b>canonical</b><br>"
              f"<span style='color:{GENE_COLOURS['exon']}'>\u2588</span> exon"
              f" &nbsp;<span style='color:{GENE_COLOURS['intron']}'><b>\u25ac</b></span> intron<br>"
              f"<span style='color:{GENE_COLOURS['utr5']}'>\u2588</span> 5'UTR"
              f" &nbsp;<span style='color:{GENE_COLOURS['utr3']}'>\u2588</span> 3'UTR<br>"
              f"<span style='color:{GENE_COLOURS['outside']}'>- -</span> other isoform(s)"),
        # Nudged up: the key grew to three lines, and centring it on the track leaves the last
        # line hanging below the gene bar. This sits the block level with the track it describes.
        xref="paper", yref="paper", x=1.01, y=sum(gene_domain) / 2 + 0.035,
        xanchor="left", yanchor="middle", showarrow=False, align="left",
        font=dict(size=12, color="#3D4653"),
    )

    # Colour limits come from the whole scan, so a shade means the same thing at any zoom, and
    # each row is one image however many candidates it covers.
    track_values = []
    for i, (column, label) in enumerate(tracks, start=3):
        values = data[column].to_numpy(dtype=float)
        figure.add_trace(
            go.Heatmap(
                x=x, z=[values], colorscale=TRACK_SCALE,
                zmin=float(np.nanmin(values)), zmax=float(np.nanmax(values)),
                hoverinfo="skip",
                name=label, showlegend=False,
                colorbar=dict(orientation="h", thickness=13, len=0.14,
                              x=1.20, xanchor="left", yanchor="middle",
                              tickfont=dict(size=8), tickangle=0, outlinewidth=0,
                              ticklabelposition="outside bottom", tickmode="array",
                              tickvals=[float(np.nanmin(values)), float(np.nanmax(values))],
                              ticktext=[f"{np.nanmin(values):.3g}", f"{np.nanmax(values):.3g}"]),
            ),
            row=i, col=1,
        )
        # The row's values are carried by the score trace instead of an invisible marker trace of
        # their own. Five such traces meant every candidate was uploaded and rendered six times
        # over, for points nobody can see.
        track_values.append(values)
        figure.update_yaxes(showticklabels=False, ticks="", showgrid=False, zeroline=False,
                            row=i, col=1)
        middle = sum(figure.layout[f"yaxis{i}"].domain) / 2
        rows_at.append({"column": len(track_values) - 1, "y": middle,
                        "label": label, "note": len(figure.layout.annotations)})
        figure.data[-1].colorbar.y = middle - 0.012
        figure.add_annotation(text=label, xref="paper", yref="paper", x=1.01, y=middle,
                              xanchor="left", yanchor="middle", showarrow=False,
                              font=dict(size=11, color="#3D4653"))

    # The whole scan drawn inside the range slider itself, so the window sits on the profile it
    # navigates rather than in an empty box beneath it. A line, and not WebGL: the slider renders
    # SVG traces only, and one path is cheaper than thousands of points anyway. The row that owns
    # it is a sliver, because the slider is the visible copy.
    figure.add_trace(
        go.Scatter(
            x=x, y=scores, mode="lines", hoverinfo="skip", showlegend=False,
            line=dict(color="#8792A2", width=1),
        ),
        row=overview_row, col=1,
    )
    figure.update_yaxes(visible=False, row=overview_row, col=1)
    # The slider draws the traces of the axis it belongs to, so the profile has to live on a real
    # row. That row is a sliver and is painted over: the slider below is the copy meant to be seen.
    sliver = figure.layout[f"yaxis{overview_row}"].domain
    figure.add_shape(type="rect", xref="paper", yref="paper", x0=0, x1=1,
                     y0=sliver[0] - 0.01, y1=sliver[1] + 0.01,
                     fillcolor="white", line_width=0, layer="above")
    figure.update_xaxes(
        rangeslider=dict(visible=True, thickness=0.13, bgcolor="#FBFCFD",
                         bordercolor="#D7DBE0", borderwidth=1,
                         yaxis=dict(rangemode="auto")),
        row=overview_row, col=1,
    )

    margin = 0.03 * (float(np.nanmax(scores)) - float(np.nanmin(scores)) or 1)
    # Nothing zooms vertically: the score axis is pinned to the whole scan, and a drag that also
    # moved it would put the same score at a different height.
    figure.update_yaxes(fixedrange=True)
    # Under the axis title, because "score" is the first thing anyone wants explained and the
    # answer is not obvious. A plotly annotation rather than anything drawn in the page: it hovers
    # inside the plot, where there is room, and moves with the axis it belongs to.
    figure.add_annotation(
        text="<b>?</b>", xref="paper", yref="paper", x=-0.055, y=0.905,
        xanchor="center", yanchor="middle", showarrow=False,
        font=dict(size=11, color="#8792A2"),
        hovertext=(
            "<b>Why a score, and not % inhibition?</b><br>"
            "The model is fitted to knockdown de-meaned within each experiment rather than to the "
            "absolute percentage.<br>Removing the offset between experiments lets it learn from all "
            "of them at once,<br>so the ranking is sharper and the biology in the features is kept."
        ),
        hoverlabel=dict(bgcolor="white", bordercolor="#D7DBE0",
                        font=dict(size=11, color="#3D4653")),
        captureevents=True,
    )
    figure.update_yaxes(title_text="score", gridcolor="#EEF0F3",
                        range=[float(np.nanmin(scores)) - margin, float(np.nanmax(scores)) + margin],
                        row=1, col=1)
    figure.update_xaxes(showgrid=False, ticks="", showline=False, zeroline=False,
                        showspikes=True, spikemode="across", spikesnap="cursor",
                        spikethickness=1, spikedash="dot", spikecolor="#3D4653",
                        range=[low, min(low + WINDOW_NT, high)])
    figure.update_xaxes(showticklabels=False, row=overview_row, col=1)
    figure.update_xaxes(showticklabels=True, side="bottom", ticks="outside", ticklen=3,
                        title_text="position in the transcript (nt)",
                        title_font=dict(size=11), title_standoff=2, tickfont=dict(size=10),
                        row=1, col=1)
    # Room under the score panel for that scale to sit in.
    bottom, top = figure.layout.yaxis.domain
    figure.layout.yaxis.domain = (bottom + 0.085, top)
    figure.update_layout(
        # "closest", not "x": with a trace a band and no distance limit, x mode reports the
        # nearest point in every band at once, so one candidate under the cursor produced five
        # scores in the tooltip. The crosshair writes the feature rows itself, so nothing is lost
        # by asking plotly for the single point actually being pointed at.
        width=1040, height=600, dragmode="pan", hovermode="closest", plot_bgcolor="white",
        paper_bgcolor="white", margin=dict(l=56, r=280, t=46, b=16),
        showlegend=True,
        legend=dict(orientation="h", yanchor="bottom", y=1.004, xanchor="left", x=0,
                    font=dict(size=10), itemsizing="constant", bgcolor="rgba(0,0,0,0)",
                    itemwidth=30, tracegroupgap=0),
        hoversubplots="axis", hoverdistance=-1, spikedistance=-1,
        hoverlabel=dict(bgcolor="white", bordercolor="#D7DBE0",
                        font=dict(size=11, color="#3D4653")),
        transition=dict(duration=0),
    )
    # One row a candidate, one column a feature row, in the order rows_at records. A plain nested
    # list, not an array: plotly serialises numpy as base64 typed data, which the crosshair cannot
    # index. Every band carries the same columns, so the crosshair can read whichever band the
    # cursor happens to be over.
    carried = ([[None if v != v else float(v) for v in column] for column in zip(*track_values)]
               if track_values else None)

    # A trace a band rather than one for everything: plotly's legend then toggles a band on and
    # off in the browser, with no round trip. Empty bands are added too -- "above +20: 0" is the
    # useful answer for a gene with no standout sites, and a missing entry would not say it.
    for adjective, span, lo, hi, colour, bonus, alpha in SCORE_BANDS:
        inside = np.flatnonzero((scores > lo) & (scores <= hi))
        # Where to go when this band is asked for, worked out here rather than in the browser:
        # the position of its best candidate, and the width of the whole scan to size the window
        # against. Reading it back out of the plotted arrays meant trusting how plotly stores
        # them, which is what the previous attempt got wrong.
        focus = float(x[inside][int(np.argmax(scores[inside]))]) if len(inside) else None
        figure.add_trace(
            go.Scattergl(
                x=x[inside], y=scores[inside], mode="markers",
                meta=dict(focus=focus, span=float(high - low)),
                marker=dict(size=6 - 2 * crowding + bonus, color=colour,
                            opacity=alpha - 0.15 * crowding),
                customdata=[carried[i] for i in inside] if carried else None,
                # The adjective sits on its own line above the range it stands for, so the legend
                # reads as a judgement with its evidence under it rather than as a bare number.
                name=f"<b>{adjective}</b><br>{span} ({len(inside):,})",
                legendgroup=adjective, showlegend=True,
                hovertemplate="%{x:,.0f} nt<br>score %{y:.2f}<extra></extra>",
            ),
            row=1, col=1,
        )

    # The cut points drawn where the eye can use them, faint enough not to compete with the marks.
    for _adjective, _span, lo, _hi, colour, _bonus, _alpha in SCORE_BANDS:
        if lo not in (float("-inf"), float("inf")) and low is not None:
            figure.add_hline(y=lo, row=1, col=1, line=dict(color=colour, width=1, dash="dot"),
                             opacity=0.55)

    return figure, rows_at


CHART_HTML = """
<style>
  /* The y axes are fixed, so plotly puts an east-west resize cursor over the entire plot. It
     reads as "drag to resize", which is not what dragging does here, so the ordinary pointer is
     restored -- panning still works, it just stops advertising itself as a resize handle. */
  .js-plotly-plot .nsewdrag,
  .js-plotly-plot .ewdrag,
  .js-plotly-plot .nsdrag,
  .js-plotly-plot .drag { cursor: default !important; }
</style>
<div id="chart"></div>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<script>
  const figure = __FIGURE__;
  const rows = __ROWS__;
  const gd = document.getElementById("chart");
  const shapes = figure.layout.shapes || [];
  const notes = figure.layout.annotations || [];

  let showing = null;
  // The crosshair repaints by relayouting the figure. A relayout that lands mid-drag interrupts
  // plotly's pan, and the matched axes can come out of it at different ranges -- which shows up
  // as one pane sliding on its own. Hovering is therefore ignored while a drag is in flight.
  let dragging = false;

  // Plotly's WebGL traces do not degrade: on a browser without WebGL they render nothing and
  // print "WebGL is not supported by your browser". Downgrading them to SVG here keeps the chart
  // working there while leaving every other browser on the fast path with all its points.
  function hasWebGL() {
    try {
      const probe = document.createElement("canvas");
      return !!(window.WebGLRenderingContext &&
                (probe.getContext("webgl") || probe.getContext("experimental-webgl")));
    } catch (e) {
      return false;
    }
  }
  if (!hasWebGL()) {
    figure.data.forEach(function (trace) {
      if (trace.type === "scattergl") { trace.type = "scatter"; }
    });
  }

  Plotly.newPlot(gd, figure.data, figure.layout, __CONFIG__).then(function () {
    // Double-clicking a legend entry isolates that band, which plotly does on its own. It is
    // only half the question though -- "which are the good ones" is usually followed by "where
    // are they" -- so the view is moved to cover them as well. The x axes are matched to the
    // rangeslider's axis, so the range has to be set on that one for the rest to follow.
    const master = Object.keys(figure.layout).filter(function (k) {
      return k.indexOf("xaxis") === 0 && figure.layout[k].rangeslider &&
             figure.layout[k].rangeslider.visible;
    })[0] || "xaxis";

    gd.on("plotly_legenddoubleclick", function (event) {
      const trace = gd.data[event.curveNumber];
      // The band's best candidate and the width of the scan, both computed server-side and
      // carried on the trace, so this does not depend on how plotly stores the plotted arrays.
      const focus = trace && trace.meta ? trace.meta.focus : null;
      if (focus === null || focus === undefined) { return true; }
      const centre = focus;
      const full = (trace.meta.span || 0);
      const half = Math.max(600, full * 0.05) / 2;
      // Every x axis, not only the master: matched axes should follow, and setting them all
      // removes any doubt about which one the group is actually driven by.
      const change = {};
      Object.keys(gd.layout).forEach(function (k) {
        if (k.indexOf("xaxis") === 0) { change[k + ".range"] = [centre - half, centre + half]; }
      });
      // After plotly's own isolate, so the two do not fight over the same relayout.
      setTimeout(function () { Plotly.relayout(gd, change); }, 120);
      return true;
    });

    gd.on("plotly_relayouting", function () { dragging = true; });
    gd.on("plotly_relayout", function () { dragging = false; });
    gd.addEventListener("mousedown", function () { dragging = true; });
    window.addEventListener("mouseup", function () { dragging = false; });

    gd.on("plotly_hover", function (event) {
      if (dragging) { return; }
      const point = event.points[0];
      const at = point.pointIndex;
      const curve = point.curveNumber;
      if (at === undefined || at === showing) { return; }
      showing = at;
      const line = {
        type: "line", xref: "x", yref: "paper", x0: point.x, x1: point.x, y0: 0, y1: 1,
        line: { color: "#3D4653", width: 1, dash: "dot" }, layer: "above",
      };
      const labelled = notes.map(function (note) { return Object.assign({}, note); });
      rows.forEach(function (row) {
        // Whichever band trace the cursor is over: every band carries the same columns, so the
        // hovered curve is always a valid place to read the feature values from.
        const cd = gd.data[curve] && gd.data[curve].customdata;
        const carried = cd ? cd[at] : null;
        const value = carried ? carried[row.column] : undefined;
        const shown = (value === null || value === undefined)
          ? "\u2014" : Number(value).toPrecision(3);
        labelled[row.note].text = row.label + "   <b>" + shown + "</b>";
        labelled[row.note].opacity = 1;
        labelled[row.note].bgcolor = "rgba(0,0,0,0)";
        labelled[row.note].font = { size: 11, color: "#3D4653" };
      });
      Plotly.relayout(gd, { shapes: shapes.concat([line]), annotations: labelled });
    });
    gd.on("plotly_unhover", function () {
      showing = null;
      Plotly.relayout(gd, { shapes: shapes, annotations: notes });
    });
  });
</script>
"""


def _chart_html(figure, rows_at) -> str:
    """The figure with a crosshair of its own: one line the height of the whole chart, and the
    value of every feature written on its row."""
    config = {"scrollZoom": True, "displaylogo": False, "doubleClick": "reset",
              "modeBarButtonsToRemove": ["select2d", "lasso2d", "autoScale2d"]}
    return (
        CHART_HTML
        .replace("__FIGURE__", figure.to_json())
        .replace("__ROWS__", json.dumps(rows_at))
        .replace("__CONFIG__", json.dumps(config))
    )


# Superscript digits, so a gene carries its mismatch count without a column of its own:
# HBD\u2070 is a perfect match, BBS9\u00b2 is two mismatches away.
SUPERSCRIPT = {0: "\u2070", 1: "\u00b9", 2: "\u00b2", 3: "\u00b3", 4: "\u2074"}


def _mark(gene, distance):
    return f"{gene}{SUPERSCRIPT.get(int(distance), '')}"


# The model has 552 features, most of them members of a family: 161 rbp_*, 80 ohe_*, 40 ribo_*.
# TreeSHAP splits credit across correlated features, so no single one looks important even when
# its family dominates -- the totals are what carry meaning, and what a reader can act on.
def _diverging(value, limit):
    """Blue for a contribution that lifted the score, red for one that pushed it down.

    Computed here rather than through pandas' background_gradient, which needs matplotlib -- a
    50 MB dependency to colour a table is a poor trade.
    """
    if value is None or value != value or not limit:
        return ""
    share = max(-1.0, min(1.0, float(value) / limit))
    if share >= 0:
        red, green, blue = 255 - int(95 * share), 255 - int(50 * share), 255
    else:
        red, green, blue = 255, 255 - int(75 * -share), 255 - int(80 * -share)
    return f"background-color: rgb({red},{green},{blue})"


SHAP_COLUMNS = 9

SHAP_FAMILIES = {
    "structure": "position", "ohe": "motif", "hybr": "duplex", "seq": "composition",
    # fold_* is the folding energy of the site and access_* how unpaired it is: two measurements
    # of the same thing, so they are one family rather than two that always move together.
    "fold": "accessibility", "access": "accessibility", "off": "off-target",
    "expr": "expression", "rbp": "RBP binding", "ribo": "ribosome", "selfaso": "self-structure",
    "mod": "chemistry", "rnase": "RNase H1", "tox": "toxicity motifs", "flank": "flanks",
    "cai": "codon usage", "enc": "codon usage", "tai": "codon usage", "on": "on-target sites",
    "halflife": "half-life", "struct": "position", "sense": "accessibility",
    "transfection": "assay", "volume": "assay", "density": "assay", "chem": "chemistry",
    "interaction": "self-structure",
}


@st.cache_data(ttl=3600, show_spinner=False)
def shap_by_family(job_id: str, sequences: tuple):
    """Per candidate, what pushed its score up or down, totalled by feature family.

    Exact TreeSHAP from the booster itself rather than the sampling approximation: xgboost can do
    it natively, and 100 candidates cost about half a second. Read from the saved features, so
    this costs the design run nothing and works on jobs that finished before it existed.
    """
    path = jobs.features_path(job_id)
    # Asked of tauso rather than hardcoded: "v1" has pointed at two different files now, and a
    # stale path here would explain a model the run never used.
    try:
        from tauso.inference.predict import DEFAULT_VERSION, MODEL_FILES

        filename = MODEL_FILES[DEFAULT_VERSION]["filename"]
    except Exception:
        filename = "tauso_score_v1.json"
    model = os.path.join(
        os.environ.get("TAUSO_DATA_DIR", "/home/mambauser/.tauso_data"), "models", filename,
    )
    if not path.exists() or not os.path.exists(model):
        return None
    try:
        import xgboost as xgb

        booster = xgb.Booster()
        booster.load_model(model)
        names = booster.feature_names
        frame = pd.read_parquet(path)
        wanted = frame[frame["aso_sequence"].isin(sequences)]
        if wanted.empty:
            return None
        # A job scored by an older feature set cannot be explained by this booster: reindex would
        # quietly fill the missing columns with NaN and produce confident nonsense.
        covered = sum(1 for n in names if n in wanted.columns)
        if covered < 0.9 * len(names):
            logger.info("SHAP skipped for %s: %d of %d features present", job_id, covered, len(names))
            return None
        matrix = wanted.reindex(columns=names)
        contributions = booster.predict(
            xgb.DMatrix(matrix, feature_names=names), pred_contribs=True
        )[:, :-1]
        families = [SHAP_FAMILIES.get(n.split("_")[0], n.split("_")[0]) for n in names]
        totals = pd.DataFrame(contributions, columns=families).T.groupby(level=0).sum().T
        totals.index = wanted["aso_sequence"].to_numpy()
        return totals
    except Exception:
        # An explanation is a nicety; never let it take the results page down with it.
        logger.exception("SHAP unavailable for %s", job_id)
        return None


def _offtarget_labels(off_targets):
    """Per candidate, the genes it hits: the two worst named with their mismatch count, rest counted.

    Ordered by severity rather than alphabetically -- fewest mismatches first, then most hits --
    so the name that shows is the one worth reacting to. A perfect match to a paralog is the point
    of this column; "3" never said that.
    """
    labels = {}
    if off_targets.empty:
        return labels
    for sequence, group in off_targets.groupby("aso_sequence"):
        order = (
            group.groupby("off_target_gene")
            .agg(best=("distance", "min"), hits=("distance", "size"))
            .sort_values(["best", "hits"], ascending=[True, False])
        )
        if order.empty:
            continue
        shown = [_mark(gene, row.best) for gene, row in order.head(2).iterrows()]
        rest = len(order) - len(shown)
        labels[sequence] = ", ".join(shown) + (
            f", and {rest} other{'s' if rest != 1 else ''}" if rest > 0 else ""
        )
    return labels


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


@st.cache_data(ttl=3600)
def _result_tables(job_id: str):
    """The three tables a finished job wrote. Finished results never change, so a rerun reads them
    from the cache rather than from disk."""
    return tuple(
        pd.read_csv(jobs.results_path(job_id, name))
        for name in ("designed_asos.csv", "safety_detail.csv", "off_targets.csv")
    )


@st.cache_data(ttl=3600)
def gene_index():
    """gene_name -> [chromosome, gene span in nt].

    Derived from the GTF once and then kept as a small JSON beside it. Scanning 1.4 GB takes a
    couple of seconds, which is fine as a one-off and far too slow to sit between choosing a gene
    and seeing what it implies. The span is what sets the candidate count: tiling runs across the
    whole pre-mRNA, so a 20-mer yields exactly span - 19 candidates.
    """
    import re

    data_dir = os.environ.get("TAUSO_DATA_DIR", "/home/mambauser/.tauso_data")
    cached = os.path.join(data_dir, "gene_index.json")
    if os.path.exists(cached):
        try:
            with open(cached) as handle:
                return json.load(handle)
        except Exception:
            pass

    gtf = os.path.join(data_dir, "GRCh38.gtf")
    table = {}
    if not os.path.exists(gtf):
        return table
    with open(gtf) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            name = re.search(r'gene_name "([^"]+)"', fields[8])
            if name:
                table[name.group(1)] = [fields[0], int(fields[4]) - int(fields[3]) + 1]
    try:
        with open(cached, "w") as handle:
            json.dump(table, handle)
    except OSError:
        # A read-only data directory just means paying the scan once per process.
        pass
    return table


def gene_chromosomes():
    """gene_name -> chromosome."""
    return {gene: entry[0] for gene, entry in gene_index().items()}


# Runtime is close to linear in the length of the target, and the rate is a property of the
# machine this happens to run on -- cores, memory, disk. It is therefore read from the data
# directory rather than committed: another host will have another number, and nobody should have
# to edit the source to correct it. See RUNTIME_CALIBRATION_FILE for how it was measured.
RUNTIME_CALIBRATION_FILE = "runtime_calibration.json"
RUNTIME_DEFAULTS = {"seconds_per_nt": 0.023, "fixed_seconds": 47.0}


@st.cache_data(ttl=300)
def runtime_calibration():
    """Seconds per nucleotide and fixed overhead for this deployment."""
    path = os.path.join(
        os.environ.get("TAUSO_DATA_DIR", "/home/mambauser/.tauso_data"),
        RUNTIME_CALIBRATION_FILE,
    )
    values = dict(RUNTIME_DEFAULTS)
    try:
        with open(path) as handle:
            stored = json.load(handle)
        for key in RUNTIME_DEFAULTS:
            if isinstance(stored.get(key), (int, float)):
                values[key] = float(stored[key])
    except Exception:
        pass
    return values


# A ceiling on what one submission may cost. Jobs run one at a time, so a very long target does
# not merely inconvenience whoever asked for it -- it blocks the queue for everyone else. Set
# TAUSO_MAX_RUNTIME_MINUTES to 0 to lift the limit entirely.
MAX_RUNTIME_MINUTES = int(os.environ.get("TAUSO_MAX_RUNTIME_MINUTES", "60"))


def runtime_estimate(gene):
    """(seconds, spoken duration) for designing against `gene`, or None if its span is unknown.

    Driven by the target's length rather than an exact candidate count: the count depends on the
    oligo length too, and quoting a precise figure would suggest a precision this does not have.
    """
    entry = gene_index().get(gene)
    if not entry or len(entry) < 2:
        return None
    calibration = runtime_calibration()
    seconds = calibration["fixed_seconds"] + calibration["seconds_per_nt"] * entry[1]
    if seconds < 90:
        spoken = f"about {round(seconds / 10) * 10:.0f} seconds"
    elif seconds < 3600:
        spoken = f"about {seconds / 60:.0f} minutes"
    else:
        spoken = f"about {seconds / 3600:.1f} hours"
    return seconds, spoken


@st.cache_data(ttl=3600)
def cell_line_sex(name: str):
    """Donor sex for a cell line, from DepMap's model table, or None if it cannot be resolved."""
    if not name:
        return None
    try:
        from tauso.data.consts import resolve_depmap_id

        depmap_id = resolve_depmap_id(name)
        if not depmap_id:
            return None
        path = os.path.join(
            os.environ.get("TAUSO_DATA_DIR", "/home/mambauser/.tauso_data"), "Model.csv"
        )
        models = pd.read_csv(path, usecols=["ModelID", "Sex"])
        row = models[models["ModelID"] == depmap_id]
        return None if row.empty else str(row.iloc[0]["Sex"])
    except Exception:
        return None


@st.cache_data(ttl=3600)
def target_expression(gene: str, cell_line: str):
    """The target's TPM in that cell line, read from the cohort table, or None if unavailable.

    Taken from processed_expression rather than a finished job's features, because this has to
    answer while the job is still running -- which is the only time the answer is any use.
    """
    if not gene or not cell_line:
        return None
    try:
        from tauso.data.consts import resolve_depmap_id

        depmap_id = resolve_depmap_id(cell_line)
        if not depmap_id:
            return None
        path = os.path.join(
            os.environ.get("TAUSO_DATA_DIR", "/home/mambauser/.tauso_data"),
            "processed_expression",
            f"{depmap_id}_expression.csv",
        )
        if not os.path.exists(path):
            return None
        table = pd.read_csv(path, usecols=["Gene", "expression_TPM"])
        row = table[table["Gene"] == gene]
        return None if row.empty else float(row.iloc[0]["expression_TPM"])
    except Exception:
        return None


# "Not detected" has to mean what it says. A cutoff at 1 TPM called 0.34 undetected and 0.058 as
# well; both are low expression, not absence. At 1e-5 the notice fires only on a reading that is
# zero, so the claim is exactly true. Anything above it is a judgement for whoever is designing.
NOT_DETECTED_TPM = 1e-5
# Named in the notice so the claim is attributable: both the expression and the donor sex come
# from this release, and a reader who doubts either knows exactly what to go and check.
EXPRESSION_SOURCE = "DepMap Public 25Q3"


def pairing_notice(gene, cell_line):
    """One line under the cell line, when the target cannot be present in the line chosen.

    Two claims, weighted differently on purpose. Expression is a measurement and measurements can
    be wrong for a given batch, so a zero is amber. A chrY gene in a female-derived line is not a
    measurement -- the sequence is absent from the genome -- so that one is red.
    """
    if not gene or not cell_line:
        return
    chrom = gene_chromosomes().get(gene)
    sex = cell_line_sex(cell_line)
    if chrom in ("chrY", "Y") and sex and sex.lower().startswith("f"):
        _notice(
            "#B42318",
            f"Warning: {gene} is on chromosome Y and {cell_line} is a female cell line "
            f"({EXPRESSION_SOURCE}) \u2014 the target is not in these cells.",
        )
        return
    tpm = target_expression(gene, cell_line)
    if tpm is not None and tpm < NOT_DETECTED_TPM:
        # Scientific notation below zero-proper: "0.000 TPM" would read as an absence that the
        # number does not actually claim.
        reading = "0 TPM" if tpm == 0 else f"{tpm:.1e} TPM"
        _notice(
            "#B54708",
            f"Warning: {gene} is not detected in {cell_line} "
            f"({reading}, {EXPRESSION_SOURCE}).",
        )


def _notice(colour, text):
    st.markdown(
        f"<div style='margin-top:-6px;font-size:12px;line-height:1.4;color:{colour}'>{text}</div>",
        unsafe_allow_html=True,
    )


def results_page(job_id: str):
    """Everything one finished job produced, opened from its own address."""
    job = jobs.get(job_id)
    if job is None:
        st.error(f"No job called {job_id}.")
        st.markdown("[Design something new](/)")
        return

    st.markdown(
        "<style>.block-container{max-width:1200px;}</style>", unsafe_allow_html=True
    )
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

    designed, safety, off_targets = _result_tables(job_id)
    # Named, not positional: the explanatory feature columns are appended after it.
    score_column = next(c for c in designed.columns if c.startswith("tauso_score_"))

    chemistry = describe_chemistry(
        parameters.get("chemical_pattern", ""), parameters.get("ps_pattern", "")
    )
    st.caption(chemistry)

    shortlist = designed[designed["aso_sequence"].isin(safety["aso_sequence"])]

    top = st.columns(4)
    top[0].metric("Candidates", len(designed))
    top[1].metric("Length", f"{len(parameters.get('chemical_pattern', ''))} nt")
    top[2].metric("Cell line", parameters.get("cell_line") or "none")
    top[3].metric("Off-target hits", len(off_targets))

    st.subheader("Score along the transcript")
    st.caption(
        "Each point is one candidate, placed where it binds. Higher is better predicted knockdown "
        "relative to the others here — it ranks candidates, it is not a percent. The tracks beneath "
        "carry the same candidates, each shaded red at the low end of its own range and green at "
        "the high end — so a red **binding dG** mark is the most negative free energy, the "
        "tightest duplex. The strip at the foot is the whole transcript: drag the shaded box on it "
        "to move the view, or its edges to widen it."
    )
    layout = jobs.get_layout(job_id)
    figure, rows_at = _position_figure(designed, score_column, layout)
    components.html(_chart_html(figure, rows_at), height=640, scrolling=False)
    if layout:
        exonic = sum(b - a for a, b in layout["exons"])
        st.caption(
            f"{job['target']} is {layout['length']:,} nt with {len(layout['exons'])} exons, "
            f"{exonic:,} nt exonic ({100 * exonic / layout['length']:.0f}%). Shaded stretches are "
            f"exonic; the gene track at the foot shows the part of the transcript drawn here."
        )

    starts = shortlist.head(10)["target_start"].sort_values().tolist()
    if len(starts) > 1 and starts[-1] - starts[0] < 2 * len(parameters.get("chemical_pattern", "x" * 20)):
        st.warning(
            f"The top 10 all start between {starts[0]} and {starts[-1]}. Tiling moves one nucleotide "
            "at a time, so these overlap heavily — they are one site rather than ten choices."
        )

    st.subheader(f"Top {len(shortlist)} candidates")
    st.caption(
        f"The chart above carries all {len(designed):,} scored candidates; this is the shortlist, "
        "which is also what the off-target search covers."
    )
    if not off_targets.empty:
        per_gene = (
            off_targets.groupby("off_target_gene")
            .agg(best=("distance", "min"), hits=("distance", "size"))
            .sort_values(["best", "hits"], ascending=[True, False])
        )
        named = ", ".join(
            f"**{_mark(gene, row.best)}** ({int(row.hits)})"
            for gene, row in per_gene.head(4).iterrows()
        )
        regions = off_targets["region"].value_counts()
        coding = int(regions.get("CDS", 0))
        exonic = int(regions.get("exon", 0))
        st.caption(
            f"Off-targets across the top {len(shortlist)}: {named}"
            + (f", and {len(per_gene) - 4} more genes" if len(per_gene) > 4 else "")
            + f". {coding} hit CDS, {exonic} hit exons, {len(off_targets) - coding - exonic} "
            "fall in introns."
        )

    merged = shortlist.merge(safety, on="aso_sequence", how="left")
    offtarget_genes = merged["aso_sequence"].map(_offtarget_labels(off_targets)).fillna("")
    families = shap_by_family(job_id, tuple(merged["aso_sequence"]))
    hits = off_targets.groupby("aso_sequence")["distance"].value_counts().unstack(fill_value=0)
    accessibility = merged.get(ACCESSIBILITY_FEATURE)
    binding = merged.get(HYBRIDIZATION_FEATURE)
    # Only worth a column when something in the shortlist actually repeats; a run with no repeated
    # candidate leaves it None and dropna below takes the column out.
    repeat_note = merged["repeat_note"] if "repeat_note" in merged else None
    if repeat_note is not None and not (repeat_note.fillna("").astype(str).str.strip() != "").any():
        repeat_note = None
    exact_match = merged["aso_sequence"].map(hits.get(0, {})).fillna(0).astype(int)
    one_mismatch = merged["aso_sequence"].map(hits.get(1, {})).fillna(0).astype(int)
    two_mismatch = merged["aso_sequence"].map(hits.get(2, {})).fillna(0).astype(int)
    table = pd.DataFrame(
        {
            "#": merged["rank"],
            "sequence (5'->3')": merged["aso_sequence"],
            "start": merged["target_start"],
            "region": merged["site_region"] if "site_region" in merged else None,
            "repeat": repeat_note,
            "score": merged[score_column].round(2),
            "liabilities": merged.apply(_liability_chips, axis=1),
            # A perfect match to another gene is the off-target that matters most, and it was
            # computed and written to the CSV without ever being shown.
            "0mm": exact_match,
            "1mm": one_mismatch,
            "2mm": two_mismatch,
            "off-target genes": offtarget_genes,
        }
    ).dropna(axis=1, how="all")
    # A perfect match to another gene disqualifies a candidate, so the cell is filled rather than
    # left as a number among numbers -- it should be findable while scrolling, not read.
    shown = table
    if "0mm" in table.columns:
        shown = table.style.map(
            lambda hits: "background-color:#F7D4D7; color:#7A1620; font-weight:600"
            if hits
            else "",
            subset=["0mm"],
        )
    st.dataframe(
        shown,
        hide_index=True,
        width="stretch",
    )
    st.caption(
        "**0mm**, **1mm** and **2mm** count genomic hits to a gene other than the target, at that "
        "many mismatches -- a 0mm hit is a perfect match elsewhere, and the one worth acting on. "
        "**off-target genes** names the two worst of them, fewest mismatches first, with the "
        "mismatch count as a superscript -- HBD\u2070 is a perfect match, BBS9\u00b2 is two away. "
        "The full list, one row per hit with its locus and region, is in the downloaded "
        "off_targets.csv. "
        "The biophysical columns (accessibility, duplex energy, MFE, RNase H1) are in the "
        "downloaded designed_asos.csv. "
        "The breakdown below says what moved each score. "
        "**region** is where that site sits in the gene model. **repeat** appears when the same "
        "sequence occurs at more than one site in the target: each copy is drawn at its own "
        "position, but only the first was scored, and the others carry that score."
    )

    if families is not None and not families.empty:
        st.subheader("SHAP breakdown")
        st.caption(
            "What moved each candidate's score, by feature family, straight from the model "
            "(exact TreeSHAP). A row sums to the distance from the model's average, so the "
            "columns explain the ranking this model produced -- not why an oligo works in a cell."
        )
        # The families are ranked by how much they actually move scores in this run, and the tail
        # is summed into one column: nineteen columns is wider than the screen, and the last ten
        # carry almost nothing.
        weight = families.abs().mean().sort_values(ascending=False)
        shown = list(weight.head(SHAP_COLUMNS).index)
        rest = [c for c in families.columns if c not in shown]
        breakdown = families[shown].round(2)
        if rest:
            breakdown["other"] = families[rest].sum(axis=1).round(2)
        breakdown.insert(0, "sequence (5'->3')", breakdown.index)
        breakdown.insert(0, "#", merged["rank"].to_numpy()[: len(breakdown)])
        numeric = shown + (["other"] if rest else [])
        limit = float(np.nanmax(np.abs(breakdown[numeric].to_numpy()))) or 1.0
        st.dataframe(
            breakdown.style
            .map(lambda v: _diverging(v, limit), subset=numeric)
            .format("{:+.2f}", subset=numeric),
            hide_index=True,
            width="stretch",
        )
        if rest:
            st.caption("**other** sums " + ", ".join(sorted(rest)) + ".")

    st.subheader("Downloads")
    for name in jobs.RESULT_FILES:
        path = jobs.results_path(job_id, name)
        st.download_button(name, path.read_bytes(), file_name=f"{job['target']}_{name}", mime="text/csv")

    with st.expander("What this run was"):
        st.json({"job": job_id, "target": job["target"], **parameters})


def contact_page():
    """A way to ask for something the server currently refuses to do."""
    st.markdown("<style>.block-container{max-width:760px;}</style>", unsafe_allow_html=True)
    st.title("TAUSO")
    st.subheader("Contact us")
    st.write(
        "If you have a business or academic request, wish to calculate ASOs for a long gene, or "
        "any other matter, please fill out your details and your request below."
    )

    # Kept in session state so the confirmation survives the rerun the button causes, and the
    # form is not left looking as though nothing happened.
    if st.session_state.get("contact_sent"):
        st.success("Thank you \u2014 your message has been sent. We will be in touch by email.")
        st.markdown("[Back to designing](/)")
        return

    with st.container(border=True):
        sender = st.text_input("Your email")
        message = st.text_area(
            "Your request", height=200, placeholder="Please describe your request here"
        )
        st.caption("We reply to the address above. Nothing else is collected.")
        submitted = st.button("Send", type="primary", width="stretch")

    if submitted:
        if "@" not in (sender or ""):
            st.error("Please leave an email address so we can reply.")
        elif not message.strip():
            st.error("Please describe your request before sending.")
        elif send_contact_message(sender.strip(), message.strip()):
            st.session_state["contact_sent"] = True
            st.rerun()
        else:
            st.error(
                "This deployment has no contact address configured, so the message was not sent. "
                "Nothing was lost \u2014 please copy your text before leaving."
            )

    st.markdown("[Back to designing](/)")


def main():
    _clear_interrupted_jobs()
    opened = st.query_params.get("job")
    if opened:
        results_page(opened)
        return
    if st.query_params.get("contact"):
        contact_page()
        return

    st.title("TAUSO")
    st.caption("Design antisense oligonucleotides against a human transcript.")

    target_name, target_sequence, source_info = target_section()

    with st.container(border=True):
        sugar, backbone, transfection, dosage, density, cell_line = conditions_section(target_name)

    email = st.text_input("Email for results")

    # Said before the button rather than after: the run is emailed, so this is the last moment the
    # number is any use in deciding whether to press it.
    estimate = runtime_estimate(target_name) if target_name and not target_sequence else None
    too_long = False
    if estimate:
        seconds, spoken = estimate
        too_long = MAX_RUNTIME_MINUTES > 0 and seconds > MAX_RUNTIME_MINUTES * 60
        if too_long:
            st.error(
                f"Scoring ASOs on the {target_name} pre-mRNA would take **{spoken}**, over the "
                f"{MAX_RUNTIME_MINUTES}-minute limit. Please [contact us](?contact=1) or bear "
                "with us until we speed up our process."
            )
        else:
            st.caption(
                f"Tiling the whole pre-mRNA of {target_name} takes **{spoken}**. "
                "You will be emailed a link when it finishes."
            )

    if not st.button("Design ASOs", type="primary", width="stretch", disabled=too_long):
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
