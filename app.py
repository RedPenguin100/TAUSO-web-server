import hashlib
import io
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit as st
import streamlit.components.v1 as components
from Bio import SeqIO

import jobs
from pipeline_runner import (
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
# How much of the transcript the chart opens on, and the width at which its marks are drawn at
# full weight. A whole scan is tens of thousands of candidates across 468 px, which is a smear
# rather than a landscape.
WINDOW_NT = 2000

# The name every row of the chart refers to its candidates by.
CANDIDATES = "candidates"

TRACK_SCALE = [(0.0, "#cf4c41"), (0.5, "#e9b23c"), (1.0, "#4aa058")]
GENE_COLOURS = {"exon": "#3D4653", "intron": "#8792A2", "utr5": "#2A78D6", "utr3": "#B4762E"}

def _position_figure(designed, score_column, layout=None):
    """Score against transcript position, with the structure and binding of each candidate on their
    own rows beneath, sharing the x axis so a column of marks is one candidate.

    The scores are a WebGL trace and each feature row is a single raster, so panning and zooming
    move work the browser has already done rather than redrawing thousands of shapes."""
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
    figure.add_trace(
        go.Scattergl(
            x=x, y=scores, mode="markers",
            marker=dict(size=6 - 2 * crowding, color="#2A78D6", opacity=0.85 - 0.45 * crowding),
            name="score", showlegend=False,
            hovertemplate="%{x:,.0f} nt<br>score %{y:.2f}<extra></extra>",
        ),
        row=1, col=1,
    )

    if layout:
        for begin, finish in layout.get("exons", []):
            if finish > low and begin < high:
                figure.add_vrect(x0=max(begin, low), x1=min(finish, high), row=1, col=1,
                                 fillcolor="#5A6473", opacity=0.09, line_width=0, layer="below")
        figure.add_shape(type="line", x0=low, x1=high, y0=0.5, y1=0.5, row=2, col=1,
                         layer="below", line=dict(color=GENE_COLOURS["intron"], width=2))
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
              f" &nbsp;<span style='color:{GENE_COLOURS['utr3']}'>\u2588</span> 3'UTR"),
        # Nudged up: the key grew to three lines, and centring it on the track leaves the last
        # line hanging below the gene bar. This sits the block level with the track it describes.
        xref="paper", yref="paper", x=1.01, y=sum(gene_domain) / 2 + 0.035,
        xanchor="left", yanchor="middle", showarrow=False, align="left",
        font=dict(size=12, color="#3D4653"),
    )

    # Colour limits come from the whole scan, so a shade means the same thing at any zoom, and
    # each row is one image however many candidates it covers.
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
        # An invisible point a candidate: the strip is a raster and has no value to hover, and
        # this is what puts a number on this row at the position under the cursor.
        figure.add_trace(
            go.Scattergl(
                x=x, y=np.zeros(len(values)), mode="markers",
                marker=dict(size=1, color="rgba(0,0,0,0)"),
                # A plain list, not the array: plotly serialises numpy as base64 typed data, which
                # the crosshair below cannot index.
                customdata=[None if v != v else float(v) for v in values],
                name=label, showlegend=False,
                hovertemplate="%{customdata:.4g}<extra>" + label + "</extra>",
            ),
            row=i, col=1,
        )
        figure.update_yaxes(showticklabels=False, ticks="", showgrid=False, zeroline=False,
                            row=i, col=1)
        middle = sum(figure.layout[f"yaxis{i}"].domain) / 2
        rows_at.append({"trace": len(figure.data) - 1, "y": middle, "label": label,
                        "note": len(figure.layout.annotations)})
        figure.data[-2].colorbar.y = middle - 0.012
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
        width=1040, height=600, dragmode="pan", hovermode="x", plot_bgcolor="white",
        paper_bgcolor="white", margin=dict(l=56, r=280, t=6, b=16),
        showlegend=False, hoversubplots="axis", hoverdistance=-1, spikedistance=-1,
        hoverlabel=dict(bgcolor="white", bordercolor="#D7DBE0",
                        font=dict(size=11, color="#3D4653")),
        transition=dict(duration=0),
    )
    return figure, rows_at


CHART_HTML = """
<div id="chart"></div>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<script>
  const figure = __FIGURE__;
  const rows = __ROWS__;
  const gd = document.getElementById("chart");
  const shapes = figure.layout.shapes || [];
  const notes = figure.layout.annotations || [];

  let showing = null;

  Plotly.newPlot(gd, figure.data, figure.layout, __CONFIG__).then(function () {
    gd.on("plotly_hover", function (event) {
      const point = event.points[0];
      const at = point.pointIndex;
      if (at === undefined || at === showing) { return; }
      showing = at;
      const line = {
        type: "line", xref: "x", yref: "paper", x0: point.x, x1: point.x, y0: 0, y1: 1,
        line: { color: "#3D4653", width: 1, dash: "dot" }, layer: "above",
      };
      const labelled = notes.map(function (note) { return Object.assign({}, note); });
      rows.forEach(function (row) {
        const value = gd.data[row.trace].customdata[at];
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


def results_page(job_id: str):
    """Everything one finished job produced, opened from its own address."""
    job = jobs.get(job_id)
    if job is None:
        st.error(f"No job called {job_id}.")
        st.page_link("app.py", label="Design something new")
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
    merged = shortlist.merge(safety, on="aso_sequence", how="left")
    hits = off_targets.groupby("aso_sequence")["distance"].value_counts().unstack(fill_value=0)
    accessibility = merged.get(ACCESSIBILITY_FEATURE)
    binding = merged.get(HYBRIDIZATION_FEATURE)
    # Only worth a column when something in the shortlist actually repeats; a run with no repeated
    # candidate leaves it None and dropna below takes the column out.
    repeat_note = merged["repeat_note"] if "repeat_note" in merged else None
    if repeat_note is not None and not (repeat_note.fillna("").astype(str).str.strip() != "").any():
        repeat_note = None
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
        "**1mm** and **2mm** count genomic hits to a gene other than the target. "
        "**region** is where that site sits in the gene model. **repeat** appears when the same "
        "sequence occurs at more than one site in the target: each copy is drawn at its own "
        "position, but only the first was scored, and the others carry that score."
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
