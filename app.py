import streamlit as st
import io
import sqlite3
import os
from tauso.genome.read_human_genome import get_locus_to_data_dict

# Configure the page layout and title
st.set_page_config(
    page_title="TAUSO | ASO Efficacy Predictor",
    page_icon="🧬",
    layout="centered"
)

@st.cache_data(ttl=3600)
def fetch_genes():
    db_path = os.path.join(os.environ.get("TAUSO_DATA_DIR", "/home/mambauser/.tauso_data"), "available_genes.json")
    if not os.path.exists(db_path):
        return ["-- Database not initialized --"]

    with open(db_path, "r") as f:
        return json.load(f)


def main():
    # Header Section
    st.title("🧬 TAUSO")
    st.markdown("**Predictive modeling for Antisense Oligonucleotide (ASO) efficacy.**")
    st.divider()

    st.write("Please select a target gene from the database OR upload a custom FASTA file to begin the analysis pipeline.")

    # Fetch the genes on load
    gene_list = fetch_genes()

    # Input Section using Streamlit columns for layout
    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Database Selection")
        selected_gene = st.selectbox(
            "Select a Target Gene:",
            options=["-- Choose a gene --"] + gene_list
        )

    with col2:
        st.subheader("Custom Upload")
        uploaded_file = st.file_uploader("Upload FASTA File:", type=['fasta', 'fa', 'txt'])

    st.divider()

    # Form submission logic
    if st.button("Initialize Pipeline", type="primary", use_container_width=True):
        target_data = None
        source_info = None

        # Determine which input method the user utilized
        if uploaded_file is not None:
            stringio = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
            target_data = stringio.read()
            source_info = f"Uploaded File: {uploaded_file.name}"

        elif selected_gene != "-- Choose a gene --" and not selected_gene.startswith("--"):
            target_data = f"Simulating database fetch for sequence: {selected_gene}..."
            source_info = f"Selected Gene: {selected_gene}"

        else:
            st.error("Please select a valid gene or upload a FASTA file before proceeding.")
            return

        # --- Processing UI ---
        st.success("Input accepted! Starting analysis...")

        with st.status("Running TAUSO Pipeline...", expanded=True) as status:
            st.write(f"**Source:** {source_info}")
            st.write("Initializing target sequence parameters...")

            # TODO: Call your actual TAUSO Python functions here.
            # Example:
            # results = tauso.predict_efficacy(target_data)

            st.write("Extracting context and finding genes...")
            st.write("Calculating mRNA half-life and cAI weights...")

            status.update(label="Analysis Complete!", state="complete", expanded=False)

        # Show a preview of the data/results
        st.subheader("Sequence Preview")
        preview_text = target_data[:200] + "..." if len(target_data) > 200 else target_data
        st.code(preview_text, language="text")

if __name__ == "__main__":
    main()