import streamlit as st
import io
import sqlite3
import os
import json
import tauso
from pipeline_runner import trigger_background_job, JobConfig
from pyfaidx import Fasta

import hashlib
from Bio import SeqIO


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


def parse_fasta_input(raw_text: str):
    """Robust FASTA parsing using Biopython's SeqIO with sequence hashing."""
    raw_text = raw_text.strip()

    # Extract the name and sequence safely
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

    # ---> ADD THIS EXACT LINE HERE <---
    cleaned_seq = cleaned_seq.upper().replace('T', 'U')

    # --- NEW: Hashing Logic ---
    # Create a deterministic MD5 hash of the actual sequence text
    seq_hash = hashlib.md5(cleaned_seq.encode('utf-8')).hexdigest()[:8]

    # Combine the user's name with the unique hash
    unique_name = f"{base_name}_{seq_hash}"

    return unique_name, cleaned_seq

def main():
    # Header Section
    st.title("🧬 TAUSO")
    st.markdown("**Predictive modeling for Antisense Oligonucleotide (ASO) efficacy.**")
    st.divider()

    st.write("Please choose your input method below to begin the analysis pipeline.")

    # Fetch the genes on load
    gene_list = fetch_genes()

    # Using tabs to cleanly separate the three input methods
    tab1, tab2, tab3 = st.tabs(["🗄️ Database Selection", "📝 Paste Sequence", "📁 Upload FASTA"])

    with tab1:
        selected_gene = st.selectbox(
            "Select a Target Gene:",
            options=["-- Choose a gene --"] + gene_list
        )

    with tab2:
        pasted_sequence = st.text_area(
            "Paste sequence or FASTA text:",
            height=150,
            placeholder=">Sequence_Name\nATGCGTACGTTAG..."
        )

    with tab3:
        uploaded_file = st.file_uploader("Upload FASTA File:", type=['fasta', 'fa', 'txt'])

    st.divider()

    selected_cell_line = st.selectbox("Cell Line (Optional)", ["None", "HEK293", "HeLa"])
    user_cell_line = None if selected_cell_line == "None" else selected_cell_line


    # <-- 2. NEW EMAIL INPUT
    user_email = st.text_input("📧 Notification Email:", placeholder="Enter your email to receive results...")

    # Form submission logic
    if st.button("Initialize Pipeline", type="primary", use_container_width=True):

        # <-- 3. NEW EMAIL VALIDATION
        if not user_email or "@" not in user_email:
            st.error("Please provide a valid email address before proceeding.")
            return

        target_data = None
        target_mrna_name = None
        source_info = None

        try:
            # Logic hierarchy to determine which input to use if they clicked around
            if uploaded_file is not None:
                raw_text = uploaded_file.getvalue().decode("utf-8")
                target_mrna_name, target_data = parse_fasta_input(raw_text)
                source_info = f"Uploaded File: {uploaded_file.name}"

            elif pasted_sequence.strip():
                target_mrna_name, target_data = parse_fasta_input(pasted_sequence.strip())
                source_info = "Pasted Sequence"

            elif selected_gene != "-- Choose a gene --" and not selected_gene.startswith("--"):
                # Defer the fetch to the worker: set the name, but leave the data empty
                target_mrna_name = selected_gene
                target_data = ""
                source_info = f"Selected Gene: {selected_gene}"

            else:
                st.error("Please select a valid gene, paste a sequence, or upload a FASTA file before proceeding.")
                return

        except ValueError as e:
            # Catches the formatting errors from parse_fasta_input and displays them safely
            st.error(f"Input Error: {str(e)}")
            return

        # Initialize the config with the parsed info
        config = JobConfig(
            target_mrna_name=target_mrna_name,
            target_data=target_data,
            source_info=source_info,
            user_email=user_email,
            cell_line=user_cell_line,
        )

        # --- Processing UI ---
        st.success("Input accepted! Starting analysis...")

        with st.status("Submitting job to TAUSO Pipeline...", expanded=True) as status:
            st.write(f"**Source:** {source_info}")
            st.write(f"**Email:** {user_email}")
            st.write("Initializing target sequence parameters...")

            # <-- 4. NEW BACKGROUND TRIGGER
            # This replaces the blocking tauso logic that was here before
            trigger_background_job(config)


            st.write("Extracting context and queuing job...")
            st.write("Job handed off to background processor.")

            status.update(label="Job successfully queued! Check your email.", state="complete", expanded=False)

        # Show a preview of the data/results
        st.subheader("Sequence Preview")
        preview_text = target_data[:200] + "..." if len(target_data) > 200 else target_data
        st.code(preview_text, language="text")

if __name__ == "__main__":
    main()