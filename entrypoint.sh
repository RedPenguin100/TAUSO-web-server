#!/bin/bash

set -e

mkdir -p "$TAUSO_DATA_DIR"

# 1. RACCESS INSTALLATION
# Because raccess cannot be legally redistributed, we compile it at runtime.
# We check for the executable explicitly because it lives outside the persistent volume.
if [ ! -f "$RACCESS_EXE" ]; then
    echo "raccess binary not found. Compiling from source..."
    tauso setup-raccess
else
    echo "raccess is already compiled."
fi

# 2. Full TAUSO Database & Weights Initialization (Persistent Volume Check)
if [ ! -f "$TAUSO_DATA_DIR/.tauso_initialized_v2" ]; then
    echo "Initial data or weights not found. Running full TAUSO setup pipeline..."

    tauso setup-genome
    tauso setup-bowtie   # Very slow, can take 1~2 hours on slow single-threaded CPUs
    tauso setup-mrna-halflife
    tauso setup-attract
    tauso setup-riboseq
    tauso setup-depmap
    tauso build-cell-context   # cohort expression + CAI weights + tGCN (default cohort)
    tauso setup-rrna

    # Create the sentinel file so this block is skipped on future reboots
    touch "$TAUSO_DATA_DIR/.tauso_initialized_v2"
    echo "Full TAUSO initialization complete!"
else
    echo "TAUSO databases and weights found. Skipping initialization."
fi

# 3. Setup the Streamlit App Cache
if [ ! -f "$TAUSO_DATA_DIR/available_genes.json" ]; then
    echo "Gene cache missing. Running cache_genes.py..."
    micromamba run -n base python /app/cache_genes.py
else
    echo "Gene cache found, skipping pre-computation."
fi

# 4. Start the Webserver
PORT="${PORT:-8501}"
echo "Starting Streamlit UI on port $PORT..."
exec streamlit run app.py --server.port=$PORT --server.address=0.0.0.0