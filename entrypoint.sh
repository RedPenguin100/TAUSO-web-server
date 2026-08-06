#!/bin/bash

set -e

mkdir -p "$TAUSO_DATA_DIR"

# 1. raccess cannot be legally redistributed, so it is built from source into the persistent
# volume. The command compares the pinned commit and the installed binary and no-ops when both
# already match, so it is cheap to run on every boot.
tauso setup-raccess

# 2. The trained booster is not shipped in the package; it is fetched from Zenodo into the
# persistent volume. The command verifies the md5 and skips the download when it is already there.
tauso setup-model

# 3. Full TAUSO Database & Weights Initialization (Persistent Volume Check)
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

# 4. Setup the Streamlit App Cache
if [ ! -f "$TAUSO_DATA_DIR/available_genes.json" ]; then
    echo "Gene cache missing. Running cache_genes.py..."
    micromamba run -n base python /app/cache_genes.py
else
    echo "Gene cache found, skipping pre-computation."
fi

# 5. Start the Webserver
PORT="${PORT:-8501}"
echo "Starting Streamlit UI on port $PORT..."
exec streamlit run app.py --server.port=$PORT --server.address=0.0.0.0