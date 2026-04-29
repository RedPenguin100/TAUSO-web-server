#!/bin/bash

set -e

mkdir -p "$TAUSO_DATA_DIR"

# 1. RACCESS INSTALLATION
# Because raccess cannot be legally redistributed, we compile it at runtime.
# We check for the executable explicitly because it lives outside the persistent volume.
if [ ! -f "$RACCESS_EXE" ]; then
    echo "raccess binary not found. Compiling from source..."
    tauso install-raccess
else
    echo "raccess is already compiled."
fi

# 2. Full TAUSO Database & Weights Initialization (Persistent Volume Check)
if [ ! -f "$TAUSO_DATA_DIR/.tauso_fully_initialized" ]; then
    echo "Initial data or weights not found. Running full TAUSO setup pipeline..."

    tauso setup-genome
    tauso setup-mrna-halflife
    tauso setup-attract
    tauso setup-depmap
    tauso setup-bowtie
    tauso add-cell \
        HEPG2 SNU449 HELA A431 SKMEL28 SHSY5Y U251MG NCIH929 \
        KMS11 NCIH460 SKNAS SKNSH KARPAS299 HEP3B217 THP1 \
        LNCAPCLONEFGC T24 A549 VCAP HUH7 JURKAT SKOV3 K562 \
        A172 PC3 MCF7 SW872 G361 HEK293 MM1S

    tauso build-omics
    tauso build-cai-weights

    tauso setup-bowtie # Very slow, can take 1~2 hours on slow single threaded CPUs

    # Create the sentinel file so this block is skipped on future reboots
    touch "$TAUSO_DATA_DIR/.tauso_fully_initialized"
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