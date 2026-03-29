#!/bin/bash

# ==========================================
# HEAVY DATA SETUP COMMENTED OUT FOR UI DEV
# ==========================================
# if [ ! -f "$TAUSO_DATA_DIR/GRCh38.db" ]; then
#     echo "Initial data not found in volume. Running setup scripts..."
#     tauso setup-genome
#     tauso setup-bowtie
#     tauso setup-mrna-halflife
#     tauso setup-attract
#     tauso setup-depmap
#     tauso add-cell HEPG2 SNU449 HELA A431 SKMEL28 SHSY5Y U251MG NCIH929 KMS11 NCIH460 SKNAS SKNSH KARPAS299 HEP3B217 THP1 LNCAPCLONEFGC T24 A549 VCAP HUH7 JURKAT SKOV3 K562 A172 PC3 MCF7 SW872 G361 HEK293 MM1S
#     tauso build-omics
#     tauso build-cai-weights
# else
#     echo "TAUSO data found in mounted volume, skipping setup."
# fi
# ==========================================

# Render injects a $PORT environment variable. We use 8501 as a fallback for local dev.
PORT="${PORT:-8501}"

echo "Starting Streamlit UI on port $PORT..."
exec streamlit run app.py --server.port=$PORT --server.address=0.0.0.0