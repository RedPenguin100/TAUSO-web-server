#!/bin/bash

# 1. Setup the Genome if the DB is missing
if [ ! -f "$TAUSO_DATA_DIR/GRCh38.db" ]; then
    echo "Initial data not found. Running setup-genome..."
    tauso setup-genome
fi

# 2. Setup the Cache if the JSON is missing (even if DB exists)
if [ ! -f "$TAUSO_DATA_DIR/available_genes.json" ]; then
    echo "Gene cache missing. Running cache_genes.py..."
    micromamba run -n base python /app/cache_genes.py
else
    echo "Gene cache found, skipping pre-computation."
fi

PORT="${PORT:-8501}"
echo "Starting Streamlit UI on port $PORT..."
exec streamlit run app.py --server.port=$PORT --server.address=0.0.0.0