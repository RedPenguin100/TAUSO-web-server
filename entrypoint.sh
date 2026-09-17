#!/bin/bash

set -e

mkdir -p "$TAUSO_DATA_DIR"

# 1. The trained booster is not shipped in the package; it is fetched into the persistent volume.
# A failure here is reported and survived rather than being fatal: without a model a design job
# fails and mails the submitter, which is a great deal better than the site never starting.
if ! tauso setup-model; then
    echo "WARNING: no scoring model. The site will run and designs will fail until one is in place."
fi

# The tRNA gene counts, seeded from the image rather than scraped. `tauso setup-tgcn` reads them
# off GtRNAdb, and when that scrape fails `set -e` takes the container down -- which Docker then
# restarts into the same failure, forever, with nothing in the log that says why. The file is 310
# bytes, its sha256 is pinned inside tauso, and it is the same for every install.
if [ ! -f "$TAUSO_DATA_DIR/human_tgcn_hsapi38.csv" ] && [ -f /app/assets/human_tgcn_hsapi38.csv ]; then
    echo "Seeding tRNA gene counts from the image."
    cp /app/assets/human_tgcn_hsapi38.csv "$TAUSO_DATA_DIR/human_tgcn_hsapi38.csv"
fi

# 2. Full TAUSO Database & Weights Initialization (Persistent Volume Check)
if [ ! -f "$TAUSO_DATA_DIR/.tauso_initialized_v2" ]; then
    echo "Initial data or weights not found. Running full TAUSO setup pipeline..."

    tauso setup-genome
    tauso setup-bowtie   # Very slow, can take 1~2 hours on slow single-threaded CPUs
    tauso setup-mrna-halflife
    tauso setup-attract
    tauso setup-depmap
    tauso build-cell-context   # cohort expression + CAI weights + tGCN (default cohort)
    # Transcript-level expression is a hard requirement now, not a warning: without it the
    # first design fails outright. It pulls a ~3 GB DepMap table the first time.
    tauso build-cohort-transcript-expression
    # The cohort-wide mean the general off-target step ranks genes by. Built once: it depends on
    # the DepMap matrix and on which genes the annotation calls valid, neither of which changes
    # from run to run, and the step raises rather than computing it if the table is absent.
    tauso build-general-expression
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