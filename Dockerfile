# Use the official micromamba image
FROM mambaorg/micromamba:1.5-jammy

# Switch to root to install git and uv
USER root
RUN apt-get update && \
    apt-get install -y git build-essential zlib1g-dev wget && \
    rm -rf /var/lib/apt/lists/*

# Grab the ultra-fast uv installer directly from Astral's official image
COPY --from=ghcr.io/astral-sh/uv:latest /uv /uvx /bin/

# Switch back to the default mamba user
USER $MAMBA_USER

# ==========================================
# 1. WORKSPACE SETUP (Code & Data)
# ==========================================
# Define the overarching workspace and the specific data subfolder
ENV TAUSO_WORKSPACE=/home/mambauser/tauso_workspace
ENV TAUSO_DATA_DIR=$TAUSO_WORKSPACE/data

# Set working directory strictly for the TAUSO source code
WORKDIR $TAUSO_WORKSPACE/code

# Pin the TAUSO source to a specific main commit for reproducible builds.
ARG TAUSO_COMMIT=9bb429edff608d17abb4543cc5fa37112e732604
RUN git init -q . && \
    git remote add origin https://github.com/RedPenguin100/TAUSO.git && \
    git config core.sparseCheckout true && \
    printf '/*\n!/notebooks/\n!/tests/\n' > .git/info/sparse-checkout && \
    git fetch --depth 1 origin ${TAUSO_COMMIT} && \
    git checkout -q FETCH_HEAD && \
    git submodule update --init --recursive

# Install dependencies and the TAUSO package natively
RUN micromamba install -y -n base -f environment.yml && \
    micromamba clean --all --yes

COPY requirements.txt ./
RUN micromamba run -n base uv pip install --system . -r requirements.txt

# Add the protobuf fallback environment variable
ENV PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python

# Expose the default Streamlit port
EXPOSE 8501

# ==========================================
# 2. ENTRYPOINT & UI APP SETUP
# ==========================================
USER root

# Move to a pristine /app directory for the Streamlit UI
WORKDIR /app

# Copy your local UI scripts into the container
COPY *.py entrypoint.sh ./

# Set permissions for the entrypoint
RUN chmod +x /app/entrypoint.sh

# Switch back to the mamba user
USER $MAMBA_USER
# ==========================================

# Set the entrypoint to run the script inside the conda environment
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "/app/entrypoint.sh"]