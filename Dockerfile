# Use the official micromamba image
FROM mambaorg/micromamba:1.5-jammy

# Switch to root to install git and uv
USER root
RUN apt-get update && \
    apt-get install -y git build-essential zlib1g-dev wget && \
    rm -rf /var/lib/apt/lists/*

# Grab the ultra-fast uv installer directly from Astral's official image
# Pinned like everything else: uv resolves the pip layer, so a floating version can change
# what gets installed even when requirements.txt has not moved.
COPY --from=ghcr.io/astral-sh/uv:0.12.10 /uv /uvx /bin/

# Switch back to the default mamba user
USER $MAMBA_USER

# ==========================================
# 1. WORKSPACE SETUP (Code & Data)
# ==========================================
# The workspace holds the TAUSO source. The data directory is the mount point the compose
# file binds the persistent volume to, and the default app.py and cache_genes.py fall back
# to, so all three name the same path.
ENV TAUSO_WORKSPACE=/home/mambauser/tauso_workspace
ENV TAUSO_DATA_DIR=/home/mambauser/.tauso_data

# Dependencies are installed BEFORE the source, so bumping TAUSO_COMMIT below does not drag the
# whole conda stack and the 1.3 GB of Streamlit/plotly/biopython through a reinstall with it.
# environment.yml is vendored here for that reason: read from the clone it would be a child of the
# commit, and every bump would re-solve the environment. The build fails loudly further down if
# the vendored copy has drifted from the pinned revision's own.
WORKDIR $TAUSO_WORKSPACE/build
COPY environment.yml requirements.txt ./
RUN micromamba install -y -n base -f environment.yml && \
    micromamba clean --all --yes
RUN micromamba run -n base uv pip install --system -r requirements.txt

# Set working directory strictly for the TAUSO source code
WORKDIR $TAUSO_WORKSPACE/code

# Pin the TAUSO source to a specific main commit for reproducible builds.
ARG TAUSO_COMMIT=79df9520f20311659415e61bd1f36794a9cd5e26
RUN git init -q . && \
    git remote add origin https://github.com/RedPenguin100/TAUSO.git && \
    git config core.sparseCheckout true && \
    printf '/*\n!/notebooks/\n!/tests/\n' > .git/info/sparse-checkout && \
    git fetch --depth 1 origin ${TAUSO_COMMIT} && \
    git checkout -q FETCH_HEAD && \
    git submodule update --init --recursive

# The vendored spec is what the environment above was built from. If this commit asks for anything
# different, stop: the alternative is an image whose conda stack silently does not match its source.
RUN cmp -s $TAUSO_WORKSPACE/build/environment.yml environment.yml || { \
        echo "ERROR: vendored environment.yml differs from TAUSO ${TAUSO_COMMIT}."; \
        echo "Refresh it:  git -C <tauso> show ${TAUSO_COMMIT}:environment.yml > environment.yml"; \
        diff $TAUSO_WORKSPACE/build/environment.yml environment.yml || true; \
        exit 1; \
    }

# Only the TAUSO package itself; its dependencies came from the layers above.
RUN micromamba run -n base uv pip install --system .

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
COPY .streamlit ./.streamlit
COPY components ./components
# The gene list is written here at boot, so the runtime user owns the directory.
RUN chown -R $MAMBA_USER /app/components

# Set permissions for the entrypoint
RUN chmod +x /app/entrypoint.sh

# Switch back to the mamba user
USER $MAMBA_USER
# ==========================================

# Set the entrypoint to run the script inside the conda environment
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "/app/entrypoint.sh"]