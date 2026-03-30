# Use the official micromamba image [cite: 2]
FROM mambaorg/micromamba:1.5-jammy

# Switch to root to install git
USER root
RUN apt-get update && \
    apt-get install -y git && \
    rm -rf /var/lib/apt/lists/*

# Switch back to the default mamba user
USER $MAMBA_USER

# Set up the data directory environment variable
ENV TAUSO_DATA_DIR=/home/mambauser/.tauso_data

WORKDIR /app

# Clone the repository directly from GitHub
ARG CACHEBUST=1
RUN git clone https://github.com/RedPenguin100/TAUSO.git .

# ==========================================
# DEPENDENCIES RE-ENABLED FOR TAUSO CLI
# ==========================================
# Restore the actual environment so `tauso setup-genome` works [cite: 3]
RUN micromamba install -y -n base -f environment.yml && \
    micromamba clean --all --yes

RUN micromamba run -n base pip install -e .

# Install the UI dependencies on top of the base environment [cite: 4]
# We use pip here to avoid conda-forge C++ library conflicts with pyarrow
RUN micromamba run -n base pip install streamlit watchdog

# Add the protobuf fallback environment variable just to be absolutely safe
ENV PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python

# Expose the default Streamlit port
EXPOSE 8501

# ==========================================
# ENTRYPOINT & APP SETUP
# ==========================================
USER root

# Copy your local UI scripts into the container
COPY app.py /app/app.py
COPY cache_genes.py /app/cache_genes.py

# Copy and set permissions for the entrypoint
COPY entrypoint.sh /app/entrypoint.sh
RUN chmod +x /app/entrypoint.sh

# Switch back to the mamba user
USER $MAMBA_USER
# ==========================================

# Set the entrypoint to run the script inside the conda environment
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "/app/entrypoint.sh"]