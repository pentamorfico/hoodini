# Multi-stage build for Hoodini
FROM condaforge/mambaforge:latest AS builder

# Install build essentials
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    make \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /tmp/build

# Copy environment file
COPY environment.yml .

# Create conda environment
RUN mamba env create -f environment.yml && \
    mamba clean -afy

# Final stage
FROM condaforge/mambaforge:latest

# Install build essentials for packages that need compilation
RUN apt-get update && apt-get install -y \
    gcc \
    g++ \
    make \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy conda environment from builder
COPY --from=builder /opt/conda/envs/hoodini /opt/conda/envs/hoodini

# Copy project files
COPY pyproject.toml README.md ./
COPY src/ ./src/

# Initialize conda for bash
RUN conda init bash

# Activate environment and install the package
RUN /bin/bash -c "source activate hoodini && \
    pip install --no-cache-dir -e ."

# Set environment variables
ENV PATH="/opt/conda/envs/hoodini/bin:${PATH}"
ENV CONDA_DEFAULT_ENV=hoodini

# Activate conda in entrypoint
SHELL ["/bin/bash", "-c"]

# Set working directory
WORKDIR /work

# Create entrypoint script
RUN echo '#!/bin/bash' > /entrypoint.sh && \
    echo 'source /opt/conda/etc/profile.d/conda.sh' >> /entrypoint.sh && \
    echo 'conda activate hoodini' >> /entrypoint.sh && \
    echo 'exec "$@"' >> /entrypoint.sh && \
    chmod +x /entrypoint.sh

ENTRYPOINT ["/entrypoint.sh"]
CMD ["hoodini", "--help"]
