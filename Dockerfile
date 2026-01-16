# Multi-stage build for Hoodini
FROM condaforge/mambaforge:latest AS builder

# Avoid interactive tzdata prompts
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

# Install build essentials
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    make \
    tzdata \
    && ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
    && echo $TZ > /etc/timezone \
    && dpkg-reconfigure -f noninteractive tzdata \
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

# Avoid interactive tzdata prompts
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

# Install build essentials (Firefox from mamba brings GUI/X11 deps)
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    g++ \
    make \
    tzdata \
    && ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
    && echo $TZ > /etc/timezone \
    && dpkg-reconfigure -f noninteractive tzdata \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# Copy conda environment from builder
COPY --from=builder /opt/conda/envs/hoodini /opt/conda/envs/hoodini

# Copy project files
COPY pyproject.toml README.md ./
COPY src/ ./src/

# Install the package in the env (no need for conda init)
RUN /bin/bash -lc "source /opt/conda/etc/profile.d/conda.sh && \
    conda activate hoodini && \
    pip install --no-cache-dir -e ."

# Install Playwright Firefox with system dependencies
RUN /bin/bash -lc "source /opt/conda/etc/profile.d/conda.sh && \
    conda activate hoodini && \
    playwright install --with-deps firefox"

# Set environment variables
ENV PATH="/opt/conda/envs/hoodini/bin:${PATH}"
ENV CONDA_DEFAULT_ENV=hoodini

# Create entrypoint script
RUN echo '#!/bin/bash' > /entrypoint.sh && \
    echo 'source /opt/conda/etc/profile.d/conda.sh' >> /entrypoint.sh && \
    echo 'conda activate hoodini' >> /entrypoint.sh && \
    echo 'exec "$@"' >> /entrypoint.sh && \
    chmod +x /entrypoint.sh

WORKDIR /work
ENTRYPOINT ["/entrypoint.sh"]
CMD ["hoodini", "--help"]
