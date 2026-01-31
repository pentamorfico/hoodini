---
title: Installation
description: How to install Hoodini on your system
---

import { Tabs, Steps, Callout } from 'nextra/components'

# Installation

Hoodini can be installed using several methods. Choose the one that best fits your workflow.

<Callout type="info" emoji="⏱️">
  **Note:** Hoodini requires downloading several databases on first run (~35GB uncompressed). 
  Initial setup takes 5+ minutes depending on your internet connection and machine specs.
</Callout>

## Installation Methods

<Tabs items={['Mamba (Recommended)', 'Pixi', 'pip', 'Docker']}>
  <Tabs.Tab>
    **Mamba** is the recommended method for most users. It handles all dependencies automatically.

    <Callout type="info">
      Hoodini is not yet on Bioconda. For now, clone the repository and install from source.
    </Callout>

    <Steps>
      ### Install Mamba
      
      If you don't have Mamba installed, get it via [Miniforge](https://github.com/conda-forge/miniforge):
      
      ```bash
      curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
      bash Miniforge3-$(uname)-$(uname -m).sh
      ```

      ### Clone the repository
      
      ```bash
      git clone https://github.com/pentamorfico/hoodini.git
      cd hoodini
      ```

      ### Create environment and install
      
      ```bash
      mamba env create -f environment.yml
      mamba activate hoodini
      pip install -e .
      ```

      ### Verify installation
      
      ```bash
      hoodini --help
      ```
    </Steps>
  </Tabs.Tab>

  <Tabs.Tab>
    **Pixi** is a fast, modern package manager that's great for reproducible environments.

    <Callout type="info">
      Hoodini is not yet on conda-forge. For now, clone the repository and use the environment file.
    </Callout>

    <Steps>
      ### Install Pixi
      
      ```bash
      curl -fsSL https://pixi.sh/install.sh | bash
      ```

      ### Clone and install
      
      ```bash
      git clone https://github.com/pentamorfico/hoodini.git
      cd hoodini
      pixi install
      ```

      ### Run Hoodini
      
      ```bash
      pixi run hoodini --help
      ```
    </Steps>
  </Tabs.Tab>

  <Tabs.Tab>
    **pip** installation works but requires external tools to be installed separately.

    <Callout type="warning">
      The pip installation requires you to manually install external dependencies like
      `mmseqs2`, `diamond`, `mafft`, `veryfasttree`, `padloc`, `defense-finder`, etc.
      **We strongly recommend using Mamba instead.**
    </Callout>

    <Steps>
      ### Clone the repository
      
      ```bash
      git clone https://github.com/pentamorfico/hoodini.git
      cd hoodini
      ```

      ### Install external tools first
      
      Make sure all bioinformatics tools from `environment.yml` are available in your PATH.

      ### Install Hoodini
      
      ```bash
      pip install -e .
      ```

      ### Verify installation
      
      ```bash
      hoodini --help
      ```
    </Steps>
  </Tabs.Tab>

  <Tabs.Tab>
    **Docker** provides an isolated environment with all dependencies included.

    <Steps>
      ### Pull the image
      
      ```bash
      docker pull pentamorfico/hoodini:latest
      ```

      ### Run Hoodini
      
      ```bash
      docker run -v $(pwd):/data pentamorfico/hoodini:latest hoodini --help
      ```

      ### Run with your data
      
      Mount your data directory to `/data` inside the container:
      
      ```bash
      docker run -v /path/to/your/data:/data pentamorfico/hoodini:latest \
        hoodini run -i /data/input.fasta -o /data/output
      ```
    </Steps>
  </Tabs.Tab>
</Tabs>

## Database Download

On first run, Hoodini will download the required databases. Here are the approximate sizes:

<Callout type="warning" emoji="💾">
  **Storage requirements:** Make sure you have at least **40GB** of free disk space for all databases.
</Callout>

| Database | Size | Description |
|----------|------|-------------|
| eggNOG-mapper | ~24 GB | Functional annotation (largest) |
| MetaCerberus HMMs | ~4.5 GB | Domain annotation (Pfam, COG, KEGG, etc.) |
| Contig lengths | ~3.8 GB | NCBI assembly metadata |
| geNomad | ~1.4 GB | Virus/plasmid detection |
| PADLOC | ~1 GB | Defense system detection |
| DefenseFinder | ~350 MB | Defense system detection |

You can skip specific databases during download:

```bash
hoodini download databases --skip-emapper --skip-genomad
```

## Verify Installation

After installation, verify everything works:

```bash
# Check version
hoodini --version

# Run with test data (will download databases on first run)
hoodini run -i test.fasta -o test_output --threads 4
```

<Callout type="info">
  **Tip:** Use the `--threads` flag to speed up processing. Hoodini scales well with multiple cores.
</Callout>
