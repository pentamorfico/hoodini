---
title: Installation
description: How to install Hoodini on your system
---

import { Tabs, Steps, Callout } from 'nextra/components'

# Installation

Hoodini can be installed using several methods. Choose the one that best fits your workflow.

<Callout type="info" emoji="⏱️">
  **Note:** Hoodini requires downloading several databases on first run (~40GB total). 
  Initial setup may take 30-60 minutes depending on your internet connection.
</Callout>

## Installation Methods

<Tabs items={['Mamba (Recommended)', 'Pixi', 'pip', 'Docker']}>
  <Tabs.Tab>
    **Mamba** is the recommended method for most users. It handles all dependencies automatically.

    <Steps>
      ### Install Mamba
      
      If you don't have Mamba installed, get it via [Miniforge](https://github.com/conda-forge/miniforge):
      
      ```bash
      curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
      bash Miniforge3-$(uname)-$(uname -m).sh
      ```

      ### Create Hoodini environment
      
      ```bash
      mamba create -n hoodini -c conda-forge -c bioconda hoodini
      ```

      ### Activate and run
      
      ```bash
      mamba activate hoodini
      hoodini --help
      ```
    </Steps>
  </Tabs.Tab>

  <Tabs.Tab>
    **Pixi** is a fast, modern package manager that's great for reproducible environments.

    <Steps>
      ### Install Pixi
      
      ```bash
      curl -fsSL https://pixi.sh/install.sh | bash
      ```

      ### Add Hoodini to your project
      
      ```bash
      pixi add hoodini
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
      The pip installation requires you to manually install external dependencies: 
      `mmseqs2`, `diamond`, `mafft`, `fasttree`, `iqtree2`, and `muscle`.
    </Callout>

    <Steps>
      ### Install external tools first
      
      Make sure `mmseqs2`, `diamond`, `mafft`, `fasttree`, `iqtree2`, and `muscle` are available in your PATH.

      ### Install Hoodini
      
      ```bash
      pip install hoodini
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
  **Storage requirements:** Make sure you have at least **50GB** of free disk space.
</Callout>

| Database | Size | Description |
|----------|------|-------------|
| MMseqs2 taxonomy | ~15GB | GTDB taxonomic classification |
| Diamond (Pfam) | ~8GB | Protein domain annotation |
| Defense Finder models | ~2GB | Defense system detection |
| CRISPRCasTyper models | ~1GB | CRISPR-Cas typing |
| Other models | ~1GB | Additional annotation databases |

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
