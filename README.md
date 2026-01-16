<h3 align="center">hoodini 🦉🎩: large-scale gene neighborhood analyses that feel like magic</h3>

---

<img src="docs/hoodini_logo_github.svg" alt="hoodini logo" width="220" align="left" style="border:0; margin:0 16px 8px 0; display:block;" />

<p align="justify">
   <b>hoodini</b> is a large-scale gene-centric comparative genomics toolkit that fetches public assemblies, extracts gene neighborhoods, runs pairwise protein and nucleotide comparisons, annotates neighborhoods with extra tools, and builds trees to help interpret genomic context.
</p>

---

<p align="center">
    <img src="docs/hoodini-viz-export - 2026-01-14T051133.869.svg" alt="Hoodini Visualization Example" width="800"/>
</p>


<div align = center>

[<kbd> <br> Click here for an interactive demo <br> </kbd>][KBD]

</div>


[KBD]: https://storage.hoodini.bio/hoodini-demo.html


## Introduction

Hoodini’s primary goal is to make large-scale genomic context analysis fast, practical, and accessible. It’s built to fetch and process thousands of gene neighborhoods in minutes, turning what used to be hours or days of manual work into an interactive workflow.

With GPU-accelerated, real-time interactive  visualization, Hoodini lets biologists across disciplines explore bacterial operons, defense systems, and mobile genetic elements at scale, revealing patterns, co-localization signals, and evolutionary signatures that drive discovery and deepen our understanding of prokaryotic genomes.

Beyond exploration, Hoodini offers extensive customization options (including flexible styling, multiple palettes, and fine-grained layout controls) to generate publication-ready, high-quality figures exported for presentations, manuscripts, and supplementary materials.

This README documents installation, the command-line interface and all available options, produced output files, and a few examples to get you started.

### Key Features

- **Automated data retrieval**: Fetches assemblies and annotations directly from NCBI using protein or nucleotide accessions
- **Neighborhood extraction**: Extracts configurable genomic windows around target genes
- **Protein clustering**: Groups homologous proteins across neighborhoods for synteny comparison  
- **Pairwise comparisons**: Computes protein (AAI) and nucleotide (ANI) similarities
- **Tree construction**: Builds phylogenetic trees from AAI or ANI distances
- **Defense system annotation**: Integrates PADLOC, DefenseFinder, CCTyper for antiphage systems
- **Mobile element detection**: Identifies prophages and plasmids via geNomad
- **Interactive visualization**: Generates self-contained HTML viewer with aligned neighborhoods and trees

### Use Cases

- Comparative analysis of gene neighborhoods across thousands of genomes
- Evolutionary analysis of genomic islands, defense systems, and mobile elements
- Identification of conserved gene clusters and syntenic regions
- Phylogenetic contextualization of protein families

## Quick Start

```bash
# Single protein query
hoodini run --input WP_012345678.1 --output results

# With protein comparisons and phylogenetic tree
hoodini run --input proteins.txt --output results --prot-links --tree-mode aai_tree

# Full analysis with annotations
hoodini run --input proteins.txt --output results \
  --prot-links --tree-mode aai_tree \
  --padloc --deffinder --cctyper --genomad \
  --num-threads 16
```

## Installation

### Mamba

```bash
mamba env create -f environment.yml
mamba activate hoodini
pip install -e .
hoodini download databases
```

### Pixi

```bash
pixi install
pixi shell
hoodini download databases
```

### Docker

```bash
docker volume create hoodini-data

# Download databases (first time only)
docker run --rm -v hoodini-data:/app/src/hoodini/data \
  pentamorfico/hoodini:latest hoodini download databases

# Run analysis
docker run --rm -v hoodini-data:/app/src/hoodini/data -v $(pwd):/work \
  pentamorfico/hoodini:latest hoodini run --input /work/proteins.txt --output /work/results
```

## Usage

```
hoodini run      Run the main pipeline
hoodini download Download required databases
```

### Input Options

| Option | Description |
|--------|-------------|
| `--input ID\|FILE` | Single accession (e.g., `WP_012345678.1`) or file with one accession per line |
| `--inputsheet FILE` | TSV with accessions and custom metadata columns |

### Output Options

| Option | Description |
|--------|-------------|
| `--output DIR` | Output directory |
| `--force` | Overwrite existing output |
| `--keep` | Retain intermediate files |

### Neighborhood Extraction

| Option | Description |
|--------|-------------|
| `--win-mode` | `win_genes` (gene count) or `win_nts` (nucleotide distance) |
| `--win INT` | Window size (default: 10 genes or 10000 nt) |
| `--min-win INT` | Minimum genes required per side |
| `--sorfs` | Re-annotate small ORFs in extracted regions |

### Pairwise Comparisons

| Option | Description |
|--------|-------------|
| `--prot-links` | Compute all-vs-all protein similarities |
| `--nt-links` | Compute pairwise nucleotide alignments |
| `--nt-aln-mode` | Alignment method: `blastn`, `fastani`, `minimap2` |
| `--clust-method` | Protein clustering algorithm |

### Tree Construction

| Option | Description |
|--------|-------------|
| `--tree-mode` | `aai_tree` (amino acid identity) or `ani_tree` (nucleotide identity) |
| `--tree-file FILE` | Use precomputed Newick tree |
| `--aai-mode` | Tree algorithm: `nj` (neighbor-joining) or `hyper` |
| `--aai-subset-mode` | Proteins for tree: `target_region`, `target_prot`, `window` |

### Functional Annotations

| Option | Description |
|--------|-------------|
| `--padloc` | Detect defense systems with PADLOC |
| `--deffinder` | Detect defense systems with DefenseFinder |
| `--cctyper` | Type CRISPR-Cas systems |
| `--genomad` | Identify mobile genetic elements |
| `--ncrna` | Predict non-coding RNAs with Infernal |
| `--domains LIST` | Search domain databases (comma-separated) |

### Performance

| Option | Description |
|--------|-------------|
| `--num-threads INT` | Parallel threads (default: 4) |
| `--max-concurrent-downloads INT` | Concurrent NCBI downloads |
| `--api-key KEY` | NCBI API key (increases rate limits) |

## Output Structure

```
results/
├── assembly_list.txt           # Downloaded assembly accessions
├── assembly_folder/            # Raw assemblies (*.gbff / *.fna, *.gff)
├── all_neigh.tsv               # All neighborhood coordinates
├── neighborhood/
│   └── neighborhoods.fasta     # Extracted neighborhood sequences
├── target_prots.fasta          # Target proteins for clustering
├── target_prots.aln            # Protein alignment (if clustering)
├── pairwise_aa.tsv             # Protein similarity hits (if --prot-links)
├── aai_matrix.tsv              # AAI distance matrix (if --tree-mode aai_tree)
├── nt_links.tsv                # Nucleotide alignments (if --nt-links)
├── ani_matrix.tsv              # ANI distance matrix (if --tree-mode ani_tree)
├── tree.nwk                    # Phylogenetic tree
├── records.csv                 # Input records with metadata
├── domains.tsv                 # Domain annotations (if --domains)
├── cctyper/                    # CRISPR-Cas results (if --cctyper)
├── ncrna/                      # ncRNA predictions (if --ncrna)
├── genomad/                    # MGE predictions (if --genomad)
└── hoodini-viz/                # Visualization bundle
    ├── hoodini-viz.html        # Self-contained interactive viewer
    ├── tree.nwk                # Newick tree copy
    ├── tsv/                    # Human-readable tables
    │   ├── gff.gff
    │   ├── hoods.txt
    │   ├── protein_metadata.txt
    │   ├── tree_metadata.txt
    │   ├── protein_links.txt
    │   ├── nucleotide_links.txt
    │   ├── domains.txt         # (if --domains)
    │   └── ncrna_metadata.txt  # (if --ncrna)
    └── parquet/                # Parquet format for viewer
```

## Database Setup

```bash
hoodini download databases        # Download all required databases
hoodini download assembly_summary # NCBI assembly index only
hoodini download metacerberus     # Domain HMM profiles
```

## Configuration

Parameters can be set via CLI flags, TOML config file, or built-in defaults.  
Priority: CLI flags > config file > defaults.

```toml
# config.toml
[run]
num_threads = 8
win_mode = "win_genes"
win = 10
prot_links = true
tree_mode = "aai_tree"
```

```bash
hoodini run --config config.toml --input proteins.txt --output results
```

## Citation

If you use hoodini in your research, please cite:

> [Citation pending publication]

## License

See LICENSE file.
