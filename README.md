<h3 align="center">hoodini 🦉🎩: large-scale gene neighborhood analyses that feel like magic</h3>

---

<img src="docs/hoodini_logo_github.svg" alt="hoodini logo" width="220" align="left" style="border:0; margin:0 16px 8px 0; display:block;" />

<p align="justify">
   <b>hoodini</b> is a large-scale gene-centric comparative genomics toolkit that fetches public assemblies, extracts gene neighborhoods, runs pairwise protein and nucleotide comparisons, annotates neighborhoods with extra tools, and builds trees to help interpret genomic context.
</p>

---


## Introduction

This README documents installation, the command-line interface and all available options, produced output files, and a few examples to get you started.

## Installation

Before using `hoodini`, ensure you have Python and required system tools installed. The project provides an example Conda environment in `environment.yml`. It is advised to create the environment with mamba instead of conda due to beter conflict resolution. 

Suggested install steps (adapt channel list to your environment):

```bash
# create and activate environment (example)
mamba env create -f environment.yml
# install python-only extras (src-layout)
pip install -e .
# download hoodini external databases
hoodini download databases
```

## Quick start

Run the core pipeline using the CLI entrypoint. The package exposes a `hoodini` console script. The primary command is `run`:

```bash
hoodini run --input path/to/accessions.txt --output results_dir --prot-links --tree-mode aai_tree
```

Or provide a TOML config and override specific values on the command line:

```bash
hoodini run --config myparams.toml --inputsheet inputs.tsv --force
```

## How configuration is resolved

hoodini merges parameters from three places (lowest to highest precedence):

- built-in defaults (`hoodini/config/defaults.toml`)
- user-provided TOML via `--config`
- explicit CLI flags

This means CLI flags override the config file which overrides packaged defaults.

## Commands and options

This section lists the available top-level commands and all options discovered in the codebase.

Top-level command group: `hoodini`

Primary pipeline command: `hoodini run`

Options for `hoodini run` (all flags may also be present in a TOML config and are merged as described above):

- `--config <file>`
   - TOML config file to load parameters from (merged with defaults; CLI overrides this file).

- `--input <path>`
   - Path to a single-column input file (mutually exclusive with `--inputsheet`). Must exist.

- `--inputsheet <path>`
   - Path to a TSV input file with additional columns (mutually exclusive with `--input`). Must exist.

- `--output <folder>`
   - Output folder name. If omitted, defaults from packaged config will be used.

- `--max-concurrent-downloads <int>`
   - Maximum concurrent downloads for fetching files from NCBI.

- `--api-key` or `--api-key <key>`
   - NCBI API key (also read from `NCBI_API_KEY` environment variable where supported).

- `--num-threads <int>`
   - Number of worker threads to use where parallelism is supported.

- `--assembly-db <path>`
   - Path to an assembly database (precomputed index / parquet file).

- `--img-db <path>` and `--img-nuc <path>`
   - Path(s) to local IMG database files (protein and nucleotide variants) used for annotation.

- `--prot-links` (flag)
   - Run all-vs-all pairwise protein comparisons (produces `protein_links.txt`).

- `--nt-links` (flag)
   - Run pairwise nucleotide comparisons (produces `nucleotide_links.txt`).

- `--ani-mode <mode>`
   - ANI calculation mode used when building `ani_tree` (only used if `--tree-mode ani_tree`).

- `--nt-aln-mode <mode>`
   - Choice of nucleotide alignment engine for pairwise comparisons. Valid choices: `blastn`, `fastani`, `minimap2`, `intergenic_blastn`.

- `--blast <file>`
   - BLAST query file to use for annotation steps.

- `--cand-mode <mode>`
   - Mode for selecting IPG (representative protein) candidates during IPG parsing.

- `--clust-method <method>`
   - Protein clustering method (controls `cluster_proteins` behavior).

- `--win-mode <win_nts|win_genes>`
   - Window mode for neighborhoods: nucleotide window (`win_nts`) or gene window (`win_genes`).

- `--win <int>`
   - Window size (number of genes or nucleotides depending on `--win-mode`).

- `--min-win <int>`
   - Minimum window size on each side of the target.

- `--min-win-type <total|upstream|downstream|both>`
   - Type of min window constraint.

- `--tree-mode <mode>`
   - Controls tree construction behavior. Examples seen in code: `ani_tree`, `aai_tree`, or other supported modes.

- `--tree-file <path>`
   - Path to a precomputed tree file (Newick). If provided, hoodini may use this instead of rebuilding a tree.

- `--aai-mode <mode>`
   - Mode for AAI tree construction (e.g., `nj` or `hyper`). Note: certain modes may be rejected by the AAI pipeline.

- `--aai-subset-mode <mode>`
   - Subset mode used for selecting sequences for AAI tree construction (e.g., `target_region`, `target_prot`, `window`).

- `--padloc` (flag)
   - Run PADLOC to detect antiphage defense systems and merge annotations into the protein table.

- `--deffinder` (flag)
   - Run DefenseFinder for antiphage defense detection and merge results.

- `--ncrna` (flag)
   - Run Infernal to predict non-coding RNAs in neighborhoods.

- `--cctyper` (flag)
   - Run CCtyper for CRISPR-Cas system detection.

- `--genomad` (flag)
   - Run GenoMAD for mobile genetic element identification.

   - Identify anti-defense (ACR) genes.

   - Annotate proteins with PHROGs.

- `--sorfs` (flag)
   - Reannotate small open reading frames (sORFs) during assembly parsing/clustering.

- `--domains <comma-separated-list>`
   - Comma-separated list of MetaCerberus domain database names for domain annotation (validated by the CLI).

- `--min-prevalence <float>`
   - Minimum prevalence threshold for coloring/highlighting genes in downstream visualizations (used when computing proteome similarity).

- `--img` and `--img-metadata` <path>
   - Paths to IMG protein DB files and optional metadata for richer annotation.

- `--keep` (flag)
   - Keep intermediate temporary files (do not delete them at the end of the run).

- `--force` (flag)
   - Overwrite existing output folder if it already exists.

Notes about `run` behavior:

- The pipeline requires either `--input` or `--inputsheet`. If neither is provided, the CLI raises a usage error.
- The code merges defaults, an optional TOML config file, and CLI flags into a single runtime configuration object.

Subcommands: `hoodini download` group

The `download` group contains helpers to fetch and update local databases used by hoodini.

- `hoodini download assembly_summary`
   - Download or update the precomputed `assembly_summary.parquet` NCBI assembly database used for mapping accessions.

- `hoodini download metacerberus [dbs] [--force]`
   - Download MetaCerberus HMM/TSV database files. Calling with no argument or `all` lists or downloads as configured.

- `hoodini download type_dive`
   - Download and normalize DSMZ BacDive and PhageDive databases.

- `hoodini download contig_lengths [--api-key <key>] [--skip-assembly-summary]`
   - Download missing NCBI contig length records and refresh the precomputed list. `--api-key` overrides `NCBI_API_KEY` env var.

Utility commands: `hoodini utils`

- `hoodini utils nuc2asmlen <input_file> [--output <file>]`
   - Fetch assembly and contig length metadata for a list of nuccore/contig accessions; prints TSV to stdout or saves to `--output`.

- `hoodini utils prefetch_links <input_file> [--output <file>] [--kinds <list>]`
   - Generate a prefetched link table for a list of assemblies (useful to avoid repeated network requests). `--kinds` is a comma-separated list of file kinds.

## Output files written by the pipeline

The pipeline writes a small set of standardized output files to the `--output` folder. The following files are produced or placeholder files are created so downstream tools can rely on them:

- `gff.gff` — combined GFF for parsed assemblies.
- `hoods.txt` — hoods table with `hood_id`, `seqid`, `start`, `end`, `align_gene` for each neighborhood.
- `protein_metadata.txt` — protein table containing `gene_id`, `cluster`, `product`, and merged annotations.
- `tree_metadata.txt` — per-leaf metadata for the tree (leaf ids, etc.).
- `tree.nwk` — Newick-formatted tree string when a tree is produced or provided.
- `nucleotide_links.txt` — pairwise nucleotide alignment links (placeholder header-only file is written if nucleotide links are not produced).
- `protein_links.txt` — pairwise protein links (placeholder header-only file is written if protein comparisons are not run).

Extra annotations (PADLOC, DefenseFinder, PHROGs, etc.) are merged into `protein_metadata.txt` when requested.

### Visualization bundle (`hoodini-viz/`)
When the pipeline renders the interactive viewer, files are organized under `results/hoodini-viz/`:

- Root: `tree.nwk` and `hoodini-viz.html`.
- `tsv/`: `gff.gff`, `hoods.txt`, `protein_metadata.txt`, `tree_metadata.txt`, `nucleotide_links.txt`, `protein_links.txt`, and optional domain tables when requested.
- `parquet/`: columnar equivalents used by the viewer (`gff.parquet`, `hoods.parquet`, `protein_metadata.parquet`, `tree_metadata.parquet`, `nucleotide_links.parquet`, `protein_links.parquet`, and optional `domains.parquet`, `domains_metadata.parquet`).

## Examples

Basic run (protein links + AAI tree):

```bash
hoodini run --input accessions.txt --output hood_results --prot-links --tree-mode aai_tree --num-threads 8
```

Run only neighborhood extraction and annotation (no pairwise comparisons):

```bash
hoodini run --input accessions.txt --output hood_neigh --num-threads 4 --padloc --ncrna --cctyper
```

Use a config file and override one value:

```bash
hoodini run --config myparams.toml --inputsheet inputs.tsv --force
```

Download helper:

```bash
hoodini download assembly_summary
hoodini download metacerberus HMMs --force
```

Utility example:

```bash
hoodini utils nuc2asmlen contigs.txt --output contig_lengths.tsv
hoodini utils prefetch_links assemblies.txt --kinds gbff,gff,fna
```

## Troubleshooting & notes

- If you see unexpected errors while downloading from NCBI, ensure your `NCBI_API_KEY` is set or pass `--api-key`.
- Network-heavy operations benefit from `--max-concurrent-downloads` tuning and adequate `--num-threads`.
- For reproducible runs, pin the versions in `environment.yml` or use the provided `pyproject.toml`/packaging metadata.

## Contributing

See `TODO.md` and project tests for guidance on contributing. The core CLI lives in `hoodini/cli.py` and additional tools are in `hoodini/extra_tools/` and `hoodini/download/`.

## License

See `hoodini.egg-info/PKG-INFO` for packaged metadata and licensing information.
