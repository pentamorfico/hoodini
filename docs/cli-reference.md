---
title: CLI Reference
description: Complete reference for all Hoodini CLI commands and options
---

import { Tabs, Callout, Cards } from 'nextra/components'

# CLI Reference

Hoodini provides a powerful command-line interface for comparative genomics analysis.

## Commands Overview

<Cards>
  <Cards.Card title="hoodini run" href="#hoodini-run" />
  <Cards.Card title="hoodini download" href="#hoodini-download" />
  <Cards.Card title="hoodini utils" href="#hoodini-utils" />
</Cards>

---

## hoodini run

The main pipeline command that orchestrates the complete analysis workflow.

```bash
hoodini run --input <accessions.txt> --output results/
```

### Input Options

<Tabs items={['--input', '--inputsheet']}>
  <Tabs.Tab>
    **Single file or literal input**
    
    ```bash
    # File with one accession per line
    hoodini run --input accessions.txt
    
    # Literal protein accession
    hoodini run --input "WP_012345678.1"
    
    # Nucleotide with coordinates
    hoodini run --input "NC_000913.3:1000-5000"
    ```
    
    <Callout type="info">
      When using a literal input, Hoodini performs a remote BLAST to expand the search set using `--remote-evalue` and `--remote-max-targets`.
    </Callout>
  </Tabs.Tab>
  
  <Tabs.Tab>
    **TSV with metadata, local files, or custom columns**
    
    ```bash
    hoodini run --inputsheet metadata.tsv
    ```
    
    **Minimum required** (at least one of):
    - `nucleotide_id` — NCBI nucleotide accession
    - `protein_id` — NCBI protein accession  
    - `uniprot_id` — UniProt accession
    
    **For local files** (GFF + FAA or GenBank):
    - `gff_path`, `faa_path` — local annotation files
    - `gbf_path` — alternative: GenBank format
    
    **Coordinates** (optional):
    - `start`, `end`, `strand` — analyze specific region
    
    **Custom columns**: Any additional columns are preserved and appear in final outputs (`hoods.txt`, `tree_metadata.txt`).
    
    <Callout type="info" emoji="📖">
      See [Input Formats](/input-formats) for complete documentation on inputsheet structure and custom columns.
    </Callout>
  </Tabs.Tab>
</Tabs>

### Output & Config

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--config` | path | — | TOML config file |
| `--output` | path | `results` | Output directory |
| `--force` | flag | `false` | Overwrite existing output |
| `--keep` | flag | `false` | Keep intermediate files |

### Performance

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--num-threads` | int | `10` | Number of threads |
| `--max-concurrent-downloads` | int | `8` | Parallel NCBI downloads |
| `--api-key` | str | — | NCBI API key (or `NCBI_API_KEY` env var) |

### Data Sources

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--assembly-folder` | path | — | Use local assemblies instead of downloading |

### Neighborhood Window

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--win-mode` | str | `win_nts` | Window mode: `win_nts` or `win_genes` |
| `--win` | int | `20000` | Window size (nucleotides or genes) |
| `--min-win` | int | `2000` | Minimum window per side |
| `--min-win-type` | str | `both` | `total`, `upstream`, `downstream`, or `both` |

### Clustering

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--cand-mode` | str | `best_id` | Candidate selection mode (see below) |
| `--clust-method` | str | `diamond_deepclust` | Clustering method |

<Callout type="default" emoji="📋">
**Candidate selection modes (`--cand-mode`)**

- **`best_id`** *(default)*: Pick single best representative per input homolog, prioritizing assembly quality and edge proximity. The protein ID must match the original query when possible.
- **`best_ipg`**: Same as `best_id`, but allows different protein IDs (e.g., GenBank vs RefSeq versions of the same protein).
- **`same_id`**: Keep all IPG records that share the same protein ID as the query. Especially relevant for non-redundant proteins (`WP_`, `YP_`, `NP_`) present in multiple assemblies.
- **`any_ipg`**: Keep all identical proteins from IPG regardless of ID. Can result in massive expansion if the protein exists in thousands of assemblies.
- **`one_id`**: Keep first IPG record per input homolog, regardless of assembly quality (order from NCBI IPG).
</Callout>

<Callout type="warning">
  **Warning:** Using `any_ipg` or `same_id` can dramatically increase the number of neighborhoods if your query protein is highly conserved across many assemblies.
</Callout>

### Pairwise Comparisons

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--prot-links` | flag | `false` | Compute protein similarity links |
| `--nt-links` | flag | `false` | Compute nucleotide links |
| `--ani-mode` | str | `fastani` | ANI calculation: `skani` or `blastn` |
| `--nt-aln-mode` | str | `blastn` | Nucleotide alignment: `blastn`, `fastani`, `minimap2`, `intergenic_blastn` |
| `--min-pident` | float | `30.0` | Minimum percent identity for AAI/wGRR |

### Tree Construction

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--tree-mode` | str | `fast_ml` | Tree building method (see below) |
| `--tree-file` | path | — | Input Newick tree for `use_input_tree` |
| `--aai-mode` | str | `wgrr` | AAI mode: `wgrr` or `aai` |
| `--aai-subset-mode` | str | `target_region` | Subset for AAI tree: `target_prot`, `target_region`, `window` |

<Callout type="default" emoji="🌳">
**Tree modes (`--tree-mode`)**

| Mode | Description |
|------|-------------|
| `taxonomy` | NCBI taxonomy distances with single-linkage clustering |
| `fast_nj` | FAMSA distance matrix → DecentTree NJ/UPGMA |
| `fast_ml` | FAMSA alignment → VeryFastTree (default) |
| `aai_tree` | AAI/wGRR pairwise distances → DecentTree |
| `ani_tree` | ANI pairwise distances → DecentTree |
| `use_input_tree` | Load from `--tree-file` |
| `foldmason_tree` | AlphaFold structures → foldmason MSA → VeryFastTree |
| `neigh_similarity_tree` | Jaccard distance on protein cluster presence/absence |
| `neigh_phylo_tree` | Weighted neighborhood similarity using gene positions |
</Callout>

### Remote BLAST

For single-query expansion:

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--remote-evalue` | float | `1e-5` | E-value for remote BLAST |
| `--remote-max-targets` | int | `100` | Max hits for remote BLAST |

### Annotations

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--padloc` | flag | `false` | Run PADLOC defense system detection |
| `--deffinder` | flag | `false` | Run DefenseFinder |
| `--cctyper` | flag | `false` | Run CCTyper CRISPR-Cas typing |
| `--genomad` | flag | `false` | Run geNomad virus/plasmid detection |
| `--ncrna` | str | — | ncRNA prediction: path to CM file OR comma-separated RFAM IDs (e.g., `RF00001,RF00234`) |
| `--domains` | str | — | Comma-separated MetaCerberus domains |
| `--emapper` | flag | `false` | Run eggNOG-mapper |
| `--blast` | path | — | FASTA file for BLAST search against neighborhood nucleotides (e.g., IS elements) |
| `--sorfs` | flag | `false` | Re-annotate small ORFs in extracted regions |

### Logging

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `--quiet` | flag | `false` | Silence non-error output |
| `--debug` | flag | `false` | Verbose debug logging |

---

## hoodini download

Download databases and resources used by Hoodini.

```bash
hoodini download databases --threads 8
```

### Subcommands

<Tabs items={['databases', 'assembly_summary', 'metacerberus', 'type_dive', 'contig_lengths']}>
  <Tabs.Tab>
    **Download all databases** (~35 GB total)
    
    Downloads all required databases for full functionality: PADLOC, DefenseFinder, geNomad, eggNOG-mapper, and supporting files.
    
    ```bash
    hoodini download databases [OPTIONS]
    ```
    
    | Option | Description |
    |--------|-------------|
    | `--force` | Re-download existing files |
    | `--threads` | Number of threads |
    | `--skip-padloc` | Skip PADLOC models |
    | `--skip-deffinder` | Skip DefenseFinder models |
    | `--skip-genomad` | Skip geNomad database |
    | `--skip-emapper` | Skip eggNOG-mapper data |
    | `--skip-parquet` | Skip parquet files |
    | `--skip-contig-lengths` | Skip contig length database |
  </Tabs.Tab>
  
  <Tabs.Tab>
    **Download NCBI assembly summary**
    
    Metadata table mapping accessions to assembly information. Used to resolve protein/nucleotide IDs to their source assemblies.
    
    ```bash
    hoodini download assembly_summary
    ```
  </Tabs.Tab>
  
  <Tabs.Tab>
    **Download MetaCerberus domain models**
    
    HMM profiles for protein domain annotation. Adds functional annotations to proteins in your neighborhoods (PFAM, TIGRFAM, etc.).
    
    ```bash
    # Download all domains
    hoodini download metacerberus all
    
    # Download specific domains
    hoodini download metacerberus PFAM,TIGRFAM
    ```
  </Tabs.Tab>
  
  <Tabs.Tab>
    **Download TypeDive database**
    
    Tables linking NCBI assemblies to BacDive and PhageDive strain collections. Useful to identify which neighborhoods come from strains available in culture collections.
    
    ```bash
    hoodini download type_dive
    ```
  </Tabs.Tab>
  
  <Tabs.Tab>
    **Download contig lengths database**
    
    Pre-computed contig lengths for NCBI assemblies. Speeds up neighborhood extraction by avoiding per-assembly lookups.
    
    ```bash
    hoodini download contig_lengths [--api-key KEY] [--skip-assembly-summary]
    ```
  </Tabs.Tab>
</Tabs>

---

## hoodini utils

Utility commands for metadata helpers.

```bash
hoodini utils nuc2asmlen --output out.tsv input.tsv
```

### Subcommands

<Tabs items={['nuc2asmlen', 'prefetch_links']}>
  <Tabs.Tab>
    **Convert nucleotide IDs to assembly lengths**
    
    ```bash
    hoodini utils nuc2asmlen --output out.tsv input.tsv
    ```
    
    Takes a TSV with nucleotide/contig accessions and adds assembly and contig length metadata.
  </Tabs.Tab>
  
  <Tabs.Tab>
    **Prefetch download links**
    
    ```bash
    hoodini utils prefetch_links --output out.tsv --kinds gbff,gff,fna,faa input.tsv
    ```
    
    Generates a table of `assembly_id`, `file_type`, and `link` for batch downloading.
  </Tabs.Tab>
</Tabs>

---

## Configuration File

You can supply a TOML config file using `--config`. CLI flags always override config values.

```bash
hoodini run --config my_config.toml --input accessions.txt
```

<Callout type="info">
  Default values are defined in `hoodini/config/defaults.toml` under sections: `general`, `window`, `tree`, `aai`, `ani`, `clustering`, `annotations`, `pairwise`, and `paths`.
</Callout>

See [outputs](/docs/hoodini/outputs) for details on output files and directory structure.
