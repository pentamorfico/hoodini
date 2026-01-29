---
title: Outputs
description: Understanding Hoodini output files and structure
---

import { FileTree, Tabs, Callout } from 'nextra/components'

# Outputs

This page describes the files and folders produced by Hoodini and how to interpret them.

## Output Directory Structure

<FileTree>
  <FileTree.Folder name="output_dir" defaultOpen>
    <FileTree.Folder name="hoodini-viz" defaultOpen>
      <FileTree.Folder name="parquet" defaultOpen>
        <FileTree.File name="gff.parquet" />
        <FileTree.File name="hoods.parquet" />
        <FileTree.File name="protein_metadata.parquet" />
        <FileTree.File name="tree_metadata.parquet" />
        <FileTree.File name="domains.parquet" />
        <FileTree.File name="nucleotide_links.parquet" />
        <FileTree.File name="protein_links.parquet" />
      </FileTree.Folder>
      <FileTree.Folder name="tsv">
        <FileTree.File name="gff.gff" />
        <FileTree.File name="hoods.txt" />
        <FileTree.File name="protein_metadata.txt" />
        <FileTree.File name="tree_metadata.txt" />
      </FileTree.Folder>
      <FileTree.File name="tree.nwk" />
      <FileTree.File name="hoodini-viz.html" />
    </FileTree.Folder>
    <FileTree.Folder name="assembly_folder">
      <FileTree.File name="GCA_000001234.1/" />
      <FileTree.File name="GCA_000005678.1/" />
    </FileTree.Folder>
    <FileTree.Folder name="neighborhood">
      <FileTree.File name="neighborhoods.fasta" />
    </FileTree.Folder>
    <FileTree.Folder name="cctyper">
      <FileTree.File name="crispr_arrays.tsv" />
      <FileTree.File name="cas_operons.tsv" />
    </FileTree.Folder>
    <FileTree.Folder name="ncrna">
      <FileTree.File name="infernal_results.tblout" />
    </FileTree.Folder>
    <FileTree.Folder name="genomad">
      <FileTree.File name="virus_summary.tsv" />
      <FileTree.File name="plasmid_summary.tsv" />
    </FileTree.Folder>
    <FileTree.File name="assembly_list.txt" />
    <FileTree.File name="all_neigh.tsv" />
    <FileTree.File name="records.csv" />
    <FileTree.File name="tree.nwk" />
    <FileTree.File name="pairwise_aa.tsv" />
    <FileTree.File name="nt_links.tsv" />
    <FileTree.File name="aai_matrix.tsv" />
    <FileTree.File name="ani_matrix.tsv" />
    <FileTree.File name="domains.tsv" />
  </FileTree.Folder>
</FileTree>

---

## Core Outputs

These files are **always** produced:

| Path | Description |
|------|-------------|
| `assembly_list.txt` | Assembly accessions selected for processing |
| `all_neigh.tsv` | Neighborhood coordinates and identifiers for each target |
| `records.csv` | Input records enriched with taxonomy and metadata |
| `neighborhood/neighborhoods.fasta` | Extracted neighborhood nucleotide sequences |

<Callout type="info">
  The `records.csv` file contains your original input columns plus all metadata added during the pipeline, including taxonomy, assembly info, and any custom columns you provided.
</Callout>

---

## Assembly Folder

Downloaded genome files (when not using local assemblies):

<FileTree>
  <FileTree.Folder name="assembly_folder" defaultOpen>
    <FileTree.Folder name="GCA_000001234.1" defaultOpen>
      <FileTree.File name="GCA_000001234.1_genomic.gbff.gz" />
      <FileTree.File name="GCA_000001234.1_genomic.fna.gz" />
      <FileTree.File name="GCA_000001234.1_genomic.gff.gz" />
      <FileTree.File name="GCA_000001234.1_protein.faa.gz" />
    </FileTree.Folder>
    <FileTree.Folder name="GCA_000005678.1">
      <FileTree.File name="..." />
    </FileTree.Folder>
  </FileTree.Folder>
</FileTree>

---

## Comparative Analysis Outputs

These files are produced depending on your configuration:

<Tabs items={['Protein Links', 'Nucleotide Links', 'Distance Matrices']}>
  <Tabs.Tab>
    **`pairwise_aa.tsv`** - Produced with `--prot-links` or `aai` tree mode

    Contains protein-protein similarity scores:

    | Column | Description |
    |--------|-------------|
    | `query_id` | Source protein ID |
    | `target_id` | Target protein ID |
    | `pident` | Percent identity |
    | `evalue` | E-value |
    | `bitscore` | Bit score |
  </Tabs.Tab>
  <Tabs.Tab>
    **`nt_links.tsv`** - Produced with `--nt-links`

    Contains nucleotide alignments for synteny visualization:

    | Column | Description |
    |--------|-------------|
    | `query_hood` | Source neighborhood ID |
    | `target_hood` | Target neighborhood ID |
    | `query_start` | Start position in query |
    | `query_end` | End position in query |
    | `target_start` | Start position in target |
    | `target_end` | End position in target |
    | `identity` | Alignment identity |
  </Tabs.Tab>
  <Tabs.Tab>
    **Distance matrices** - Produced with `aai` or `ani` tree modes

    - `aai_matrix.tsv` - Average Amino acid Identity matrix
    - `ani_matrix.tsv` - Average Nucleotide Identity matrix

    Square matrices with pairwise distances between all neighborhoods.
  </Tabs.Tab>
</Tabs>

---

## Annotation Outputs

<Tabs items={['Domains', 'CRISPR-Cas', 'ncRNA', 'Mobile Elements']}>
  <Tabs.Tab>
    **`domains.tsv`** - Produced with `--domains`

    <FileTree>
      <FileTree.Folder name="output_dir">
        <FileTree.File name="domains.tsv" />
      </FileTree.Folder>
    </FileTree>

    | Column | Description |
    |--------|-------------|
    | `protein_id` | Protein accession |
    | `domain` | Domain name |
    | `database` | Source database (Pfam, TIGRfam, COG) |
    | `start` | Domain start position |
    | `end` | Domain end position |
    | `evalue` | E-value |
  </Tabs.Tab>
  <Tabs.Tab>
    **`cctyper/`** - Produced with `--cctyper`

    <FileTree>
      <FileTree.Folder name="cctyper" defaultOpen>
        <FileTree.File name="crispr_arrays.tsv" />
        <FileTree.File name="cas_operons.tsv" />
        <FileTree.File name="spacers.fasta" />
        <FileTree.File name="repeats.fasta" />
      </FileTree.Folder>
    </FileTree>

    Contains CRISPR array predictions and Cas protein classifications.
  </Tabs.Tab>
  <Tabs.Tab>
    **`ncrna/`** - Produced with `--ncrna`

    <FileTree>
      <FileTree.Folder name="ncrna" defaultOpen>
        <FileTree.File name="infernal_results.tblout" />
        <FileTree.File name="ncrna_summary.tsv" />
      </FileTree.Folder>
    </FileTree>

    Non-coding RNA predictions from Infernal (Rfam).
  </Tabs.Tab>
  <Tabs.Tab>
    **`genomad/`** - Produced with `--genomad`

    <FileTree>
      <FileTree.Folder name="genomad" defaultOpen>
        <FileTree.File name="virus_summary.tsv" />
        <FileTree.File name="plasmid_summary.tsv" />
        <FileTree.File name="provirus_summary.tsv" />
      </FileTree.Folder>
    </FileTree>

    Mobile genetic element predictions (viruses, plasmids, proviruses).
  </Tabs.Tab>
</Tabs>

---

## Visualization Bundle

The `hoodini-viz/` folder contains everything needed for the interactive viewer:

<FileTree>
  <FileTree.Folder name="hoodini-viz" defaultOpen>
    <FileTree.Folder name="parquet" defaultOpen>
      <FileTree.File name="gff.parquet" />
      <FileTree.File name="hoods.parquet" />
      <FileTree.File name="protein_metadata.parquet" />
      <FileTree.File name="tree_metadata.parquet" />
      <FileTree.File name="domains.parquet" />
      <FileTree.File name="nucleotide_links.parquet" />
      <FileTree.File name="protein_links.parquet" />
    </FileTree.Folder>
    <FileTree.Folder name="tsv">
      <FileTree.File name="gff.gff" />
      <FileTree.File name="hoods.txt" />
      <FileTree.File name="protein_metadata.txt" />
      <FileTree.File name="tree_metadata.txt" />
    </FileTree.Folder>
    <FileTree.File name="tree.nwk" />
    <FileTree.File name="hoodini-viz.html" />
  </FileTree.Folder>
</FileTree>

### Parquet Files

Efficient binary format used by the viewer:

| File | Contents |
|------|----------|
| `gff.parquet` | Gene annotations for all neighborhoods |
| `hoods.parquet` | Neighborhood metadata and sequences |
| `protein_metadata.parquet` | Protein info including clusters and domains |
| `tree_metadata.parquet` | Tree leaf metadata (taxonomy, custom columns) |
| `domains.parquet` | Domain annotations (if `--domains` used) |
| `nucleotide_links.parquet` | Synteny links (if `--nt-links` used) |
| `protein_links.parquet` | Protein similarity links (if `--prot-links` used) |

### TSV Files

Human-readable versions for inspection:

| File | Contents |
|------|----------|
| `gff.gff` | Standard GFF3 format annotations |
| `hoods.txt` | Neighborhood coordinates |
| `protein_metadata.txt` | Protein annotations |
| `tree_metadata.txt` | Leaf metadata |

### Viewer HTML

**`hoodini-viz.html`** - Self-contained HTML file that loads the parquet files and displays the interactive visualization.

<Callout type="info" emoji="💡">
  **Tip:** You can share the entire `hoodini-viz/` folder. The HTML file loads data from the parquet files, so keep them together!
</Callout>

---

## Tree File

**`tree.nwk`** - Newick format phylogenetic tree

Produced when `tree_mode` is not `none`. The tree type depends on your setting:

- **taxonomy**: Tree based on NCBI taxonomy hierarchy
- **aai**: Tree based on Average Amino acid Identity distances
- **ani**: Tree based on Average Nucleotide Identity distances

The leaf names correspond to neighborhood IDs (`uid` column in other files).
