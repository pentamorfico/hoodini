---
title: Input Formats
description: Complete guide to Hoodini input formats, inputsheet structure, and custom columns
---

import { Tabs, Callout, FileTree } from 'nextra/components'

# Input Formats

Hoodini supports multiple input formats depending on your data source and use case. This guide covers all supported formats, required columns, and how to add custom metadata that propagates to the final outputs.

## Quick Reference

| Input Type | Command | Use Case |
|------------|---------|----------|
| Simple text file | `--input proteins.txt` | List of NCBI/UniProt IDs |
| Single ID | `--input WP_012345678.1` | Single query (triggers BLAST) |
| Inputsheet (TSV) | `--inputsheet samples.tsv` | Custom metadata, local files, regions |

---

## Simple Input File (`--input`)

A plain text file with one accession per line. Hoodini auto-detects the ID type.

```txt title="proteins.txt"
WP_000000001.1
WP_000000002.1
NP_414542.1
```

```txt title="nucleotides.txt"
NC_000913.3
NZ_CP012345.1
MZ501047.1
```

### Supported ID Formats

| Format | Example | Description |
|--------|---------|-------------|
| NCBI Protein | `WP_000000001.1`, `NP_414542.1` | RefSeq/GenBank protein IDs |
| NCBI Nucleotide | `NC_000913.3`, `NZ_CP012345.1` | RefSeq/GenBank contig/chromosome |
| UniProt | `P12345`, `Q9Y6K9` | Auto-converted to NCBI protein ID |
| Region format | `NC_000913.3:1000-5000` | Specific genomic coordinates |

<Callout type="info">
  **Region format**: Use `NucID:start-end` to analyze a specific genomic region. If `start > end`, the strand is set to `-`.
</Callout>

---

## Inputsheet Format (`--inputsheet`)

A **tab-separated file (TSV)** that allows you to:

- Specify genomic regions with coordinates
- Use local annotation files (GFF, FAA, GenBank)
- Add custom metadata columns that propagate to outputs
- Mix different input types in one analysis

### Minimum Required Columns

You must provide **at least one** of these ID columns:

| Column | Description |
|--------|-------------|
| `nucleotide_id` | NCBI nucleotide accession (e.g., `NC_000913.3`) |
| `protein_id` | NCBI protein accession (e.g., `WP_000000001.1`) |
| `uniprot_id` | UniProt accession (e.g., `P12345`) |

<Callout type="warning">
  The priority order is: `nucleotide_id` > `protein_id` > `uniprot_id`. If multiple columns have values, the higher priority one is used.
</Callout>

### Basic Example

```tsv title="basic_inputsheet.tsv"
nucleotide_id
NC_000913.3
MZ501047.1
MZ501048.1
```

### With Coordinates

```tsv title="regions_inputsheet.tsv"
nucleotide_id	start	end	strand
NC_000913.3	1000000	1050000	+
NC_000913.3	2000000	2050000	-
MZ501047.1			
```

<Callout type="info">
  Leave `start` and `end` empty to analyze the **full contig**. This is useful for phage genomes or plasmids where you want the entire sequence.
</Callout>

---

## Using Local Files

When you have your own genome annotations (not from NCBI), provide paths to local files.

### Option 1: GFF + FASTA (Recommended)

| Column | Required | Description |
|--------|:--------:|-------------|
| `nucleotide_id` | ✅ | Sequence ID (must match seqid in GFF) |
| `gff_path` | ✅ | Path to GFF3 annotation file |
| `faa_path` | ✅ | Path to protein FASTA file |
| `fna_path` | ❌ | Path to nucleotide FASTA (for NT analysis) |

```tsv title="local_gff_inputsheet.tsv"
nucleotide_id	gff_path	faa_path	fna_path
contig_001	/data/genome1/annotation.gff	/data/genome1/proteins.faa	/data/genome1/genome.fna
contig_002	/data/genome2/annotation.gff	/data/genome2/proteins.faa	/data/genome2/genome.fna
```

### Option 2: GenBank Format

| Column | Required | Description |
|--------|:--------:|-------------|
| `nucleotide_id` | ✅ | Sequence ID |
| `gbf_path` | ✅ | Path to GenBank file (.gbf, .gbk, .gb) |

```tsv title="local_genbank_inputsheet.tsv"
nucleotide_id	gbf_path
contig_001	/data/genome1.gbk
contig_002	/data/genome2.gbk
```

### File Format Requirements

<Callout type="warning" emoji="⚠️">
  **Important**: The `nucleotide_id` must match the sequence identifier in your files:
  - In GFF: the first column (seqid) of CDS features
  - In GenBank: the LOCUS name or ACCESSION
  - In FASTA headers: the sequence ID before the first space
</Callout>

**GFF3 Requirements:**
- Must have `CDS` features
- Each CDS must have an `ID=` attribute (or `locus_tag=`)
- The seqid (column 1) must match your `nucleotide_id`

**Protein FASTA Requirements:**
- Headers must contain the CDS ID from the GFF
- Example: `>gene_001 hypothetical protein` where `gene_001` is the `ID` in the GFF

---

## Custom Columns (Extra Metadata)

Any column in your inputsheet that is **not a reserved column** will be treated as **custom metadata** and automatically propagated to the final outputs.

### How It Works

1. Add any columns you want to your inputsheet
2. Hoodini preserves them through the pipeline
3. They appear in:
   - `hoods.txt` / `hoods.parquet` — neighborhood data
   - `tree_metadata.txt` / `tree_metadata.parquet` — tree leaf metadata

### Example with Custom Columns

```tsv title="samples_with_metadata.tsv"
nucleotide_id	sample_name	host	isolation_source	collection_year	experiment_id
MZ501047.1	Phage_Alpha	E. coli	Wastewater	2023	EXP001
MZ501048.1	Phage_Beta	S. enterica	Soil	2022	EXP002
MZ501049.1	Phage_Gamma	K. pneumoniae	Hospital	2024	EXP003
NC_000913.3:1000000-1050000	Region_A	E. coli K-12	Lab strain	2020	EXP004
```

### Output with Custom Columns

**hoods.txt:**
```tsv
hood_id	seqid	start	end	align_gene	sample_name	host	isolation_source	collection_year	experiment_id
0	MZ501047.1	1	45678	gene_001	Phage_Alpha	E. coli	Wastewater	2023	EXP001
1	MZ501048.1	1	43210	gene_042	Phage_Beta	S. enterica	Soil	2022	EXP002
```

**tree_metadata.txt:**
```tsv
leaf_id	og_index	superkingdom	phylum	...	sample_name	host	isolation_source	collection_year	experiment_id
0	0	Viruses	...	...	Phage_Alpha	E. coli	Wastewater	2023	EXP001
1	1	Viruses	...	...	Phage_Beta	S. enterica	Soil	2022	EXP002
```

<Callout type="info" emoji="💡">
  **Use case**: Add sample metadata, experimental conditions, or any annotation you want to visualize alongside your genomic neighborhoods in the HTML viewer.
</Callout>

---

## Reserved Columns

These columns have special meaning in the pipeline and **should not be used for custom data**:

<details>
<summary>Click to see all reserved columns</summary>

**Input identification:**
- `og_index`, `unique_id`, `protein_id`, `nucleotide_id`, `uniprot_id`, `input_type`

**File paths:**
- `gff_path`, `faa_path`, `fna_path`, `gbf_path`

**Coordinates:**
- `start`, `end`, `strand`

**Assembly/taxonomy:**
- `taxid`, `assembly_id`

**Status flags:**
- `failed`, `failed_reason`, `premade`, `is_full_contig`

**Query info (added by pipeline):**
- `query_protein_id`, `is_refseq_query`, `sequence_length`, `group`
- `species_taxid`, `organism_name`, `infraspecific_name`, `assembly_level`
- `nucleotide_id_no_prefix`

**DSMZ columns:**
- `dive_id`, `collection_id`, `dive_type`

</details>

---

## Complete Example

Here's a complete inputsheet combining multiple features:

```tsv title="complete_example.tsv"
nucleotide_id	start	end	gff_path	faa_path	sample_name	condition	replicate
NC_000913.3	1000000	1050000			Sample_A	treatment	1
NC_000913.3	2000000	2050000			Sample_B	treatment	2
MZ501047.1				Phage_X	control	1
contig_local	10000	50000	/data/local.gff	/data/local.faa	Local_Sample	experimental	1
```

This inputsheet:
- Analyzes two regions from `NC_000913.3` (downloaded from NCBI)
- Analyzes the full genome of phage `MZ501047.1` (downloaded from NCBI)
- Uses local files for `contig_local`
- Adds `sample_name`, `condition`, and `replicate` as custom metadata

---

## Tips and Best Practices

1. **Use TSV, not CSV**: Hoodini expects tab-separated values
2. **Empty cells**: Leave cells empty (not `NA` or `NULL`) for missing values
3. **Paths**: Can be absolute or relative to your working directory
4. **Consistent IDs**: Ensure `nucleotide_id` matches your GFF seqid exactly
5. **Column names**: Avoid spaces and special characters in custom column names
6. **Full contigs**: Leave `start`/`end` empty to analyze entire sequences (great for phages/plasmids)
