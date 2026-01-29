# CLI Reference

This page documents the Hoodini CLI as implemented in [hoodini/src/hoodini/cli.py](hoodini/src/hoodini/cli.py) and the pipeline orchestration in [hoodini/src/hoodini/pipeline/runner.py](hoodini/src/hoodini/pipeline/runner.py).

## Commands

```
hoodini run
hoodini download
hoodini utils
```

## Inputs

### --input

Accepts one of the following:

- A literal protein accession or FASTA string (for a single query)
- A text file with one accession per line (must be multi‑line and single‑column; no commas or tabs)

If you use a literal input, Hoodini will perform a remote BLAST step to expand the search set using --remote-evalue and --remote-max-targets.

### --inputsheet

TSV with required columns:

- nucleotide_id
- protein_id
- gff_path
- fna_path
- faa_path

Additional optional columns that will be carried forward:

- uniprot_id, gbf_path, taxid, assembly_id, input_type, premade, failed

## Config file

You can supply a TOML config file using --config. Defaults are defined in [hoodini/src/hoodini/config/defaults.toml](hoodini/src/hoodini/config/defaults.toml) under sections general, window, tree, aai, ani, clustering, annotations, pairwise, and paths. CLI flags override config values.

## Pipeline stages and outputs

The CLI runs the following stages in order:

1. Initialization
	- Creates the output directory and reads inputs.
2. IPG parsing
	- Resolves IPG relationships for candidate selection.
3. Assembly parsing and neighborhood extraction
	- Downloads assemblies or uses local files, then extracts neighborhoods.
4. Protein comparisons (optional; required for AAI trees)
	- All‑vs‑all protein similarity table.
5. Nucleotide comparisons (optional; required for ANI trees)
	- ANI matrix and nucleotide links.
6. Protein clustering
	- Clusters neighborhood proteins for comparative visualization.
7. Proteome similarity (AAI)
	- AAI matrix for tree building.
8. Tree construction
	- Builds Newick tree for visualization.
9. Extra annotations (optional)
	- PADLOC, DefenseFinder, CCTyper, geNomad, ncRNA, eggNOG‑mapper, domains.
10. Visualization outputs
	- Parquet/TSV tables and a standalone HTML viewer.

For a full breakdown of output files, see [outputs](outputs).

## hoodini run options

### Output and config

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| --config | path | none | TOML config file |
| --output | path | results | Output directory |
| --force | flag | false | Overwrite existing output |
| --keep | flag | false | Keep intermediate files |

### Performance

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| --num-threads | int | 10 | Number of threads |
| --max-concurrent-downloads | int | 8 | Parallel NCBI downloads |
| --api-key | str | empty | NCBI API key (or NCBI_API_KEY env var) |

### Data sources

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| --assembly-folder | path | empty | Use local assemblies instead of downloading |
| --blast | path | empty | BLAST query file to use as region annotations |

### Neighborhood window

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| --win-mode | str | win_nts | Window mode: win_nts or win_genes |
| --win | int | 20000 | Window size (nucleotides or genes) |
| --min-win | int | 2000 | Minimum window per side |
| --min-win-type | str | both | total, upstream, downstream, or both |
| --sorfs | flag | false | Re‑annotate small ORFs |

### Clustering

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| --cand-mode | str | best_id | Candidate selection mode for IPG records (see details below) |
| --clust-method | str | diamond_deepclust | Protein clustering method (diamond_deepclust, deepmmseqs, jackhmmer, blastp) |

### Pairwise comparisons

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| --prot-links | flag | false | Compute protein similarity links |
| --nt-links | flag | false | Compute nucleotide links |
| --ani-mode | str | fastani | ANI calculation mode (skani or blastn) |
| --nt-aln-mode | str | blastn | Nucleotide alignment for links (blastn, fastani, minimap2, intergenic_blastn) |
| --min-pident | float | 30.0 | Minimum percent identity for AAI/wGRR |

### Tree construction

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| --tree-mode | str | fast_ml | Tree mode (see below) |
| --tree-file | path | target_prots.nwk | Input Newick tree for use_input_tree |
| --aai-mode | str | wgrr | AAI mode (wgrr or aai). hyper is rejected for AAI trees |
| --aai-subset-mode | str | target_region | Subset for AAI tree (target_prot, target_region, window) |

#### Candidate selection modes (cand_mode)

- any_ipg: keep any IPG entry; if RefSeq and GenBank versions overlap, prefer RefSeq for identical coordinates.
- best_ipg: pick a single representative per og_index using assembly level and edge proximity.
- best_id: same as best_ipg, but requires that the protein_id matches the original query when possible.
- one_id: keep the first IPG record per og_index.
- same_id: require match to the query protein_id and deduplicate by nucleotide_id and coordinates.

#### Tree modes (detailed)

- taxonomy: builds a tree from NCBI taxonomy distances using single‑linkage clustering.
- fast_nj: uses FAMSA to compute a distance matrix and DecentTree to build a neighbor‑joining or UPGMA tree.
- fast_ml: builds a protein alignment with FAMSA and infers a tree with VeryFastTree (removing rare columns at 0.10).
- aai_tree: builds a tree from AAI or wGRR pairwise tables using DecentTree; missing distances are filled with max+2std.
- ani_tree: builds a tree from ANI pairwise tables using DecentTree; missing distances are filled with max+2std.
- use_input_tree: loads the Newick tree from the path in --tree-file.
- foldmason_tree: maps target proteins to UniProt, downloads AlphaFold structures, builds an MSA with foldmason, and infers a tree with VeryFastTree. Falls back to fast_ml if mapping fails.
- neigh_similarity_tree: Jaccard distance on presence/absence of protein clusters across neighborhoods.
- neigh_phylo_tree: weighted neighborhood similarity using relative gene positions (cosine distance); falls back to neigh_similarity_tree if positions are missing.

### Remote BLAST (single‑query expansion)

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| --remote-evalue | float | 1e-5 | E‑value for remote BLAST |
| --remote-max-targets | int | 100 | Max hits for remote BLAST |

### Annotations

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| --padloc | flag | false | Run PADLOC |
| --deffinder | flag | false | Run DefenseFinder |
| --cctyper | flag | false | Run CCTyper |
| --genomad | flag | false | Run geNomad |
| --ncrna | flag | false | Infernal ncRNA prediction |
| --sorfs | flag | false | Reannotate small ORFs |
| --domains | str | empty | Comma‑separated MetaCerberus domains |
| --emapper | flag | false | Run eggNOG‑mapper |

### Logging

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| --quiet | flag | false | Silence non‑error output |
| --debug | flag | false | Verbose debug logging |

## hoodini download

Download resources used by Hoodini.

### Subcommands

- hoodini download assembly_summary
- hoodini download metacerberus [all|comma‑separated]
- hoodini download type_dive
- hoodini download contig_lengths [--api-key] [--skip-assembly-summary]
- hoodini download databases [--force] [--skip-padloc] [--skip-deffinder] [--skip-genomad] [--skip-emapper] [--skip-parquet] [--skip-contig-lengths] [--threads]

## hoodini utils

Utility commands for metadata helpers.

### Subcommands

- hoodini utils nuc2asmlen --output out.tsv input.tsv
- hoodini utils prefetch_links --output out.tsv --kinds gbff,gff,fna,faa,sequence_report input.tsv

## Defaults

Defaults come from [hoodini/src/hoodini/config/defaults.toml](hoodini/src/hoodini/config/defaults.toml). CLI flags override config values, which override defaults.
