# Hoodini CLI

Hoodini is a gene‑centric comparative genomics toolkit for microbial genomes. It automates sequence retrieval, neighborhood extraction, optional comparative analyses, and interactive visualization.

This documentation is written for biologists and bioinformaticians working in microbial genomics. It explains how to run the CLI, how to structure inputs, and what each pipeline stage does.

## Inputs at a glance

You can run Hoodini with either:

- --input: a literal protein accession or FASTA string, or a multi‑line text file with one ID per line
- --inputsheet: a TSV with required columns nucleotide_id, protein_id, gff_path, fna_path, faa_path

## Pipeline at a glance

Hoodini performs the following stages, with optional steps enabled by flags:

1. Initialize inputs and output folder.
2. Parse IPG relationships for candidate selection.
3. Download assemblies and extract neighborhoods.
4. Protein comparisons (optional; required for AAI trees).
5. Nucleotide comparisons (optional; required for ANI trees).
6. Protein clustering.
7. Proteome similarity for AAI trees.
8. Tree construction.
9. Extra annotations (domains, PADLOC, DefenseFinder, CCTyper, geNomad, ncRNA, eggNOG‑mapper).
10. Visualization bundle and tables.

See [CLI Reference](cli-reference) for stage‑by‑stage details and outputs.

## When to use Hoodini

- Compare neighborhoods across many genomes for a protein of interest.
- Explore synteny, operons, and genomic islands at scale.
- Add defense and mobile element annotations.
- Generate publication‑ready figures.

## Next steps

- [Installation](installation)
- [Quick Start](quickstart)
- [CLI Reference](cli-reference)
si