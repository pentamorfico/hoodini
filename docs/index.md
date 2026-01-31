# Quick Start

This guide shows minimal runs and then summarizes the pipeline so you know what happens after you launch a job.

## 1) Single protein accession

```bash
hoodini run --input WP_012345678.1 --output results
```

Hoodini will download assemblies, extract neighborhoods, and generate a viewer in the output directory.

## 2) Batch input

Create a file with one accession per line:

```bash
cat > proteins.txt << EOF
WP_012345678.1
WP_087654321.1
WP_111222333.1
EOF

hoodini run --input proteins.txt --output results
```

## 3) Add protein links and a tree

```bash
hoodini run --input proteins.txt --output results \
  --prot-links \
  --tree-mode aai_tree
```

## 4) Add annotations

```bash
hoodini run --input proteins.txt --output results \
  --padloc --deffinder --cctyper --genomad
```

## 5) Provide a TSV with metadata

Use a TSV with the required columns and local file paths:

```bash
hoodini run --inputsheet my_samples.tsv --output results
```

## Pipeline summary (what happens under the hood)

1. Initialize inputs and output directory.
2. Resolve IPG records and select candidates.
3. Download assemblies and extract neighborhoods.
4. Optional protein comparisons and nucleotide comparisons.
5. Protein clustering and (optional) AAI/ANI tree construction.
6. Optional annotation tools (PADLOC, DefenseFinder, CCTyper, geNomad, ncRNA, eggNOG‑mapper).
7. Write visualization bundle and tabular outputs.

## Tips

- Add an NCBI API key to increase rate limits:
  ```bash
  export NCBI_API_KEY="your_key"
  ```
- Use --assembly-folder to reuse local GenBank/GFF files.
- If you only want neighborhood extraction, skip --prot-links and --tree-mode.

Next: read [CLI Reference](cli-reference) for full options and pipeline details.
