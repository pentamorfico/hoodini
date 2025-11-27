## Preparing data

...

## Example basic run

The mandatory arguments for a basic run of Hoodini are an input file containing a list of accessions to be downloaded (`--input`) and an output directory where results will be stored (`--output`). The input file should contain one accession per line.

```bash
hoodini run --input example/accessions.txt --output example/results
```

??? folder "**Output**" 

	| Method   | Description                          |
	| :------- | :----------------------------------- |
	| `defaultGFF.gff`    | combined GFF for parsed assemblies.  |
	| `defaultBaselines.txt`    | baseline table with hood_id, seqid, start, end, align_gene for each neighborhood. |
	| `defaultProteinMetadata.txt` | protein table containing gene_id, cluster, product, and merged annotations. |
    | `defaultTreeMetadata.txt` | per-leaf metadata for the tree (leaf ids, etc.). |
    | `defaultNewick.txt` | Newick-formatted tree string when a tree is produced or provided. |
    | `defaultNucleotideLinks.txt` | pairwise nucleotide alignment links (placeholder header-only file is written if nucleotide links are not produced). |
    | `defaultProteinLinks.txt` | pairwise protein links (placeholder header-only file is written if protein comparisons are not run). |
