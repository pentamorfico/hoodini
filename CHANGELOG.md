# Changelog

User-facing changes for each Hoodini release are recorded here. Entries under
Unreleased describe changes implemented on the update branch, not a published release.

## Unreleased

### Fixed

- Local GFF/FAA neighborhoods now include features only from the selected contig,
  preventing unrelated proteins from entering downstream annotations. Protein
  queries can infer a unique contig; ambiguous inputs report a selection error.
  The matching FNA sequence is used when the contig is inferred. (#83)
- Mixed GBFF and GFF/FAA inputs no longer fail because numeric GFF score/phase
  values conflict with GenBank strings. GFF comments are skipped during parsing.
  (#84)
- Contig metadata lookups tolerate incrementally updated Parquet partitions with
  different optional accession columns, including GenBank-only and RefSeq-only
  datasets. Both DuckDB queries and the Polars fallback normalize their schemas.
  New partitions preserve optional metadata appearing after the first row and
  consistently type the core lookup columns. (#86)
- NCBI lineages using `domain` now populate the legacy `superkingdom` output field,
  avoiding false unclassified labels for Bacteria and Archaea. An explicit
  `superkingdom` takes priority when both ranks are present. (#81)

See the [release work log](docs/releases/next-release.md) for validation evidence,
remaining release checks and the full issue inventory.
