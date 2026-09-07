"""Lazy, schema-stable access to incrementally updated NCBI contig tables."""

from __future__ import annotations

from glob import glob
from typing import Any

import polars as pl

CONTIG_SCHEMA = {
    "genbankAccession": pl.Utf8,
    "refseqAccession": pl.Utf8,
    "assemblyAccession": pl.Utf8,
    "length": pl.Int64,
}
REQUIRED_COLUMNS = {"assemblyAccession", "length"}


def _validate_columns(columns: set[str]) -> None:
    missing = REQUIRED_COLUMNS - columns
    if missing:
        raise ValueError(f"Contig metadata is missing required columns: {sorted(missing)}")


def register_contig_table(connection: Any, parquet_path: str) -> None:
    """Expose a projected DuckDB view without materializing the contig database.

    GenBank-only and RefSeq-only partitions are valid. Optional accession columns
    must also exist when absent from every partition, which union_by_name alone
    cannot guarantee. The dataset must contain assembly and length columns.
    """
    relation = connection.read_parquet(parquet_path, union_by_name=True)
    columns = set(relation.columns)
    _validate_columns(columns)
    expressions = []
    for name in CONTIG_SCHEMA:
        value = f'"{name}"' if name in columns else "NULL"
        dtype = "BIGINT" if name == "length" else "VARCHAR"
        expressions.append(f'CAST({value} AS {dtype}) AS "{name}"')
    relation.project(", ".join(expressions)).create_view("hoodini_contigs", replace=True)


def scan_contig_table(parquet_path: str) -> pl.LazyFrame:
    """Project and normalize each partition before a lazy Polars concatenation.

    This supports older Polars versions without relying on missing_columns or
    extra_columns scan options. Only the four lookup columns enter the query.
    """
    paths = sorted(glob(parquet_path))
    if not paths:
        raise FileNotFoundError(f"No contig metadata files match {parquet_path}")
    partitions = []
    for path in paths:
        frame = pl.scan_parquet(path)
        columns = set(pl.read_parquet_schema(path))
        _validate_columns(columns)
        partitions.append(
            frame.select(
                [
                    (pl.col(name) if name in columns else pl.lit(None)).cast(dtype).alias(name)
                    for name, dtype in CONTIG_SCHEMA.items()
                ]
            )
        )
    return pl.concat(partitions, how="vertical")
