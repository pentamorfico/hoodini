"""Centralized table schemas and column helpers for Polars dataframes.

These schemas define the expected columns and dtypes for core tables. They should be
used at module boundaries to validate inputs/outputs and to make pandas→Polars
migration explicit.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass

import polars as pl


@dataclass(frozen=True)
class TableSchema:
    name: str
    required: Mapping[str, pl.DataType]
    optional: Mapping[str, pl.DataType] | None = None

    def ensure(self, df: pl.DataFrame, *, allow_extra: bool = True) -> pl.DataFrame:
        """Validate presence and cast columns to declared dtypes.

        - Ensures required columns exist.
        - Optionally casts to declared dtypes when possible.
        - If allow_extra is False, drops undeclared columns.
        """
        missing = [c for c in self.required if c not in df.columns]
        if missing:
            raise ValueError(f"{self.name}: missing required columns: {missing}")

        cast_exprs = []
        for col, dtype in self.required.items():
            cast_exprs.append(_cast_column(col, dtype))
        if self.optional:
            for col, dtype in self.optional.items():
                if col in df.columns:
                    cast_exprs.append(_cast_column(col, dtype))

        casted = df.with_columns(cast_exprs)
        if not allow_extra:
            keep = set(self.required) | set(self.optional or {})
            casted = casted.select([c for c in casted.columns if c in keep])
        return casted


TRUE_STRINGS = {"true", "1", "yes", "y", "t"}
FALSE_STRINGS = {"false", "0", "no", "n", "f"}


def _cast_column(col: str, dtype: pl.DataType) -> pl.Expr:
    """Robustly cast columns, with special handling for booleans."""

    if dtype == pl.Boolean:
        norm = pl.col(col).cast(pl.Utf8, strict=False).str.to_lowercase()
        return (
            pl.when(pl.col(col).is_null())
            .then(None)
            .otherwise(
                pl.when(norm.is_in(list(TRUE_STRINGS)))
                .then(True)
                .when(norm.is_in(list(FALSE_STRINGS)))
                .then(False)
                .otherwise(pl.col(col).cast(pl.Boolean, strict=False))
            )
            .alias(col)
        )

    return pl.col(col).cast(dtype, strict=False)


RECORDS = TableSchema(
    name="records",
    required={
        "og_index": pl.Int64,
        "input_type": pl.Utf8,
    },
    optional={
        "protein_id": pl.Utf8,
        "nucleotide_id": pl.Utf8,
        "uniprot_id": pl.Utf8,
        "failed": pl.Boolean,
        "failed_reason": pl.Utf8,
        "gff_path": pl.Utf8,
        "faa_path": pl.Utf8,
        "fna_path": pl.Utf8,
        "strand": pl.Utf8,
        "start": pl.Int64,
        "end": pl.Int64,
        "gbf_path": pl.Utf8,
        "taxid": pl.Int64,
        "assembly_id": pl.Utf8,
        "premade": pl.Boolean,
        "is_full_contig": pl.Boolean,
        # DSMZ BacDive/PhageDive columns
        "dive_id": pl.Utf8,
        "collection_id": pl.Utf8,
        "dive_type": pl.Utf8,
    },
)

PROTEINS = TableSchema(
    name="proteins",
    required={"id": pl.Utf8, "sequence": pl.Utf8},
    optional={
        "protein_id": pl.Utf8,
        "target_prot": pl.Utf8,
        "target_nuc": pl.Utf8,
        "fam_cluster": pl.Utf8,
        "product": pl.Utf8,
        "attributes": pl.Utf8,
        "unique_id": pl.Utf8,
    },
)

NEIGHBORHOODS = TableSchema(
    name="neighborhoods",
    required={
        "seqid": pl.Utf8,
        "start_win": pl.Int64,
        "end_win": pl.Int64,
        "sequence": pl.Utf8,
        "unique_id": pl.Utf8,
    },
    optional={
        "start_target": pl.Int64,
        "end_target": pl.Int64,
        "target_prot": pl.Utf8,
        "temp_seqid": pl.Utf8,
        "strand_win": pl.Utf8,
        "is_full_contig": pl.Boolean,
    },
)

GFF = TableSchema(
    name="gff",
    required={
        "seqid": pl.Utf8,
        "source": pl.Utf8,
        "type": pl.Utf8,
        "start": pl.Int64,
        "end": pl.Int64,
        "score": pl.Utf8,
        "strand": pl.Utf8,
        "phase": pl.Utf8,
        "attributes": pl.Utf8,
    },
    optional={"id": pl.Utf8, "protein_id": pl.Utf8},
)

PAIRWISE_AA = TableSchema(
    name="pairwise_aa",
    required={"qseqid": pl.Utf8, "sseqid": pl.Utf8, "pident": pl.Float64},
    optional={
        "length": pl.Int64,
        "evalue": pl.Float64,
        "bitscore": pl.Float64,
    },
)

PAIRWISE_NT = TableSchema(
    name="pairwise_nt",
    required={
        "query": pl.Utf8,
        "ref": pl.Utf8,
        "ani": pl.Float64,
        "query_start": pl.Int64,
        "query_end": pl.Int64,
        "ref_start": pl.Int64,
        "ref_end": pl.Int64,
    },
)

TREE_META = TableSchema(
    name="tree_metadata",
    required={"leaf_id": pl.Utf8},
    optional={"taxon": pl.Utf8},
)


def ensure_schema(
    df: pl.DataFrame, schema: TableSchema, *, allow_extra: bool = True
) -> pl.DataFrame:
    return schema.ensure(df, allow_extra=allow_extra)


def select_existing(df: pl.DataFrame, cols: Sequence[str]) -> pl.DataFrame:
    """Return a DataFrame with only existing columns from cols."""
    keep = [c for c in cols if c in df.columns]
    return df.select(keep)


def ensure_columns(df: pl.DataFrame, required: Iterable[str]) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
