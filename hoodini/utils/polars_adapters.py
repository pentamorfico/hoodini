"""Helpers to bridge pandas↔Polars at the edges of the codebase.

Goal: keep Polars as the internal dataframe representation. Use these adapters at
integration boundaries that still emit pandas (external libs, legacy code) until
migration is complete.
"""

from __future__ import annotations

from typing import Any, Iterable, Mapping

import pandas as pd
import polars as pl

from hoodini.models.schemas import TableSchema, ensure_schema


def to_polars(df: pl.DataFrame | pl.LazyFrame | pd.DataFrame | Any, *, schema: TableSchema | None = None) -> pl.DataFrame:
    """Convert incoming data to a Polars DataFrame and optionally enforce schema."""
    if isinstance(df, pl.DataFrame):
        out = df
    elif isinstance(df, pl.LazyFrame):
        out = df.collect()
    elif isinstance(df, pd.DataFrame):
        out = pl.from_pandas(df)
    else:
        raise TypeError(f"Unsupported dataframe type: {type(df)}; expected Polars")

    if schema:
        out = ensure_schema(out, schema)
    return out


def to_pandas(df: pl.DataFrame | pl.LazyFrame | Any) -> pd.DataFrame:
    """Convert Polars data to pandas for downstream libraries that require it."""
    if isinstance(df, pl.LazyFrame):
        df = df.collect()
    if isinstance(df, pl.DataFrame):
        return df.to_pandas()
    raise TypeError(f"Unsupported dataframe type: {type(df)}; expected Polars")


def ensure_required(df: pl.DataFrame, cols: Iterable[str]) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")


def rename_if_present(df: pl.DataFrame, mapping: Mapping[str, str]) -> pl.DataFrame:
    """Rename columns that exist; ignore missing keys."""
    present = {old: new for old, new in mapping.items() if old in df.columns}
    return df.rename(present) if present else df
