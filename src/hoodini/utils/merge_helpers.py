"""Lightweight dataframe merge helpers."""

from __future__ import annotations

import polars as pl


def merge_cluster_result(all_prots: pl.DataFrame, merged: pl.DataFrame) -> pl.DataFrame:
    if "fam_cluster" in merged.columns:
        merged = merged.rename({"fam_cluster": "fam_cluster"})
    return all_prots.join(merged, on="id", how="left")
