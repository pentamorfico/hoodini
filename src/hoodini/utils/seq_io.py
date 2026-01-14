"""FASTA helpers and simple file transforms."""

from __future__ import annotations

from pathlib import Path

import polars as pl
from Bio import SeqIO


def to_fasta(df: pl.DataFrame, id_col: str, seq_col: str, path: str | Path) -> None:
    """Write a Polars DataFrame to FASTA."""
    dest = Path(path)
    with dest.open("w") as f:
        for row in df.iter_rows(named=True):
            f.write(f">{row[id_col]}\n{row[seq_col]}\n")


def _df_to_fasta(self: pl.DataFrame, id_col: str, seq_col: str, path: str | Path) -> None:
    to_fasta(self, id_col, seq_col, path)


pl.DataFrame.to_fasta = _df_to_fasta  # type: ignore[attr-defined]


def read_fasta(filename: str | Path) -> pl.DataFrame:
    seqs = []
    for record in SeqIO.parse(str(filename), "fasta"):
        seqs.append({"id": record.id, "sequence": str(record.seq)})
    return pl.DataFrame(seqs)
