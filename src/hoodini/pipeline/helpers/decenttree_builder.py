"""Run DecentTree on a pairwise distance table and return Newick string.

The table should contain rows with (query_id, target_id, distance).
This module will build a relaxed PHYLIP-like distance matrix (tab-separated,
no name length limits), run DecentTree (bundled in `hoodini/extra_tools`) and
return the resulting Newick string.

The function `run_decenttree_from_table` accepts column-name arguments so it
can be used with different input DataFrame schemas.
"""

from __future__ import annotations

import contextlib
import csv
import platform
import subprocess
import tempfile
from collections.abc import Iterable
from importlib.resources import files
from pathlib import Path

import polars as pl


def _choose_decenttree_binary() -> Path:
    sysname = platform.system().lower()
    if sysname.startswith("darwin"):
        p = files("hoodini").joinpath("extra_tools", "decenttree_macos")
    else:
        p = files("hoodini").joinpath("extra_tools", "decenttree_linux")
    if not p.exists():
        raise FileNotFoundError(f"DecentTree binary not found at {p}")
    return p


def _normalize_algorithm(name: str) -> str:
    """Map common algorithm aliases to DecentTree accepted names.

    Accepts case-insensitive inputs like 'nj', 'nj-r', 'bionj' and returns
    the proper DecentTree token (e.g. 'NJ', 'NJ-R', 'BIONJ').
    """
    if not name:
        return "NJ"
    n = str(name).strip().lower()
    mapping = {
        "nj": "NJ",
        "nj-r": "NJ-R",
        "nj_r": "NJ-R",
        "njr": "NJ-R",
        "nj-v": "NJ-V",
        "bionj": "BIONJ",
        "bionj-r": "BIONJ-R",
        "rapidnj": "RapidNJ",
        "unj": "UNJ",
        "upgma": "UPGMA",
        "auction": "AUCTION",
        "bionj-v": "BIONJ-V",
        "nj-r-d": "NJ-R-D",
    }
    if n in mapping:
        return mapping[n]
    cand = str(name).upper()
    return cand


def _build_relaxed_phylip(
    df: pl.DataFrame, qcol: str, tcol: str, dcol: str, out_path: Path, fill_missing=None
) -> None:
    """Create a relaxed PHYLIP-like matrix (tab-separated) at out_path.

    The matrix will include all unique sequence ids found in qcol and tcol,
    and distances will be placed in the appropriate cells. Missing entries will
    be set to 0 (or a large value if you prefer).
    """
    if hasattr(fill_missing, "__iter__") and not isinstance(fill_missing, str):
        fill_value, ids = fill_missing
    else:
        ids = (
            pl.concat(
                [
                    df.select(pl.col(qcol).cast(pl.Utf8)).to_series(),
                    df.select(pl.col(tcol).cast(pl.Utf8)).to_series(),
                ]
            )
            .unique()
            .to_list()
        )
        fill_value = fill_missing
    ids = list(ids)
    n = len(ids)

    if fill_value is None:
        fill_value = 0.0
    elif isinstance(fill_value, str):
        sval = fill_value.lower()
        vals_series = df.select(pl.col(dcol).cast(pl.Float64, strict=False)).to_series()
        vals = vals_series.drop_nulls().to_numpy()
        if vals.size == 0:
            fill_value = 0.0
        elif sval == "max":
            fill_value = float(vals.max())
        elif sval in ("max+2std", "max_plus_2std", "max+2*std"):
            fill_value = float(vals.max()) + 2.0 * float(vals.std())
        elif sval == "min":
            fill_value = float(vals.min())
        else:
            fill_value = 0.0
    else:
        try:
            fill_value = float(fill_value)
        except Exception:
            fill_value = 0.0

    mat = [[fill_value] * n for _ in range(n)]

    df_lookup = dict(
        zip(
            zip(df[qcol].to_list(), df[tcol].to_list()),
            df[dcol].to_list(),
        )
    )

    for i, q in enumerate(ids):
        for j, t in enumerate(ids):
            v = None
            if (q, t) in df_lookup:
                v = df_lookup.get((q, t))
            elif (t, q) in df_lookup:
                v = df_lookup.get((t, q))
            if v is not None and v is not None:
                try:
                    d = float(v)
                    mat[i][j] = d
                except Exception:
                    pass

    out_path_parent = out_path.parent
    out_path_parent.mkdir(parents=True, exist_ok=True)

    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_MINIMAL)
        fh.write(f"{n}\n")
        for i, idv in enumerate(ids):
            row = [idv] + [str(x) for x in mat[i]]
            writer.writerow(row)


def run_decenttree_from_table(
    df: pl.DataFrame,
    qcol: str = "query",
    tcol: str = "target",
    dcol: str = "distance",
    algorithm: str = "nj",
    threads: int = 1,
    decenttree_bin: Path | None = None,
    fill_missing=None,
    ids: Iterable[str] | None = None,
) -> str:
    """Run DecentTree on a Polars pairwise distance table and return a Newick string."""
    if decenttree_bin is None:
        decenttree_bin = _choose_decenttree_binary()

    dt_df = df.clone() if hasattr(df, "clone") else df.select(list(df.columns))

    if "AAI" in dt_df.columns and dcol not in dt_df.columns:
        with contextlib.suppress(Exception):
            dt_df = dt_df.with_columns(pl.col("AAI").cast(pl.Float64))
        aai_max = dt_df["AAI"].max()
        if aai_max is not None and aai_max <= 1.0:
            dt_df = dt_df.with_columns((pl.col("AAI") * 100.0).alias("AAI"))
        dt_df = dt_df.with_columns((100.0 - pl.col("AAI")).alias("distance"))
        dcol_use = "distance"
    elif "ANI" in dt_df.columns and dcol not in dt_df.columns:
        with contextlib.suppress(Exception):
            dt_df = dt_df.with_columns(pl.col("ANI").cast(pl.Float64))
        ani_max = dt_df["ANI"].max()
        if ani_max is not None and ani_max <= 1.0:
            dt_df = dt_df.with_columns((pl.col("ANI") * 100.0).alias("ANI"))
        dt_df = dt_df.with_columns((100.0 - pl.col("ANI")).alias("distance"))
        dcol_use = "distance"
    else:
        dcol_use = dcol

    if qcol not in dt_df.columns or tcol not in dt_df.columns:
        raise ValueError(f"Input DataFrame must contain columns: {qcol}, {tcol}")

    if ids is not None:
        valid_uids = [str(x) for x in ids]
    else:
        valid_uids = (
            pl.concat(
                [
                    dt_df.select(pl.col(qcol).cast(pl.Utf8)).to_series(),
                    dt_df.select(pl.col(tcol).cast(pl.Utf8)).to_series(),
                ]
            )
            .unique()
            .to_list()
        )

    with tempfile.TemporaryDirectory() as td:
        tdpth = Path(td)
        phylip_path = tdpth / "matrix.phylip"
        out_newick = tdpth / "tree.nwk"

    _build_relaxed_phylip(
        dt_df, qcol, tcol, dcol_use, phylip_path, fill_missing=(fill_missing, valid_uids)
    )
    alg = _normalize_algorithm(algorithm)
    cmd = [
        str(decenttree_bin),
        "-in",
        str(phylip_path),
        "-out",
        str(out_newick),
        "-t",
        alg,
        "-nt",
        str(int(threads)),
    ]
    proc = subprocess.run(cmd, check=False, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"DecentTree failed: rc={proc.returncode}\nstderr={proc.stderr}\nstdout={proc.stdout}"
        )

    if not out_newick.exists():
        raise RuntimeError("DecentTree did not produce output tree")

    newick = out_newick.read_text(encoding="utf-8")
    return newick.strip()


def run_decenttree_from_matrix(
    matrix_df: pl.DataFrame,
    algorithm: str = "nj",
    threads: int = 1,
    decenttree_bin: Path | None = None,
) -> str:
    """Run DecentTree on a square distance matrix (Polars) and return Newick."""
    if not isinstance(matrix_df, pl.DataFrame):
        raise TypeError("matrix_df must be a Polars DataFrame with an 'id' column")

    if "id" not in matrix_df.columns:
        raise ValueError("matrix_df must contain an 'id' column for row labels")

    cols = [c for c in matrix_df.columns if c != "id"]
    if len(cols) != matrix_df.height:
        raise ValueError(
            "matrix_df must be square with one 'id' column and the remaining columns as samples"
        )

    if decenttree_bin is None:
        decenttree_bin = _choose_decenttree_binary()

    with tempfile.TemporaryDirectory() as td:
        tdpth = Path(td)
        phylip_path = tdpth / "matrix.phylip"
        out_newick = tdpth / "tree.nwk"

        ids = [str(x) for x in matrix_df["id"].to_list()]
        with open(phylip_path, "w", newline="") as fh:
            fh.write(f"{len(ids)}\n")
            writer = csv.writer(fh, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_MINIMAL)
            for i, idv in enumerate(ids):
                row_vals = [matrix_df[c][i] for c in cols]
                row = [idv] + [str(x) for x in row_vals]
                writer.writerow(row)
        alg = _normalize_algorithm(algorithm)
        cmd = [
            str(decenttree_bin),
            "-in",
            str(phylip_path),
            "-out",
            str(out_newick),
            "-t",
            alg,
            "-nt",
            str(int(threads)),
        ]
        proc = subprocess.run(cmd, check=False, capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(
                f"DecentTree failed: rc={proc.returncode}\nstderr={proc.stderr}\nstdout={proc.stdout}"
            )

        if not out_newick.exists():
            raise RuntimeError("DecentTree did not produce output tree")

        newick = out_newick.read_text(encoding="utf-8")
        return newick.strip()
