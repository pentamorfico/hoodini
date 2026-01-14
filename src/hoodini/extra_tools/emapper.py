import subprocess
from importlib.resources import files
from pathlib import Path
from shutil import copyfile

import polars as pl

from hoodini.utils.logging_utils import info, success, warn


def run_emapper(all_prots: pl.DataFrame, output: str | Path, num_threads: int = 1) -> pl.DataFrame:
    """
    Run mmseqs easy-search, pick best hit per query directly in Polars,
    join to eggNOG metadata, pick the deepest OG per query,
    and return one row per input protein as a pandas DataFrame.
    """

    info("🧾\tRunning eggNOG-mapper (mmseqs + eggNOG, best+deepest OG via Polars) ...")

    output = Path(output)
    emapper_dir = output / "emapper"
    emapper_dir.mkdir(parents=True, exist_ok=True)

    fasta_path = output / "results.faa"
    fasta_fallback = output / "results.fasta"

    if not fasta_path.exists():
        if fasta_fallback.exists():
            copyfile(fasta_fallback, fasta_path)
            info(f"Copied {fasta_fallback} -> {fasta_path}")
        else:
            seq_df = all_prots[["id", "sequence"]].drop_nulls().drop_duplicates("id")
            seq_df.to_fasta("id", "sequence", fasta_path)
            success(f"Generated {fasta_path}")

    mmseqs_dir = files("hoodini").joinpath("data", "emapper", "mmseqs")
    mmseqs_db_padded = str(mmseqs_dir.joinpath("mmseqs.db_pad"))
    mmseqs_db_unpadded = str(mmseqs_dir.joinpath("mmseqs.db"))

    results_m8 = emapper_dir / "results.m8"
    tmpdir = emapper_dir / "mmseqs_tmp"
    tmpdir.mkdir(parents=True, exist_ok=True)

    def run_mmseqs(db_prefix: str, use_gpu: bool):
        cmd = [
            "mmseqs",
            "easy-search",
            str(fasta_path),
            db_prefix,
            str(results_m8),
            str(tmpdir),
            "--threads",
            str(max(1, int(num_threads or 1))),
        ]
        if use_gpu:
            cmd += ["--gpu", "1"]
        info(f"Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)

    try:
        run_mmseqs(mmseqs_db_padded, use_gpu=True)
    except subprocess.CalledProcessError:
        warn("GPU search failed — retrying CPU (padded DB)...")
        try:
            run_mmseqs(mmseqs_db_padded, use_gpu=False)
        except subprocess.CalledProcessError:
            warn("CPU (padded DB) failed — trying unpadded DB GPU...")
            try:
                run_mmseqs(mmseqs_db_unpadded, use_gpu=True)
            except subprocess.CalledProcessError:
                warn("GPU (unpadded DB) failed — CPU unpadded DB...")
                run_mmseqs(mmseqs_db_unpadded, use_gpu=False)

    if not results_m8.exists():
        warn(f"mmseqs results not found at {results_m8}")
        return pl.DataFrame()

    hits_all = pl.read_csv(
        results_m8,
        has_header=False,
        separator="\t",
        new_columns=[
            "qseqid",
            "sseqid",
            "pident",
            "alnlen",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
        ],
    )

    hits_best = (
        hits_all.sort(["qseqid", "bitscore"], descending=[False, True])
        .group_by("qseqid")
        .head(1)
        .select(["qseqid", "sseqid"])
    )

    eggnog_prots_path = str(files("hoodini").joinpath("data", "emapper", "eggnog_prots.parquet"))
    eggnog_og_path = str(files("hoodini").joinpath("data", "emapper", "eggnog_og.parquet"))

    prots = (
        pl.scan_parquet(eggnog_prots_path)
        .with_columns(pl.col("ogs").fill_null("").str.split(",").alias("ogs_list"))
        .explode("ogs_list")
        .filter((pl.col("ogs_list") != "") & pl.col("ogs_list").str.contains("@"))
        .with_columns(pl.col("ogs_list").str.split_exact("@", 1).alias("og_split"))
        .with_columns(
            [
                pl.col("og_split").struct.field("field_0").alias("og"),
                pl.col("og_split").struct.field("field_1").cast(pl.Utf8).alias("level"),
            ]
        )
        .drop(["ogs_list", "og_split"])
    )

    hits_prots = hits_best.lazy().join(prots, left_on="sseqid", right_on="name", how="left")

    og = pl.scan_parquet(eggnog_og_path).with_columns(pl.col("level").cast(pl.Utf8))

    annotated = hits_prots.join(og, on=["og", "level"], how="left", suffix="_og").collect()

    annotated = annotated.with_columns(pl.col("level").cast(pl.Int64, strict=False))
    annotated = annotated.sort(["qseqid", "level"], descending=[False, True])
    deepest = annotated.group_by("qseqid").head(1)
    deepest = deepest.with_columns(pl.col("level").cast(pl.Utf8))

    deepest = deepest.rename({"qseqid": "id"})
    exclude_cols = {"level", "nm", "og", "ogs", "orthoindex", "sseqid", "name"}
    lead = ["id", "pname", "description", "COG_categories", "pfam"]
    lead_present = [c for c in lead if c in deepest.columns]
    rest = [c for c in deepest.columns if c not in lead_present and c not in exclude_cols]
    deepest = deepest.select(lead_present + rest)

    info("🔎 Head of annotated Polars DF (deepest OG per best hit):")
    info(deepest.head(10))
    info(f"shape: {deepest.shape}")

    success(f"mmseqs annotations ready: {deepest.height} queries annotated")
    return deepest
