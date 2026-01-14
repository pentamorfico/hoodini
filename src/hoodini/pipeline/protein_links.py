import subprocess
import tempfile
from pathlib import Path

import polars as pl

from hoodini.utils.logging_utils import info, warn


def run_protein_links(
    output_dir: str, all_prots: pl.DataFrame, threads: int = 4, evalue: float = 1e-5
) -> pl.DataFrame:
    """Build Diamond DB from `all_prots` (DataFrame) and run all-vs-all blastp.

    Parameters
    - output_dir: base output folder where results.fasta will be read/written
    - all_prots: polars DataFrame containing at least columns ['id', 'sequence'] or ['protein_id', 'sequence']
    - threads: number of threads for Diamond
    - evalue: evalue cutoff for blastp

    Returns
    - polars.DataFrame with columns [qseqid, sseqid, pident, length, evalue, bitscore]
      with self-hits removed (same genome prefix in qseqid and sseqid).
    """
    df = all_prots.clone()
    if "sequence" not in df.columns:
        raise ValueError("all_prots must contain a 'sequence' column")

    if "id" in df.columns:
        id_col = "id"
    elif "protein_id" in df.columns:
        id_col = "protein_id"
    else:
        raise ValueError("all_prots must contain an 'id' or 'protein_id' column")

    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    fasta_path = out_dir / "results.fasta"
    db_path = out_dir / "results.dmnd"

    if not fasta_path.exists():
        info(f"Writing protein FASTA to {fasta_path}")
        with open(fasta_path, "w") as fh:
            for row in df.iter_rows(named=True):
                fh.write(f">{row[id_col]}\n")
                seq = row["sequence"]
                for i in range(0, len(seq), 80):
                    fh.write(seq[i : i + 80] + "\n")

    info("Building Diamond database...")
    makedb_cmd = [
        "diamond",
        "makedb",
        "--in",
        str(fasta_path),
        "--db",
        str(db_path),
        "--threads",
        str(threads),
    ]
    subprocess.run(makedb_cmd, check=True, capture_output=True)
    info(f"Diamond database built: {db_path}")

    info("Running all-vs-all blastp (Diamond)...")
    with tempfile.NamedTemporaryFile(mode="w+", suffix=".tsv", delete=False) as tmp_out:
        blastp_cmd = [
            "diamond",
            "blastp",
            "--db",
            str(db_path),
            "--query",
            str(fasta_path),
            "--out",
            tmp_out.name,
            "--outfmt",
            "6",
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "evalue",
            "bitscore",
            "--evalue",
            str(evalue),
            "--threads",
            str(threads),
        ]
        subprocess.run(blastp_cmd, check=True, capture_output=True)

        try:
            results_df = pl.read_csv(
                tmp_out.name,
                separator="\t",
                has_header=False,
                new_columns=["qseqid", "sseqid", "pident", "length", "evalue", "bitscore"],
            )
        except Exception:
            results_df = pl.DataFrame(
                schema={
                    "qseqid": pl.Utf8,
                    "sseqid": pl.Utf8,
                    "pident": pl.Float64,
                    "length": pl.Int64,
                    "evalue": pl.Float64,
                    "bitscore": pl.Float64,
                }
            )
        finally:
            Path(tmp_out.name).unlink(missing_ok=True)

    try:
        if db_path.exists():
            db_path.unlink()
    except Exception:
        warn(f"Could not remove diamond DB {db_path}")

    if results_df.height == 0:
        warn("No protein pairwise comparisons found")
        return results_df

    results_df = results_df.with_columns(
        [
            pl.col("qseqid").str.split("|").list.first().alias("qgenome"),
            pl.col("sseqid").str.split("|").list.first().alias("sgenome"),
        ]
    )

    filtered = results_df.filter(pl.col("qgenome") != pl.col("sgenome")).drop(
        ["qgenome", "sgenome"]
    )

    info(f"Protein pairwise comparisons complete: {filtered.height} hits (self-hits removed)")
    return filtered
