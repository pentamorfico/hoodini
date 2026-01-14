import subprocess
from pathlib import Path

import polars as pl

from hoodini.utils.logging_utils import info


def run_padloc(all_gff, all_prots, output: str | Path, num_threads):
    info("🧬\tRunning PADLOC...")

    if not isinstance(all_gff, pl.DataFrame):
        all_gff = pl.from_pandas(all_gff)
    if not isinstance(all_prots, pl.DataFrame):
        all_prots = pl.from_pandas(all_prots)

    temp_gff = all_gff.clone()
    if "id" not in temp_gff.columns:
        if "attributes" in temp_gff.columns:
            temp_gff = temp_gff.with_columns(
                pl.col("attributes").str.extract(r"ID=([^;]+)").cast(pl.Utf8).alias("id")
            )
        elif "protein_id" in temp_gff.columns:
            temp_gff = temp_gff.with_columns(pl.col("protein_id").cast(pl.Utf8).alias("id"))
        else:
            temp_gff = temp_gff.with_columns(pl.lit(None).cast(pl.Utf8).alias("id"))
    temp_gff = temp_gff.join(all_prots.select(["id", "unique_id", "sequence"]), on="id", how="left")
    output = Path(output)

    if temp_gff.filter(pl.col("sequence").is_not_null()).is_empty():
        fasta_path = output / "results.fasta"
        if fasta_path.exists():
            ids = []
            seqs = []
            with open(fasta_path) as fh:
                curr_id = None
                for line in fh:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith(">"):
                        curr_id = line[1:]
                        ids.append(curr_id)
                        seqs.append("")
                    elif seqs:
                        seqs[-1] += line
            fasta_df = pl.DataFrame({"id": ids, "sequence": seqs})
            temp_gff = temp_gff.join(fasta_df, on="id", how="left")
    temp_gff = temp_gff.with_columns(
        (
            pl.col("id").cast(pl.Utf8)
            + "-"
            + pl.when(pl.col("unique_id").is_not_null())
            .then(pl.col("unique_id").cast(pl.Utf8))
            .otherwise(pl.col("id").cast(pl.Utf8))
            + "-"
            + pl.col("start").cast(pl.Utf8)
        ).alias("temp_id")
    )
    temp_gff = temp_gff.with_columns(("ID=" + pl.col("temp_id")).alias("attributes"))
    temp_gff = temp_gff.filter(pl.col("temp_id").is_not_null() & pl.col("sequence").is_not_null())
    temp_gff = temp_gff.unique(subset=["attributes", "seqid"])
    temp_gff = temp_gff.unique(subset=["seqid", "start", "end"])

    temp_gff_out = temp_gff.select(
        ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
    )
    temp_gff_out.write_csv(output / "temp.gff", separator="\t", include_header=False)

    temp_fasta = (
        temp_gff.select(["temp_id", "sequence"])
        .filter(pl.col("temp_id").is_not_null() & pl.col("sequence").is_not_null())
        .unique(subset=["temp_id"])
    )
    temp_fasta_path = output / "temp.fasta"
    with open(temp_fasta_path, "w") as f:
        for row in temp_fasta.iter_rows(named=True):
            f.write(f">{row['temp_id']}\n{row['sequence']}\n")

    padloc_dir = output / "padloc"
    padloc_dir.mkdir(parents=True, exist_ok=True)
    command = [
        "padloc",
        "--faa",
        str(temp_fasta_path),
        "--gff",
        str(output / "temp.gff"),
        "--cpu",
        str(num_threads),
        "--outdir",
        str(padloc_dir),
    ]
    subprocess.run(command, check=True)
    result_path = padloc_dir / "temp.fasta_padloc.csv"
    if result_path.exists():
        padloc_df = pl.read_csv(result_path)
        if padloc_df.height > 0:
            padloc_df = padloc_df.rename({"system": "padloc_system", "protein.name": "padloc_gene"})
            padloc_df = padloc_df.with_columns(
                pl.col("target.name").str.split("-").list.first().alias("target.name")
            )
            padloc_df = padloc_df.with_columns(pl.col("target.name").alias("id"))
            padloc_df = padloc_df.select(["id", "padloc_system", "padloc_gene"])
            padloc_df = padloc_df.unique(subset=["id"])
    else:
        padloc_df = pl.DataFrame()

    return padloc_df
