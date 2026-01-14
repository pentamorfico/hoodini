import subprocess
import tempfile
from pathlib import Path

import polars as pl

from hoodini.utils.logging_utils import info


def run_blast(all_neigh, output, blast, num_threads, valid_unique_ids):
    if blast:
        info("🔍\tRunning BLAST annotation...")
        neighborhood_fasta = f"{output}/neighborhood/neighborhoods.fasta"
        query = blast

        info("Creating BLAST database...")
        makeblastdb_cmd = [
            "makeblastdb",
            "-in",
            neighborhood_fasta,
            "-dbtype",
            "nucl",
            "-parse_seqids",
        ]
        subprocess.run(makeblastdb_cmd, check=True, capture_output=True)

        info("Running BLAST search...")
        with tempfile.NamedTemporaryFile(mode="w+", suffix=".tsv", delete=False) as tmp_out:
            blast_cmd = [
                "blastn",
                "-query",
                query,
                "-db",
                neighborhood_fasta,
                "-out",
                tmp_out.name,
                "-outfmt",
                "6 qseqid sseqid sstart send evalue",
                "-word_size",
                "8",
                "-evalue",
                "1e-5",
                "-dust",
                "no",
                "-reward",
                "1",
                "-penalty",
                "-2",
                "-gapopen",
                "6",
                "-gapextend",
                "2",
                "-num_threads",
                str(num_threads),
            ]
            subprocess.run(blast_cmd, check=True, capture_output=True)

            try:
                results_blast = pl.read_csv(
                    tmp_out.name,
                    separator="\t",
                    has_header=False,
                    new_columns=["qseqid", "sseqid", "sstart", "send", "evalue"],
                )
            except Exception:
                return pl.DataFrame()
            finally:
                Path(tmp_out.name).unlink(missing_ok=True)

        if results_blast.height == 0:
            return pl.DataFrame()

        if not isinstance(all_neigh, pl.DataFrame):
            all_neigh = pl.from_pandas(all_neigh)

        valid = all_neigh.filter(
            pl.col("unique_id").cast(pl.Utf8).is_in([str(n) for n in valid_unique_ids])
        ).select(
            [
                "seqid",
                "start_target",
                "end_target",
                "start_win",
                "end_win",
                "strand_win",
                "unique_id",
                "length",
                "temp_seqid",
            ]
        )

        results_blast = results_blast.join(
            valid, left_on="sseqid", right_on="temp_seqid", how="left"
        )

        results_blast = results_blast.with_columns(
            (pl.col("sstart") + pl.col("start_win")).alias("start"),
            (pl.col("send") + pl.col("start_win")).alias("end"),
        )

        results_blast = results_blast.rename({"qseqid": "nc_feature"})
        results_blast = results_blast.with_columns(
            ("BLAST " + pl.col("nc_feature")).alias("nc_feature"),
            pl.col("unique_id").cast(pl.Utf8).alias("unique_id"),
        )

        return results_blast

    return pl.DataFrame()
