import subprocess
from pathlib import Path

import polars as pl

from hoodini.utils.logging_utils import info, warn


def run_defensefinder(all_gff, all_prots, output):
    info("🛡️\tRunning DefenseFinder...")
    deffinder_df = pl.DataFrame()
    try:
        if not isinstance(all_gff, pl.DataFrame):
            all_gff = pl.from_pandas(all_gff)
        if not isinstance(all_prots, pl.DataFrame):
            all_prots = pl.from_pandas(all_prots)

        temp_gff = all_gff.clone()

        temp_gff = temp_gff.with_columns(
            pl.col("attributes").str.extract(r"ID=([^;]+)").alias("id")
        )
        temp_gff = temp_gff.join(
            all_prots.select(["id", "unique_id", "sequence"]), on="id", how="left"
        )

        temp_gff = temp_gff.select(
            [
                "seqid",
                "source",
                "type",
                "start",
                "end",
                "score",
                "strand",
                "phase",
                "id",
                "sequence",
            ]
        ).with_row_count("temp_index")

        temp_gff = temp_gff.with_columns(
            (pl.col("seqid") + "_" + pl.col("temp_index").cast(pl.Utf8)).alias("fasta_id")
        )

        output = Path(output)
        with open(output / "proteome.fasta", "w") as f:
            for row in temp_gff.iter_rows(named=True):
                if row["sequence"] is not None:
                    f.write(f">{row['fasta_id']}\n{row['sequence']}\n")

        command = [
            "defense-finder",
            "run",
            str(output / "proteome.fasta"),
            "-a",
            "--db-type",
            "gembase",
            "-o",
            str(output / "defense_finder"),
        ]
        subprocess.run(command, check=True)

        result_file = output / "defense_finder" / "proteome_defense_finder_genes.tsv"
        if result_file.exists():
            deffinder_result = pl.read_csv(result_file, separator="\t")
            if deffinder_result.height > 0:
                deffinder_result = deffinder_result.with_columns(
                    pl.col("hit_id").str.split("_").list.last().cast(pl.Int32).alias("temp_index")
                )

                deffinder_result = deffinder_result.with_columns(
                    pl.col("type").alias("deffinder_type"),
                    pl.col("subtype").alias("deffinder_subtype"),
                )

                def extract_gene(gene_name, def_type):
                    if gene_name is None:
                        return ""
                    try:
                        if def_type == "Cas":
                            return gene_name.split("_")[0]
                        else:
                            parts = gene_name.split("__")
                            return parts[1] if len(parts) > 1 else gene_name
                    except (IndexError, AttributeError):
                        return ""

                deffinder_result = deffinder_result.with_columns(
                    pl.struct(["gene_name", "deffinder_type"])
                    .map_elements(lambda s: extract_gene(s["gene_name"], s["deffinder_type"]))
                    .alias("deffinder_gene")
                )

                temp_gff = temp_gff.join(deffinder_result, on="temp_index", how="left")

                deffinder_df = (
                    temp_gff.select(["id", "deffinder_type", "deffinder_subtype", "deffinder_gene"])
                    .filter(pl.col("deffinder_type").is_not_null())
                    .unique(subset=["id"])
                )
            else:
                warn("No defensefinder results found")
        else:
            warn("DefenseFinder output file not found")

    except Exception as e:
        warn(f"DefenseFinder failed: {e}")

    return deffinder_df
