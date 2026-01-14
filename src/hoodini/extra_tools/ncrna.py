import re
import subprocess
from importlib.resources import files
from pathlib import Path

import polars as pl

from hoodini.utils.logging_utils import info, warn


def run_ncrna(all_neigh, den_data, output, num_threads, valid_unique_ids):
    info("🔬\tRunning Infernal for ncRNA annotation...")
    output = Path(output)
    ncrna_dir = output / "ncrna"
    ncrna_dir.mkdir(parents=True, exist_ok=True)
    cm_path = files("hoodini").joinpath("data", "all.cm")
    stockholm_file = ncrna_dir / "results.sto"
    tblout_file = ncrna_dir / "results.txt"
    command = [
        "cmsearch",
        "--tblout",
        str(tblout_file),
        "-A",
        str(stockholm_file),
        "-E",
        "1e-5",
        "--incE",
        "1e-5",
        "--cpu",
        str(num_threads),
        cm_path,
        str(output / "neighborhood" / "neighborhoods.fasta"),
    ]
    subprocess.run(command, check=True, stdout=subprocess.DEVNULL)
    column_names = [
        "nucid",
        "-",
        "nc_feature",
        "--",
        "cm",
        "mdlfrom",
        "mdlto",
        "seqfrom",
        "seqto",
        "strand_ncrna",
        "trunc",
        "pass",
        "gc",
        "bias",
        "score",
        "E-value",
        "inc",
        "desc",
    ]
    if stockholm_file.stat().st_size > 0:
        # Parse tblout file manually (whitespace-separated, comments start with #)
        rows = []
        with open(tblout_file) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = re.split(r"\s+", line.strip(), maxsplit=17)
                if len(parts) >= 17:
                    rows.append(parts[:18] if len(parts) >= 18 else parts + [""])

        if not rows:
            warn(f"No ncRNA found by Infernal (no valid rows in {tblout_file})")
            empty_df = pl.DataFrame()
            empty_df.write_csv(
                ncrna_dir / "ncrna_results.tsv", separator="\t", include_header=False
            )
            return empty_df

        cmdf = pl.DataFrame(rows, schema=column_names, orient="row")
        cmdf = cmdf.with_columns(
            [
                pl.col("seqfrom").cast(pl.Int64),
                pl.col("seqto").cast(pl.Int64),
            ]
        )

        # Build sequence and structure lookup from stockholm file
        seq_lookup = {}
        structure_lookup = {}

        from Bio import AlignIO

        for alignment in AlignIO.parse(stockholm_file, "stockholm"):
            # Get consensus secondary structure if available
            ss_cons = None
            if (
                hasattr(alignment, "column_annotations")
                and "secondary_structure" in alignment.column_annotations
            ):
                ss_cons = alignment.column_annotations["secondary_structure"]

            for record in alignment:
                # Parse sequence ID: seqid/start-end
                parts = record.id.split("/")
                seqid = parts[0]
                coords = parts[1].split("-")
                seqfrom = int(coords[0])
                seqto = int(coords[1])

                # Clean sequence (remove gaps)
                sequence = str(record.seq).replace(".", "").replace("-", "")
                seq_lookup[(seqid, seqfrom, seqto)] = sequence

                # Map structure to sequence (remove positions with gaps in sequence)
                # Convert to Vienna RNA format: . for unpaired, () for base pairs
                # Stockholm WUSS notation: https://en.wikipedia.org/wiki/Stockholm_format
                # Unpaired: . , ; : _ - ~
                # Base pairs (nested): <> () [] {}
                # Pseudoknots: Aa Bb Cc ... Zz (uppercase 5', lowercase 3')
                if ss_cons:
                    structure = ""
                    for i, char in enumerate(str(record.seq)):
                        if char not in ".-" and i < len(ss_cons):
                            ss_char = ss_cons[i]
                            # Convert Stockholm/WUSS to Vienna format
                            if ss_char in ".,;:_-~":
                                # Unpaired characters -> .
                                structure += "."
                            elif ss_char in "<([{" or ss_char.isupper():
                                # Opening base pairs (including pseudoknot 5' end) -> (
                                structure += "("
                            elif ss_char in ">)]}" or ss_char.islower():
                                # Closing base pairs (including pseudoknot 3' end) -> )
                                structure += ")"
                            else:
                                # Unknown character -> unpaired
                                structure += "."
                    structure_lookup[(seqid, seqfrom, seqto)] = structure

        # Add sequences and structures to dataframe
        sequences = []
        structures = []
        for row in cmdf.iter_rows(named=True):
            key = (row["nucid"], row["seqfrom"], row["seqto"])
            sequences.append(seq_lookup.get(key, ""))
            structures.append(structure_lookup.get(key, ""))
        cmdf = cmdf.with_columns(
            [
                pl.Series("sequence", sequences),
                pl.Series("structure", structures),
            ]
        )

        valid = all_neigh.filter(pl.col("unique_id").is_in([str(n) for n in valid_unique_ids]))[
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
        ]
        info(f"Parsed {cmdf.height} ncRNA hits from Infernal.")
        cmdf = cmdf.join(valid, left_on="nucid", right_on="temp_seqid", how="left")
        cmdf = cmdf.with_columns(
            (pl.col("seqfrom") + pl.col("start_win")).alias("start"),
            (pl.col("seqto") + pl.col("start_win")).alias("end"),
            pl.col("seqid").alias("nucid"),
            pl.col("unique_id").cast(pl.Utf8),
        )
        cmdf.write_csv(ncrna_dir / "ncrna_results.tsv", separator="\t", include_header=True)
        return cmdf

    else:
        warn(f"No ncRNA found by Infernal (empty {stockholm_file})")
        empty_df = pl.DataFrame()
        empty_df.write_csv(ncrna_dir / "ncrna_results.tsv", separator="\t", include_header=False)
        return empty_df
