import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

import polars as pl

from hoodini.utils.logging_utils import error


def _run_command(command):
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        error(f"Error executing command: {command}")
        error(result.stderr)
        sys.exit(1)
    return result.stdout


def read_output_to_df(output_file):
    df = pl.read_csv(output_file, separator="\t", header=None)
    return df


def cluster_with_mmseqs(
    fasta,
    temp_folder,
    max_steps=5,
    sensitivity=15,
    cluster_mode=1,
    cluster_steps=9,
    cov_mode=0,
    coverage=0.7,
    output=None,
):
    fasta = Path(fasta)
    temp_folder = Path(temp_folder)
    output = Path(output) if output is not None else None

    dbname = "temp_db"
    _run_command(f"mmseqs createdb {fasta}  {temp_folder / dbname}")
    _run_command(
        f"mmseqs cluster {temp_folder / dbname} {temp_folder / (dbname + '_clu')} {temp_folder} -s {sensitivity} --cluster-mode {cluster_mode} --cluster-steps {cluster_steps} --cov-mode {cov_mode} -c {coverage}"
    )
    _run_command(
        f"mmseqs createsubdb {temp_folder / (dbname + '_clu')} {temp_folder / dbname} {temp_folder / (dbname + '_clu_repseq')}"
    )

    input_db = f"{dbname}_clu_repseq"
    for step in range(1, int(max_steps) + 1):
        _run_command(
            f"mmseqs search {temp_folder / input_db} {temp_folder / input_db} {temp_folder / f'search_step_{step}'} {temp_folder} --add-self-matches"
        )
        _run_command(
            f"mmseqs result2profile {temp_folder / (dbname + '_clu_repseq')} {temp_folder / (dbname + '_clu_repseq')} {temp_folder / f'search_step_{step}'} {temp_folder / f'search_step_{step}_profile'}"
        )
        _run_command(
            f"mmseqs search {temp_folder / f'search_step_{step}_profile'} {temp_folder / input_db} {temp_folder / f'search_step_{step}_pp'} {temp_folder} --add-self-matches"
        )
        _run_command(
            f"mmseqs clust {temp_folder / f'search_step_{step}_profile'} {temp_folder / f'search_step_{step}_pp'} {temp_folder / f'search_step_{step}_pp_clu'}"
        )
        _run_command(
            f"mmseqs createsubdb {temp_folder / f'search_step_{step}_pp_clu'} {temp_folder / dbname} {temp_folder / f'search_step_{step}_pp_clu_repseq'}"
        )
        input_db = f"search_step_{step}_pp_clu_repseq"

    cluster_files = " ".join(
        [f"{temp_folder}/search_step_{step}_pp_clu" for step in range(1, int(max_steps) + 1)]
    )
    _run_command(
        f"mmseqs mergeclusters {temp_folder / dbname} {temp_folder}/deep_cluster_db {temp_folder / (dbname + '_clu')} {cluster_files}"
    )
    _run_command(
        f"mmseqs createtsv {temp_folder / dbname} {temp_folder / dbname} {temp_folder}/deep_cluster_db {output}"
    )


def main():
    parser = argparse.ArgumentParser(
        description="Perform deep clustering on protein sequences using MMseqs2"
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Input FASTA file containing protein sequences"
    )
    parser.add_argument(
        "-t", "--temp-folder", required=True, help="Folder to save intermediary files"
    )
    parser.add_argument(
        "-m",
        "--max-steps",
        type=int,
        default=5,
        help="Maximum number of iterative steps (default: 5)",
    )
    parser.add_argument(
        "-s", "--sensitivity", type=int, default=15, help="Sensitivity for clustering (default: 15)"
    )
    parser.add_argument("--cluster-mode", type=int, default=1, help="Clustering mode (default: 1)")
    parser.add_argument(
        "--cluster-steps", type=int, default=9, help="Number of clustering steps (default: 9)"
    )
    parser.add_argument(
        "--cov-mode", type=int, default=0, help="Coverage mode for clustering (default: 0)"
    )
    parser.add_argument(
        "-c", "--coverage", type=float, default=0.7, help="Coverage for clustering (default: 0.7)"
    )
    parser.add_argument("-o", "--output", required=True, help="Output file name")

    args = parser.parse_args()

    os.makedirs(args.temp_folder, exist_ok=True)
    cluster_with_mmseqs(
        args.input,
        args.temp_folder,
        args.max_steps,
        args.sensitivity,
        args.cluster_mode,
        args.cluster_steps,
        args.cov_mode,
        args.coverage,
        args.output,
    )
    shutil.rmtree(args.temp_folder)


if __name__ == "__main__":
    main()
