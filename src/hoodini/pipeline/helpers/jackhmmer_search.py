import argparse
import multiprocessing
from datetime import datetime
from functools import partial

import polars as pl
import pyhmmer
from networkx.utils.union_find import UnionFind
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)


def process_sequence(seq, sequences, max_iterations=3):
    alphabet = pyhmmer.easel.Alphabet.amino()
    pli = pyhmmer.plan7.Pipeline(alphabet)
    name = seq.name.decode("utf-8")
    iterator = pli.iterate_seq(seq, sequences)
    for iteration in range(max_iterations):
        iteration = next(iterator)
        if iteration.converged:
            break
    return name, iteration.hits


def run_jackhmmer(faa, cpus=multiprocessing.cpu_count(), max_iterations=3):
    with pyhmmer.easel.SequenceFile(faa, format="fasta", digital=True) as seq_file:
        sequences = seq_file.read_block()

    process_sequence_partial = partial(
        process_sequence, sequences=sequences, max_iterations=max_iterations
    )

    dicc_hits = {}
    with (
        multiprocessing.Pool(cpus) as pool,
        Progress(
            TextColumn(
                f"[grey53][[/grey53][light_slate_grey]{datetime.now():%H:%M:%S}[/light_slate_grey][grey53]][/grey53]"
            ),
            SpinnerColumn(),
            TextColumn("{task.description}"),
            BarColumn(bar_width=40),
            TimeElapsedColumn(),
            TimeRemainingColumn(),
        ) as progress,
    ):
        task = progress.add_task("Running jackhmmer clustering...", total=len(sequences))
        for name, hits in pool.imap_unordered(process_sequence_partial, sequences):
            dicc_hits[name] = hits
            progress.update(task, advance=1)

    return dicc_hits


def cluster_jackhmmer_results(dicc_hits, min_evalue=1e-10):
    ds = UnionFind()
    for _, hits in dicc_hits.items():
        ds.union(*[hit.name.decode() for hit in hits if hit.evalue <= min_evalue])
    data = []
    for group in ds.to_sets():
        representative = next(iter(group))
        for member in group:
            data.append([representative, member])
    df = pl.DataFrame(data, columns=["clu_rep_seq", "member"])
    return df


def parallel_jackhmmer(faa, cpus=multiprocessing.cpu_count(), max_iterations=3, min_evalue=1e-10):
    dicc_hits = run_jackhmmer(faa, cpus, max_iterations)
    df = cluster_jackhmmer_results(dicc_hits, min_evalue)
    return df


def main():
    parser = argparse.ArgumentParser(description="Run parallel JackHMMER and cluster results")
    parser.add_argument("-f", "--fasta", required=True, help="FASTA file path")
    parser.add_argument(
        "-c", "--cpus", type=int, default=multiprocessing.cpu_count(), help="Number of CPUs to use"
    )
    parser.add_argument(
        "-m", "--max_iterations", type=int, default=3, help="Maximum number of iterations"
    )
    parser.add_argument("-e", "--min_evalue", type=float, default=1e-10, help="Minimum e-value")
    parser.add_argument("-o", "--output", required=True, help="Output file path")
    args = parser.parse_args()

    df = parallel_jackhmmer(args.fasta, args.cpus, args.max_iterations, args.min_evalue)

    df.write_csv(args.output, include_header=False)


if __name__ == "__main__":
    main()
