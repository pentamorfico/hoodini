"""tRNA and tmRNA gene detection using pyaragorn.

Scans gene neighborhood nucleotide sequences for tRNA and tmRNA genes,
analogous to how ncrna.py uses pyinfernal for ncRNA detection. The output
DataFrame follows the same schema as ncRNA results so that both can be
merged seamlessly into the GFF and ncrna_metadata outputs.
"""

from __future__ import annotations

from pathlib import Path

import polars as pl
import pyaragorn
import pyhmmer

from hoodini.utils.logging_utils import info, warn


def run_trna(
    all_neigh: pl.DataFrame,
    den_data: pl.DataFrame,
    output: str | Path,
    num_threads: int,
    valid_unique_ids: list,
    *,
    translation_table: int = 11,
) -> pl.DataFrame:
    """Run pyaragorn tRNA/tmRNA detection on neighborhood sequences.

    Parameters
    ----------
    all_neigh : pl.DataFrame
        Neighborhood metadata.
    den_data : pl.DataFrame
        Taxonomy / tree metadata (unused, kept for API parity).
    output : str | Path
        Pipeline output directory.
    num_threads : int
        Number of threads (unused; pyaragorn's finder runs in-process per
        sequence, kept for API parity).
    valid_unique_ids : list
        Unique IDs of non-failed neighbourhoods to include.
    translation_table : int, optional
        Genetic code to use (default 11 = bacterial).

    Returns
    -------
    pl.DataFrame
        DataFrame with columns compatible with ncRNA results.
    """
    output = Path(output)
    trna_dir = output / "trna"
    trna_dir.mkdir(parents=True, exist_ok=True)

    fasta_path = output / "neighborhood" / "neighborhoods.fasta"
    if not fasta_path.exists():
        warn(f"Neighborhood FASTA not found: {fasta_path}.  Skipping tRNA/tmRNA search.")
        return pl.DataFrame()

    info("🔬\tRunning pyaragorn for tRNA/tmRNA detection...")

    # Neighborhoods are linear excerpts of a genome, not closed replicons.
    finder = pyaragorn.RNAFinder(translation_table, trna=True, tmrna=True, linear=True)

    rows = []
    with pyhmmer.easel.SequenceFile(fasta_path, format="fasta") as seq_file:
        for record in seq_file:
            seqid = record.name
            sequence = record.sequence
            for gene in finder.find_rna(sequence):
                if isinstance(gene, pyaragorn.TRNAGene):
                    nc_feature = f"tRNA-{gene.amino_acid}({gene.anticodon})"
                    tag_peptide = ""
                else:
                    nc_feature = "tmRNA"
                    tag_peptide = gene.peptide().rstrip("*")

                rows.append(
                    {
                        "nucid": seqid,
                        "nc_feature": nc_feature,
                        "seqfrom": gene.begin,
                        "seqto": gene.end,
                        "strand_ncrna": "-" if gene.strand < 0 else "+",
                        "score": gene.energy,
                        "E-value": ".",
                        "sequence": gene.sequence(),
                        "structure": "",
                        "tag_peptide": tag_peptide,
                    }
                )

    if not rows:
        info("   No tRNA/tmRNA genes found.")
        empty = pl.DataFrame()
        empty.write_csv(trna_dir / "trna_results.tsv", separator="\t", include_header=False)
        return empty

    info(f"   Found {len(rows)} tRNA/tmRNA genes across all neighborhoods")

    trna_df = pl.DataFrame(rows)

    # Map neighbourhood-local coords -> absolute genomic coords
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

    trna_df = trna_df.join(valid, left_on="nucid", right_on="temp_seqid", how="left")
    trna_df = trna_df.with_columns(
        (pl.col("seqfrom") + pl.col("start_win")).alias("start"),
        (pl.col("seqto") + pl.col("start_win")).alias("end"),
        pl.col("seqid").alias("nucid"),
        pl.col("unique_id").cast(pl.Utf8),
    )

    info(f"   Mapped {trna_df.height} tRNA/tmRNA hits to genomic coordinates")
    trna_df.write_csv(trna_dir / "trna_results.tsv", separator="\t", include_header=True)
    return trna_df
