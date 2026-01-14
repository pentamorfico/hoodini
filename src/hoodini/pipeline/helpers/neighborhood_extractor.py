"""Modular and robust neighborhood extraction for GBFF and GFF+FAA inputs.

Architecture:
1. Parse input files → standardized feature DataFrame
2. Compute window coordinates based on mode
3. Extract features in window
4. Optionally annotate sORFs
5. Return structured results or raise clear exceptions
"""

from __future__ import annotations

import traceback
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import gb_io
import orfipy_core
import polars as pl
import pyrodigal
from Bio.Seq import Seq


@dataclass
class NeighborhoodResult:
    """Structured result from extract_neighborhood."""

    proteins: pl.DataFrame
    neighborhood: pl.DataFrame
    unique_id: str


class NeighborhoodExtractionError(Exception):
    """Base exception for neighborhood extraction failures."""


def parse_gbff(gbff_path: str, nucleotide_id: str | None = None) -> tuple[pl.DataFrame, object]:
    """Parse GBFF and return features DataFrame + record object.

    Returns:
        (features_df, record): Polars DataFrame with CDS features and the gb_io record

    Raises:
        NeighborhoodExtractionError: if file missing, parsing fails, or nucleotide_id not found
    """
    if not Path(gbff_path).exists():
        raise NeighborhoodExtractionError(f"GenBank file not found: {gbff_path}")

    try:
        records = gb_io.load(gbff_path)
        records_list = list(records)
    except Exception as e:
        raise NeighborhoodExtractionError(f"Failed to load GBFF {gbff_path}: {e}")

    if not records_list:
        raise NeighborhoodExtractionError(f"No records found in {gbff_path}")

    record = None
    if nucleotide_id:
        for rec in records_list:
            rec_version = getattr(rec, "version", "")

            if nucleotide_id in rec_version:
                record = rec
                break
        if record is None:
            raise NeighborhoodExtractionError(
                f"Nucleotide {nucleotide_id} not found in {gbff_path}. Available: {[getattr(r, 'version', '?') for r in records_list[:3]]}"
            )
    else:
        record = records_list[0] if records_list else None
        if record is None:
            raise NeighborhoodExtractionError(f"No records found in {gbff_path}")

    record_version = getattr(record, "version", "unknown")

    feature_data = _process_gbff_features(record.features, record_version)

    if not feature_data:
        raise NeighborhoodExtractionError(f"No CDS features found in {gbff_path}")

    features_df = pl.DataFrame(feature_data)

    if "attributes" in features_df.columns:
        attrs_list = features_df["attributes"].to_list()
        attrs_df = pl.DataFrame(attrs_list)

        if "protein_id" in attrs_df.columns:
            attrs_df = attrs_df.drop("protein_id")

        features_df = features_df.drop("attributes")

        features_df = features_df.with_row_count("_join_idx")
        attrs_df = attrs_df.with_row_count("_join_idx")

        features_df = features_df.join(attrs_df, on="_join_idx", how="inner")
        features_df = features_df.drop("_join_idx")

        if "translation" in features_df.columns:
            features_df = features_df.rename({"translation": "sequence"})

    if features_df.height == 0:
        raise NeighborhoodExtractionError("GenBank file is not annotated")

    return features_df, record


def _process_gbff_features(features, record_accession: str) -> list[dict]:
    """Extract CDS features from gb_io features list."""
    data = []
    for feature in features:
        if feature.kind != "CDS":
            continue

        location = feature.location
        qualifiers = {q.key: q.value for q in feature.qualifiers}
        protein_id = qualifiers.get("protein_id", "")

        if not protein_id:
            continue

        class_type = type(location).__name__

        if class_type == "Range":
            start, end, strand = location.start, location.end, "+"
        elif class_type == "Complement":
            start, end, strand = location.start, location.end, "-"
        elif class_type == "Join":
            start = location.locations[0].start
            end = location.locations[0].end
            if location.locations[-1].start in range(end - 5, end + 6):
                start = location.locations[0].start
                end = location.locations[-1].end
            strand = "-" if type(location.locations[0]).__name__ == "Complement" else "+"
        else:
            continue

        if end < start:
            start, end, strand = end, start, "-"

        data.append(
            {
                "seqid": record_accession,
                "source": "hoodini",
                "type": "CDS",
                "start": start,
                "end": end,
                "score": ".",
                "strand": strand,
                "phase": ".",
                "attributes": qualifiers,
                "protein_id": protein_id,
            }
        )

    return data


def parse_gff_faa(gff_path: str, faa_path: str) -> pl.DataFrame:
    """Parse GFF+FAA and return features DataFrame.

    Raises:
        NeighborhoodExtractionError: if files missing or parsing fails
    """
    if not Path(gff_path).exists():
        raise NeighborhoodExtractionError(f"GFF file not found: {gff_path}")
    if not Path(faa_path).exists():
        raise NeighborhoodExtractionError(f"FAA file not found: {faa_path}")

    gff_header = [
        "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes",
    ]

    try:
        gff = pl.read_csv(
            gff_path,
            separator="\t",
            has_header=False,
            new_columns=gff_header,
            skip_rows_after_header=0,
        )
        if gff.height > 0 and "seqid" in gff.columns:
            gff = gff.filter(~pl.col("seqid").str.starts_with("#"))
    except Exception as e:
        raise NeighborhoodExtractionError(f"Failed to read GFF file {gff_path}: {e}")

    try:
        faa_df = _read_fasta(faa_path)
    except Exception as e:
        raise NeighborhoodExtractionError(f"Failed to read FAA file {faa_path}: {e}")

    gff = gff.with_columns(pl.col("attributes").str.extract(r"ID=([^;]+)", 1).alias("protein_id"))

    features_df = gff.join(faa_df, left_on="protein_id", right_on="id", how="left")

    if features_df.height == 0:
        raise NeighborhoodExtractionError("No features found in GFF file")

    return features_df


def _read_fasta(filename: str) -> pl.DataFrame:
    """Read FASTA file into DataFrame with id and sequence columns."""
    with open(filename) as file:
        records = file.read().split(">")[1:]
        records = [record.split("\n", 1) for record in records]
        records = [(t[0].split(" ")[0], "".join(t[1].split())) for t in records]
    return pl.DataFrame(records, schema=["id", "sequence"])


def compute_window(
    features_df: pl.DataFrame,
    *,
    input_type: Literal["protein", "nucleotide"],
    mode: Literal["win_nts", "win_ngen"],
    window: int,
    protein_id: str | None = None,
    start: int | None = None,
    end: int | None = None,
    strand: str | None = None,
    nucleotide_id: str | None = None,
    sequence_length: int | None = None,
) -> tuple[int, int, str, int, int]:
    """Compute window coordinates (start_win, end_win, strand, start_target, end_target).

    Raises:
        NeighborhoodExtractionError: if protein_id or coordinates not found
    """
    if input_type == "protein":
        if not protein_id:
            raise NeighborhoodExtractionError("protein_id required for input_type=protein")

        match = features_df.filter(pl.col("protein_id") == protein_id)
        if match.height == 0:
            raise NeighborhoodExtractionError(f"Protein {protein_id} not found in features")

        prot_row = match.row(0, named=True)
        start = prot_row["start"]
        end = prot_row["end"]
        strand = prot_row["strand"]

        if mode == "win_nts":
            start_win = max(0, start - window)
            end_win = end + window
            if sequence_length:
                end_win = min(end_win, sequence_length)

        elif mode == "win_ngen":
            indexed = features_df.with_row_count("_idx")
            prot_idx = indexed.filter(pl.col("protein_id") == protein_id).row(0, named=True)["_idx"]
            slice_start = max(0, prot_idx - window)
            slice_end = min(indexed.height, prot_idx + window + 1)
            window_features = indexed.slice(slice_start, slice_end - slice_start)
            start_win = window_features["start"].min()
            end_win = window_features["end"].max()
        else:
            raise NeighborhoodExtractionError(f"Unknown mode: {mode}")

    elif input_type == "nucleotide":
        if start is None or end is None:
            start_win = 0
            end_win = sequence_length or features_df["end"].max()
            start = start_win
            end = end_win
            strand = strand or "+"
        else:
            start, end = int(start), int(end)
            if not strand:
                strand = "+" if end > start else "-"

            if mode == "win_nts":
                start_win = max(0, start - window)
                end_win = end + window
                if sequence_length:
                    end_win = min(end_win, sequence_length)

            elif mode == "win_ngen":
                indexed = features_df.with_row_count("_idx")
                in_range = indexed.filter((pl.col("start") >= start) & (pl.col("end") <= end))
                if in_range.height == 0:
                    raise NeighborhoodExtractionError("No features in specified nucleotide range")
                start_idx = in_range.row(0, named=True)["_idx"]
                end_idx = in_range.row(-1, named=True)["_idx"]
                slice_start = max(0, start_idx - window)
                slice_end = min(indexed.height, end_idx + window + 1)
                window_features = indexed.slice(slice_start, slice_end - slice_start)
                start_win = window_features["start"].min()
                end_win = window_features["end"].max()
            else:
                raise NeighborhoodExtractionError(f"Unknown mode: {mode}")
    else:
        raise NeighborhoodExtractionError(f"Unknown input_type: {input_type}")

    return start_win, end_win, strand, start, end


def extract_features_in_window(
    features_df: pl.DataFrame,
    start_win: int,
    end_win: int,
) -> pl.DataFrame:
    """Filter features to those within the window."""
    subset = features_df.filter((pl.col("start") >= start_win) & (pl.col("end") <= end_win))

    if "protein_id" in subset.columns:
        subset = subset.with_columns(pl.col("protein_id").alias("id"))

    if "product" not in subset.columns:
        subset = subset.with_columns(pl.lit(None).alias("product"))

    columns = [
        "seqid",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "protein_id",
        "id",
        "sequence",
    ]
    if "product" in subset.columns:
        columns.append("product")

    subset = subset.select([c for c in columns if c in subset.columns])

    return subset


def annotate_sorfs(
    features_df: pl.DataFrame,
    sequence: str,
    start_win: int,
    nucleotide_id: str,
    unique_id: str,
) -> pl.DataFrame:
    """Annotate sORFs using pyrodigal and orfipy, return updated features."""
    orf_finder = pyrodigal.GeneFinder(meta=True, min_gene=10, max_overlap=9)
    new_genes = []

    for i, pred in enumerate(orf_finder.find_genes(sequence)):
        overlap_flag = False
        for row in features_df.iter_rows(named=True):
            overlap_pct = _calculate_overlap(
                row["start"], row["end"], pred.begin + start_win, pred.end + start_win
            )
            if overlap_pct > 10:
                overlap_flag = True
                break

        if not overlap_flag:
            new_genes.append(
                {
                    "seqid": nucleotide_id,
                    "source": "pyrodigal",
                    "type": "CDS",
                    "start": pred.begin + start_win,
                    "end": pred.end + start_win,
                    "score": pred.score,
                    "strand": "-" if pred.strand == "-1" else "+",
                    "phase": ".",
                    "protein_id": f"sORF_{unique_id}_{i}",
                    "id": f"sORF_{unique_id}_{i}",
                    "sequence": pred.translate(),
                }
            )

    seq_upper = sequence.upper()
    for i, (start, stop, strand, _description) in enumerate(
        orfipy_core.orfs(seq_upper, minlen=100, maxlen=1000, partial3=False, between_stops=False)
    ):
        overlap_flag = False
        for row in features_df.iter_rows(named=True):
            overlap_pct = _calculate_overlap(
                row["start"], row["end"], start + start_win, stop + start_win
            )
            if overlap_pct > 0:
                overlap_flag = True
                break

        if not overlap_flag:
            orf_sequence = Seq(sequence[start:stop])
            if strand == "-":
                orf_sequence = orf_sequence.reverse_complement()
            protein_sequence = orf_sequence.translate(table=11, to_stop=True)

            new_genes.append(
                {
                    "seqid": nucleotide_id,
                    "source": "orfipy",
                    "type": "CDS",
                    "start": start + start_win,
                    "end": stop + start_win,
                    "score": ".",
                    "strand": "-" if strand == "-" else "+",
                    "phase": ".",
                    "protein_id": f"sORF_orfipy_{unique_id}_{i}",
                    "id": f"sORF_orfipy_{unique_id}_{i}",
                    "sequence": str(protein_sequence),
                }
            )

    if new_genes:
        new_genes_df = pl.DataFrame(new_genes)
        features_df = pl.concat([features_df, new_genes_df], how="vertical")

    return features_df


def _calculate_overlap(coord1A: int, coord1B: int, coord2A: int, coord2B: int) -> float:
    """Calculate overlap percentage of coord2 with coord1."""
    start1, end1 = sorted([coord1A, coord1B])
    start2, end2 = sorted([coord2A, coord2B])

    max_start = max(start1, start2)
    min_end = min(end1, end2)

    overlap = max(0, min_end - max_start)
    length_second = end2 - start2

    if length_second == 0:
        return 0.0

    return (overlap / length_second) * 100


def extract_neighborhood(
    protein_id: str | None,
    nucleotide_id: str | None,
    gbf_file: str | None,
    gff_file: str | None,
    faa_file: str | None,
    fna_file: str | None,
    mode: str = "win_nts",
    window: int = 20000,
    strand: str | None = None,
    start: int | None = None,
    end: int | None = None,
    unique_id: str | None = None,
    input_type: str | None = None,
    sorfs: bool = False,
) -> tuple[pl.DataFrame | None, pl.DataFrame | None, str, str | None]:
    """Extract neighborhood from GBFF or GFF+FAA. Returns (proteins_df, neighborhood_df, unique_id, error_msg).

    This is the main entry point called by multiprocessing pool.
    """
    try:
        if gbf_file:
            features_df, record = parse_gbff(gbf_file, nucleotide_id)
            record_version = getattr(record, "version", "unknown")
            sequence_bytes = record.sequence
            sequence_length = len(sequence_bytes)
        elif gff_file and faa_file:
            features_df = parse_gff_faa(gff_file, faa_file)
            record_version = nucleotide_id or "unknown"
            if fna_file and Path(fna_file).exists():
                fna_df = _read_fasta(fna_file)
                seq_match = (
                    fna_df.filter(pl.col("id") == nucleotide_id)
                    if nucleotide_id
                    else fna_df.head(1)
                )
                if seq_match.height > 0:
                    sequence_str = seq_match.row(0, named=True)["sequence"]
                    sequence_bytes = sequence_str.encode("utf-8")
                    sequence_length = len(sequence_bytes)
                else:
                    raise NeighborhoodExtractionError("Nucleotide sequence not found in FNA file")
            else:
                sequence_bytes = b""
                sequence_length = features_df["end"].max() if features_df.height > 0 else 0
        else:
            raise NeighborhoodExtractionError(
                "Must provide either gbf_file or (gff_file and faa_file)"
            )

        start_win, end_win, strand, start_target, end_target = compute_window(
            features_df,
            input_type=input_type,
            mode=mode,
            window=window,
            protein_id=protein_id,
            start=start,
            end=end,
            strand=strand,
            nucleotide_id=nucleotide_id,
            sequence_length=sequence_length,
        )

        subgff = extract_features_in_window(features_df, start_win, end_win)

        if subgff.height == 0:
            raise NeighborhoodExtractionError("No features found in computed window")

        subgff = subgff.with_columns(pl.lit(unique_id).alias("unique_id"))

        if sequence_bytes:
            window_sequence = sequence_bytes[start_win:end_win].decode("utf-8")
        else:
            window_sequence = ""

        if sorfs and window_sequence:
            subgff = annotate_sorfs(
                subgff,
                window_sequence,
                start_win,
                nucleotide_id or record_version,
                unique_id or "unknown",
            )

        neighborhood = pl.DataFrame(
            {
                "seqid": [record_version],
                "start_target": [start_target],
                "end_target": [end_target],
                "start_win": [start_win],
                "end_win": [end_win],
                "strand_win": [strand],
                "sequence": [window_sequence],
                "unique_id": [unique_id or "unknown"],
            }
        )

        return subgff, neighborhood, unique_id, None

    except NeighborhoodExtractionError as e:
        return None, None, unique_id, str(e)
    except Exception as e:
        tb = traceback.format_exc()
        error_msg = f"Unexpected error: {e}\n{tb}"
        return None, None, unique_id, error_msg
