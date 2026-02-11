"""tRNA and tmRNA gene detection using ARAGORN CLI.

Scans gene neighborhood nucleotide sequences for tRNA and tmRNA genes,
analogous to how ncrna.py uses Infernal/cmsearch for ncRNA detection.
The output DataFrame follows the same schema as ncRNA results so that
both can be merged seamlessly into the GFF and ncrna_metadata outputs.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import polars as pl

from hoodini.utils.logging_utils import info, warn

# Regex for batch-mode gene lines produced by ``aragorn -w -seq``
# Examples:
#   1   tRNA-Arg               [49558,49635]      36      (tct)
#   3   tRNA-???               [136143,136241]     37      (atcc)
#   2   tRNA-Met              c[145324,145397]     35      (cat)
#   1   tmRNA(p)               [10234,10567]      234,300
_GENE_RE = re.compile(
    r"^\s*\d+\s+"  # gene number
    r"(tRNA-\S+|tmRNA\S*)"  # feature name
    r"\s+c?\[(\d+),(\d+)\]"  # optional complement + coords
    r"\s+(\S+)"  # anticodon pos or tag offsets
    r"(?:\s+\((\S+)\))?"  # optional anticodon triplet
)
_COMPLEMENT_RE = re.compile(r"c\[")


def run_trna(
    all_neigh: pl.DataFrame,
    den_data: pl.DataFrame,
    output: str | Path,
    num_threads: int,
    valid_unique_ids: list,
    *,
    translation_table: int = 11,
) -> pl.DataFrame:
    """Run ARAGORN tRNA/tmRNA detection on neighborhood sequences.

    Parameters
    ----------
    all_neigh : pl.DataFrame
        Neighborhood metadata.
    den_data : pl.DataFrame
        Taxonomy / tree metadata (unused, kept for API parity).
    output : str | Path
        Pipeline output directory.
    num_threads : int
        Number of threads (unused by ARAGORN, kept for API parity).
    valid_unique_ids : list
        Unique IDs of non-failed neighbourhoods to include.
    translation_table : int, optional
        Genetic code to use (default 11 = bacterial).

    Returns
    -------
    pl.DataFrame
        DataFrame with columns compatible with ncRNA results.
    """
    if not shutil.which("aragorn"):
        warn(
            "aragorn is not installed.  Install it with: mamba install -c bioconda aragorn\n"
            "Skipping tRNA/tmRNA detection."
        )
        return pl.DataFrame()

    output = Path(output)
    trna_dir = output / "trna"
    trna_dir.mkdir(parents=True, exist_ok=True)

    fasta_path = output / "neighborhood" / "neighborhoods.fasta"
    if not fasta_path.exists():
        warn(f"Neighborhood FASTA not found: {fasta_path}.  Skipping tRNA/tmRNA search.")
        return pl.DataFrame()

    info("🔬\tRunning ARAGORN for tRNA/tmRNA detection...")

    # Run aragorn CLI
    # -t  Search for tRNA genes
    # -m  Search for tmRNA genes
    # -l  Assume linear topology (neighborhoods are linear excerpts)
    # -w  Batch output mode (easy to parse)
    # -seq  Print primary sequence of each gene
    # -gc<n>  Genetic code
    cmd = [
        "aragorn",
        "-t",
        "-m",
        "-l",
        "-w",
        "-seq",
        f"-gc{translation_table}",
        str(fasta_path),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

    if result.returncode != 0:
        warn(f"ARAGORN failed (exit {result.returncode}): {result.stderr[:500]}")
        return pl.DataFrame()

    # Parse batch output
    rows = _parse_batch_output(result.stdout)

    if not rows:
        info("   No tRNA/tmRNA genes found.")
        empty = pl.DataFrame()
        empty.write_csv(trna_dir / "trna_results.tsv", separator="\t", include_header=False)
        return empty

    info(f"   Found {len(rows)} tRNA/tmRNA genes across all neighborhoods")

    # Build DataFrame
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


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------


def _parse_batch_output(stdout: str) -> list[dict]:
    """Parse ARAGORN ``-w -seq`` batch output.

    The format alternates between::

        >seqid
        N genes found
        1  tRNA-Xxx  [start,end]  Apos  (anticodon)
        <sequence across one or more lines>
        ...

    For tmRNA the anticodon field is replaced by tag offsets and a
    tag-peptide line follows.
    """
    rows: list[dict] = []
    current_seqid: str | None = None
    pending_gene: dict | None = None
    collecting_seq = False
    seq_buf: list[str] = []

    for line in stdout.splitlines():
        line_stripped = line.strip()

        # New sequence header
        if line.startswith(">") and not line.startswith(">end"):
            # Flush any pending gene
            if pending_gene is not None:
                pending_gene["sequence"] = "".join(seq_buf)
                rows.append(pending_gene)
                pending_gene = None
                seq_buf = []
                collecting_seq = False
            current_seqid = line_stripped.lstrip(">").split()[0]
            continue

        # Summary line ("N genes found" or ">end ...")
        if "genes found" in line_stripped or "gene found" in line_stripped:
            continue
        if line_stripped.startswith(">end"):
            # Flush last gene
            if pending_gene is not None:
                pending_gene["sequence"] = "".join(seq_buf)
                rows.append(pending_gene)
                pending_gene = None
                seq_buf = []
                collecting_seq = False
            continue

        # Nothing found line
        if "nothing found" in line_stripped.lower():
            continue

        # Try matching a gene line
        m = _GENE_RE.match(line_stripped)
        if m and current_seqid is not None:
            # Flush previous gene if any
            if pending_gene is not None:
                pending_gene["sequence"] = "".join(seq_buf)
                rows.append(pending_gene)
                seq_buf = []

            feature_raw = m.group(1)  # e.g. "tRNA-Arg", "tRNA-???", "tmRNA", "tmRNA(p)"
            start = int(m.group(2))
            end = int(m.group(3))
            anticodon = m.group(5) or ""
            is_complement = bool(_COMPLEMENT_RE.search(line_stripped))
            strand = "-" if is_complement else "+"

            # Build feature name
            if feature_raw.startswith("tRNA"):
                nc_feature = f"{feature_raw}({anticodon})" if anticodon else feature_raw
            else:
                nc_feature = feature_raw

            pending_gene = {
                "nucid": current_seqid,
                "nc_feature": nc_feature,
                "seqfrom": start,
                "seqto": end,
                "strand_ncrna": strand,
                "score": ".",
                "E-value": ".",
                "sequence": "",
                "structure": "",
            }
            collecting_seq = True
            seq_buf = []
            continue

        # tmRNA tag peptide line (e.g. "ANDNFAEEFAV*")
        if (
            pending_gene is not None
            and pending_gene["nc_feature"].startswith("tmRNA")
            and line_stripped
            and re.match(r"^[A-Z*]+$", line_stripped)
        ):
            pending_gene["tag_peptide"] = line_stripped.rstrip("*")
            continue

        # Sequence lines (lowercase nucleotides)
        if (
            collecting_seq
            and line_stripped
            and re.match(r"^[acgturykmswbdhvnACGTURYKMSWBDHVN\s]+$", line_stripped)
        ):
            seq_buf.append(line_stripped.replace(" ", ""))
            continue

    # Flush final pending gene
    if pending_gene is not None:
        pending_gene["sequence"] = "".join(seq_buf)
        rows.append(pending_gene)

    return rows
