import shutil
from pathlib import Path
from tempfile import mkdtemp

import polars as pl

from hoodini.utils.downloader import download_urls
from hoodini.utils.logging_utils import logger

# NCBI fields are typed by name, never by their position or first few values.
# All other fields (including future additions) remain strings. In particular,
# PubMed IDs may contain semicolon-separated lists and isolates can start with 0.
NUMERIC_COLUMNS = dict.fromkeys(
    (
        "taxid",
        "species_taxid",
        "genome_size",
        "genome_size_ungapped",
        "replicon_count",
        "scaffold_count",
        "contig_count",
        "total_gene_count",
        "protein_coding_gene_count",
        "non_coding_gene_count",
    ),
    pl.Int64,
)
NUMERIC_COLUMNS["gc_percent"] = pl.Float64


def generate_summary_urls(dbs: list[str], include_historical: bool = True) -> list[str]:
    """Generate NCBI assembly summary URLs based on database and historical flag."""
    suffixes = ["", "_historical"] if include_historical else [""]
    return [
        f"https://ftp.ncbi.nlm.nih.gov/genomes/{db}/assembly_summary_{db}{suf}.txt"
        for db in dbs
        for suf in suffixes
    ]


def get_ncbi_header(file_path: Path) -> list[str]:
    """Find the named TSV header, accepting spacing and column-order changes."""
    with open(file_path, encoding="utf-8-sig") as f:
        for line in f:
            if not line.startswith("#"):
                continue
            header = [value.strip() for value in line[1:].rstrip("\r\n").split("\t")]
            if "assembly_accession" in header:
                if not all(header) or len(header) != len(set(header)):
                    raise ValueError(f"Empty or duplicate header fields in {file_path}")
                return header
    raise ValueError(f"No header found in {file_path}")


def _validate_row_widths(file_path: Path, expected: int) -> None:
    """Reject ragged TSV rows before the CSV reader can fill missing fields.

    NCBI summary files use literal tabs without CSV quoting. This streaming pass
    keeps memory bounded and gives physical line numbers for malformed rows.
    """
    rows = 0
    with file_path.open("rb") as source:
        for line_number, line in enumerate(source, start=1):
            if line_number == 1:
                line = line.removeprefix(b"\xef\xbb\xbf")
            if line.startswith(b"#") or not line.strip():
                continue
            fields = line.count(b"\t") + 1
            if fields != expected:
                raise ValueError(
                    f"{file_path}: line {line_number} has {fields} fields; "
                    f"the header declares {expected}. Check for a truncated or shifted row."
                )
            rows += 1
    if not rows:
        raise ValueError(f"No assembly records found in {file_path}")


def read_assembly_summary(
    file_path: Path, columns_to_keep: list[str] | None = None
) -> pl.DataFrame:
    """Read a summary with an explicit schema and validate accession identities."""
    header = get_ncbi_header(file_path)
    _validate_row_widths(file_path, len(header))
    schema = {name: NUMERIC_COLUMNS.get(name, pl.Utf8) for name in header}
    try:
        frame = pl.read_csv(
            file_path,
            separator="\t",
            comment_prefix="#",
            has_header=False,
            schema=schema,
            quote_char=None,
            null_values="na",
        )
    except pl.exceptions.PolarsError as exc:
        raise ValueError(f"Invalid assembly summary {file_path}: {exc}") from exc

    valid_ids = pl.col("assembly_accession").str.contains(r"^GC[AF]_\d+\.\d+$").fill_null(False)
    if frame.filter(~valid_ids).height:
        raise ValueError(f"Invalid or missing assembly_accession in {file_path}")

    if columns_to_keep:
        frame = frame.select(
            [
                (
                    pl.col(name)
                    if name in frame.columns
                    else pl.lit(None).cast(NUMERIC_COLUMNS.get(name, pl.Utf8)).alias(name)
                )
                for name in columns_to_keep
            ]
        )
    return frame


def download_assembly_db(
    dbs: list[str],
    output_path: Path,
    columns_to_keep: list[str] | None = None,
    include_historical: bool = True,
) -> None:
    """Validate every requested source before atomically replacing the database.

    Each attempt downloads into its own directory on the destination filesystem.
    Failed downloads/inputs remain there for diagnosis, with the path logged.
    """
    urls = generate_summary_urls(dbs, include_historical)
    if not urls:
        raise ValueError("At least one assembly summary database must be requested")
    if columns_to_keep and "assembly_accession" not in columns_to_keep:
        raise ValueError("columns_to_keep must include assembly_accession")

    out_names = [Path(url).name for url in urls]
    data_dir = output_path.parent
    data_dir.mkdir(parents=True, exist_ok=True)
    staging_dir = Path(mkdtemp(prefix="assembly-summary-update-", dir=data_dir))
    try:
        logger.info(f"Downloading {len(urls)} assembly summary files...")
        result_files = download_urls(urls, staging_dir, show_progress=True, out_names=out_names)
        downloaded = {Path(path).resolve() for path in result_files}
        expected = [staging_dir / name for name in out_names]
        missing = [
            path.name for path in expected if path.resolve() not in downloaded or not path.is_file()
        ]
        if missing:
            raise FileNotFoundError(f"Missing requested assembly summaries: {', '.join(missing)}")

        dfs = []
        for file_path in expected:
            logger.info(f"Parsing {file_path}")
            dfs.append(read_assembly_summary(file_path, columns_to_keep))

        combined = (
            pl.concat(dfs, how="diagonal_relaxed")
            .unique(subset=["assembly_accession"], keep="first", maintain_order=True)
            .rechunk()
        )
        staged_output = staging_dir / output_path.name
        combined.write_parquet(staged_output)
        staged_output.replace(output_path)
    except BaseException:
        logger.error(
            f"Assembly summary update failed; destination was not replaced. "
            f"Downloaded files retained in {staging_dir}"
        )
        raise

    try:
        shutil.rmtree(staging_dir)
    except OSError as exc:
        logger.warning(f"Database updated, but could not remove {staging_dir}: {exc}")
    logger.info(f"Saved combined parquet to {output_path}")


def download_assembly_summary_db(output_path: Path | None = None) -> Path:
    """Convenience wrapper to download and merge RefSeq + GenBank assembly summaries."""
    from importlib.resources import files

    if output_path is None:
        output_path = files("hoodini").joinpath("data", "assembly_summary.parquet")

    columns = [
        "assembly_accession",
        "refseq_category",
        "taxid",
        "species_taxid",
        "organism_name",
        "infraspecific_name",
        "isolate",
        "assembly_level",
        "genome_rep",
        "gbrs_paired_asm",
        "paired_asm_comp",
        "group",
        "ftp_path",
        "genome_size",
        "gc_percent",
        "replicon_count",
        "scaffold_count",
        "contig_count",
        "seq_rel_date",
    ]

    output_path.parent.mkdir(parents=True, exist_ok=True)

    download_assembly_db(
        dbs=["refseq", "genbank"],
        output_path=output_path,
        columns_to_keep=columns,
        include_historical=True,
    )

    return output_path
