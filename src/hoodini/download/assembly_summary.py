import contextlib
from pathlib import Path

import polars as pl

from hoodini.utils.downloader import download_with_aria2c
from hoodini.utils.logging_utils import logger


def generate_summary_urls(dbs: list[str], include_historical: bool = True) -> list[str]:
    """Generate NCBI assembly summary URLs based on database and historical flag."""
    suffixes = ["", "_historical"] if include_historical else [""]
    return [
        f"https://ftp.ncbi.nlm.nih.gov/genomes/{db}/assembly_summary_{db}{suf}.txt"
        for db in dbs
        for suf in suffixes
    ]


def get_ncbi_header(file_path: Path) -> list[str]:
    """Extract header from an NCBI assembly summary file (line starting with #assembly_accession)."""
    with open(file_path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#assembly_accession"):
                return line.lstrip("#").strip().split("\t")
    raise ValueError(f"No header found in {file_path}")


def download_assembly_db(
    dbs: list[str],
    output_path: Path,
    columns_to_keep: list[str] | None = None,
    include_historical: bool = True,
) -> None:
    """Download multiple assembly summary files using aria2c and save combined table."""
    urls = generate_summary_urls(dbs, include_historical)

    out_names = [Path(url).name for url in urls]
    data_dir = output_path.parent
    logger.info(f"Downloading {len(urls)} assembly summary files with aria2c...")

    result_files = download_with_aria2c(urls, data_dir, show_progress=True, out_names=out_names)

    dfs = []
    for file_path in result_files:
        logger.info(f"Parsing {file_path}")
        try:
            header = get_ncbi_header(file_path)
            df = pl.read_csv(
                file_path,
                separator="\t",
                comment_prefix="#",
                has_header=False,
                new_columns=header,
                quote_char=None,
                null_values="na",
                infer_schema_length=1000,
            )
        except Exception as e:
            logger.error(f"Failed to parse {file_path}: {e}")
            with contextlib.suppress(Exception):
                Path(file_path).unlink()
            continue

        if "assembly_accession" not in df.columns:
            logger.error(f"No 'assembly_accession' in {file_path}!")
            continue

        if columns_to_keep:
            keep = [c for c in columns_to_keep if c in df.columns]
            df = df.select(keep)

        dfs.append(df)
        logger.info(f"Done {file_path}")
        with contextlib.suppress(Exception):
            Path(file_path).unlink()

    if not dfs:
        logger.warning("No dataframes were downloaded.")
        return

    combined = (
        pl.concat(dfs, how="vertical_relaxed").unique(subset=["assembly_accession"]).rechunk()
    )

    combined.write_parquet(output_path)
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
