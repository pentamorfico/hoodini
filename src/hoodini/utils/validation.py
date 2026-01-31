"""Input validation and ID parsing utilities."""

from __future__ import annotations

import csv
import re
from pathlib import Path

import polars as pl
import rich_click as click

from hoodini.utils.id_parsing import categorize_id
from hoodini.utils.logging_utils import warn


def validate_literal_id(query: str) -> None:
    """
    Validate that a literal query (not a file) is a valid ID format.

    Accepts:
    - NCBI protein IDs: WP_*, NP_*, XP_*, YP_*, ZP_*, or 3-letter + 5-8 digits
    - NCBI nucleotide IDs: NC_*, NZ_*, etc.
    - UniProt IDs: e.g., P12345, Q9Y6K9
    - FASTA sequences (start with > or are pure amino acid sequences)

    Raises:
        ValueError: If the ID format is not recognized
    """
    query = query.strip()

    # Allow FASTA sequences
    if query.startswith(">"):
        return

    # Allow pure amino acid sequences (all uppercase letters, possibly with *)
    if re.fullmatch(r"[A-Z*]+", query.replace("\n", ""), re.I):
        return

    # Check against known ID patterns
    result = categorize_id(query)
    if result["type"] == "unmatched":
        raise ValueError(
            f"Unrecognized ID format: '{query}'. "
            "Expected NCBI protein ID (e.g., WP_000000001, NP_414542), "
            "UniProt ID (e.g., P12345), nucleotide ID (e.g., NC_000913), "
            "or a FASTA sequence."
        )


def validate_input_file(ctx, param, value):
    """Validate the input file based on its type, either single-column or TSV format."""
    if value is None:
        return None

    if not Path(value).is_file():
        # Allow literal IDs/FASTA strings; runner will handle them.
        return value

    try:
        with open(value) as file:
            if param.name == "input_path":
                lines = [line.strip() for line in file.readlines() if line.strip()]
                if len(lines) <= 1:
                    raise click.BadParameter("File must contain multiple lines.")

                for line in lines:
                    if "," in line or "\t" in line:
                        raise click.BadParameter(
                            "File must be a single-column text file without delimiters like commas or tabs."
                        )

                special_char_pattern = re.compile(r"[^A-Za-z0-9\s._]+")
                for i, line in enumerate(lines):
                    match = special_char_pattern.search(line)
                    if match:
                        raise click.BadParameter(
                            f"Invalid character '{match.group()}' found in line {i + 1}: \"{line}\""
                        )

            elif param.name == "inputsheet":
                first_line = file.readline()
                delimiter = "\t" if "\t" in first_line else None
                if delimiter is None:
                    raise click.BadParameter("The file is not in TSV format.")
                file.seek(0)
                reader = csv.DictReader(file, delimiter=delimiter)
                required_columns = [
                    "nucleotide_id",
                    "protein_id",
                    "gff_path",
                    "fna_path",
                    "faa_path",
                ]
                if not all(col in reader.fieldnames for col in required_columns):
                    raise click.BadParameter(
                        "TSV file must contain columns: nucleotide_id, protein_id, gff_path, fna_path, faa_path"
                    )

                found_valid_row = False
                for row in reader:
                    if row["nucleotide_id"].strip() or row["protein_id"].strip():
                        found_valid_row = True
                        break
                if not found_valid_row:
                    raise click.BadParameter(
                        "The TSV file must contain at least one valid row with required data."
                    )

    except Exception as e:
        raise click.BadParameter(f"Error reading file: {e}")

    return value


def validate_domains(ctx, param, value):
    """Validate domain database names against MetaCerberus availability."""
    if not value or (isinstance(value, str) and value.strip() == ""):
        return None

    db_names = [d.strip().lower() for d in value.split(",") if d.strip()]
    if not db_names:
        raise click.BadParameter("Domain parameter must contain at least one database name.")

    for db in db_names:
        if not re.match(r"^[a-zA-Z0-9_-]+$", db):
            raise click.BadParameter(
                f"Invalid database name '{db}'. Only letters, numbers, underscores, and hyphens are allowed."
            )

    try:
        from hoodini.download.metacerberus import check_downloaded, get_db_groups, list_db_files

        files_list = list_db_files()
        groups = get_db_groups(files_list)
        status = check_downloaded(groups)

        valid_dbs = []
        invalid_dbs = []
        missing_files_dbs = []

        for db in db_names:
            if db not in groups:
                invalid_dbs.append(db)
                continue

            file_statuses = status.get(db, [])
            hmm_present = any(
                present for f, present in file_statuses if f["name"].endswith(".hmm.gz")
            )
            tsv_present = any(present for f, present in file_statuses if f["name"].endswith(".tsv"))

            has_hmm_file = any(f["name"].endswith(".hmm.gz") for f, _ in file_statuses)
            has_tsv_file = any(f["name"].endswith(".tsv") for f, _ in file_statuses)

            if has_tsv_file and not tsv_present:
                missing_files_dbs.append(f"{db} (missing TSV)")
            elif has_hmm_file and not hmm_present:
                missing_files_dbs.append(f"{db} (missing HMM)")
            elif not has_hmm_file and not has_tsv_file:
                invalid_dbs.append(db)
            elif (has_hmm_file and hmm_present and has_tsv_file and tsv_present) or (
                not has_hmm_file and has_tsv_file and tsv_present
            ):
                valid_dbs.append(db)
            else:
                missing_files_dbs.append(db)

        if invalid_dbs:
            raise click.BadParameter(
                f"Unknown MetaCerberus databases: {', '.join(invalid_dbs)}. Run 'hoodini download metacerberus' to see available databases."
            )

        if missing_files_dbs:
            raise click.BadParameter(
                f"MetaCerberus databases not downloaded: {', '.join(missing_files_dbs)}. Run 'hoodini download metacerberus {','.join(db_names)}' to download them."
            )

        if not valid_dbs:
            raise click.BadParameter("No valid MetaCerberus databases found.")

        return valid_dbs

    except ImportError as e:
        warn(f"Could not validate MetaCerberus databases: {e}")
        return value
    except Exception as e:
        warn(f"Could not validate MetaCerberus databases: {e}")
        return value


def switch_assembly_prefix(asm_id):
    if not isinstance(asm_id, str):
        return asm_id
    if asm_id.startswith("GCA_"):
        return "GCF_" + asm_id[4:]
    if asm_id.startswith("GCF_"):
        return "GCA_" + asm_id[4:]
    return asm_id


def is_refseq_nuccore(nuc_id):
    """Return True if the nuccore accession is a RefSeq accession else False."""
    refseq_prefixes = ("NC_", "NZ_", "NM_", "NR_", "XM_", "XR_", "AP_", "YP_", "XP_", "WP_")
    return isinstance(nuc_id, str) and nuc_id.startswith(refseq_prefixes)


def read_input_sheet(filename):
    df = pl.read_csv(filename, separator="\t", dtype=str)
    df = df.with_row_count("og_index").with_columns(pl.col("og_index").cast(pl.Utf8))

    expected_columns = [
        "protein_id",
        "nucleotide_id",
        "uniprot_id",
        "gff_path",
        "faa_path",
        "fna_path",
        "gbf_path",
        "taxid",
        "assembly_id",
        "failed",
        "input_type",
        "premade",
    ]
    for col in expected_columns:
        if col not in df.columns:
            df = df.with_columns(pl.lit(None).alias(col))
    return df


def read_input_list(filename):
    input_list = Path(filename).read_text().splitlines()
    data = []
    for index, id_ in enumerate(input_list):
        if not id_ or id_.strip() == "":
            continue

        category = categorize_id(id_)
        record = {
            "og_index": index,
            "protein_id": None,
            "nucleotide_id": None,
            "uniprot_id": None,
            "failed": None,
            "failed_reason": None,
            "gff_path": None,
            "faa_path": None,
            "fna_path": None,
            "strand": None,
            "start": None,
            "end": None,
            "gbf_path": None,
            "taxid": None,
            "assembly_id": None,
            "input_type": None,
            "premade": None,
        }

        if category["type"] == "protein":
            record["protein_id"] = category["id"]
            record["input_type"] = "protein"
        elif category["type"] == "nucleotide":
            record["nucleotide_id"] = category["id"]
            record["protein_id"] = category.get("protein_id")
            record["input_type"] = "nucleotide"
        elif category["type"] == "uniprot":
            record["uniprot_id"] = category["id"]
            record["input_type"] = "protein"
        elif category["type"] == "unmatched":
            record["failed"] = True
            record["failed_reason"] = "not valid ID"
        else:
            if ":" in id_:
                category = categorize_id(id_.split(":")[0])
                if category["type"] == "nucleotide":
                    if (
                        "-" in id_.split(":")[1]
                        and id_.split(":")[1].split("-")[0].isdigit()
                        and id_.split(":")[1].split("-")[1].isdigit()
                    ):
                        record["nucleotide_id"] = id_.split(":")[0]
                        record["start"] = id_.split(":")[1].split("-")[0]
                        record["end"] = id_.split(":")[1].split("-")[1]
                        if record["start"] > record["end"]:
                            record["start"], record["end"] = record["end"], record["start"]
                            record["strand"] = "-"
                        else:
                            record["strand"] = "+"
                        record["input_type"] = "nucleotide"
                    else:
                        record["failed"] = True
                        record["failed_reason"] = "not valid ID"
                    record["nucleotide_id"] = category["id"]
                    record["input_type"] = "nucleotide"

            if record["protein_id"] is None and category.get("id"):
                record["protein_id"] = category["id"]
                record["input_type"] = record.get("input_type") or "protein"
        data.append(record)
    df = pl.DataFrame(data, infer_schema_length=len(data))
    return df


def uniprot2ncbi(df: pl.DataFrame) -> pl.DataFrame:
    """
    Map UniProt accessions to NCBI protein IDs using UniProtMapper.

    - Only attempts mapping for rows with `uniprot_id` present, `protein_id` null/empty,
      and no local files (`gff_path`/`faa_path` missing).
    - If mapping returns multiple hits, merged `To` values are applied directly.
    - On failure, sets the `failed` column with a descriptive message.
    """
    required_cols = {"uniprot_id", "protein_id", "gff_path", "faa_path"}
    if not required_cols.issubset(set(df.columns)):
        return df

    # Cast protein_id to Utf8 to handle Null-type columns safely
    protein_id_str = df["protein_id"].cast(pl.Utf8, strict=False)
    mask = (
        df["uniprot_id"].is_not_null()
        & (protein_id_str.is_null() | (protein_id_str == ""))
        & df["gff_path"].is_null()
        & df["faa_path"].is_null()
    )

    if mask.sum() == 0:
        return df

    from UniProtMapper import ProtMapper

    df_pd = df.to_pandas()
    to_map = df_pd.loc[mask.to_pandas(), "uniprot_id"].dropna().unique().tolist()
    if not to_map:
        return df

    mapper = ProtMapper()
    try:
        mapped_df, failed_ids = mapper.get(
            ids=to_map, from_db="UniProtKB_AC-ID", to_db="EMBL-GenBank-DDBJ_CDS"
        )
    except Exception:
        failed_ids = to_map
        mapped_df = None

    if mapped_df is not None and not mapped_df.empty:
        df_pd = df_pd.merge(
            mapped_df[["From", "To"]], left_on="uniprot_id", right_on="From", how="left"
        )
        df_pd["protein_id"] = df_pd["protein_id"].fillna(df_pd["To"])
        df_pd = df_pd.drop(columns=["From", "To"])

    if failed_ids:
        failed_mask = df_pd["uniprot_id"].isin(failed_ids)
        df_pd.loc[failed_mask, "failed"] = True
        df_pd.loc[failed_mask, "failed_reason"] = "No associated NCBI found for the UniProt entry."

    return pl.from_pandas(df_pd)
