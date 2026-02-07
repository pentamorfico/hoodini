"""Input validation and ID parsing utilities."""

from __future__ import annotations

import csv
import re
from pathlib import Path

import polars as pl
import rich_click as click

from hoodini.utils.id_parsing import categorize_id
from hoodini.utils.logging_utils import info, warn

# Reserved columns used internally by the pipeline - these should not be overwritten by user data
RESERVED_COLUMNS = {
    # Input identification
    "og_index",
    "unique_id",
    "protein_id",
    "nucleotide_id",
    "uniprot_id",
    "input_type",
    # Paths
    "gff_path",
    "faa_path",
    "fna_path",
    "gbf_path",
    # Coordinates
    "start",
    "end",
    "strand",
    # Assembly/taxonomy
    "taxid",
    "assembly_id",
    # Status
    "failed",
    "failed_reason",
    "premade",
    "is_full_contig",
    # Query info (added by pipeline)
    "query_protein_id",
    "is_refseq_query",
    "sequence_length",
    "group",
    "species_taxid",
    "organism_name",
    "infraspecific_name",
    "assembly_level",
    "nucleotide_id_no_prefix",
    # DSMZ dive columns
    "dive_id",
    "collection_id",
    "dive_type",
}


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

                # Allow alphanumeric, dots, underscores, colons and hyphens for NucID:start-end format
                special_char_pattern = re.compile(r"[^A-Za-z0-9\s._:\-]+")
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
                fieldnames = reader.fieldnames or []

                # Check for at least one ID column
                id_columns = ["nucleotide_id", "protein_id", "uniprot_id"]
                has_id_col = any(col in fieldnames for col in id_columns)
                if not has_id_col:
                    raise click.BadParameter(
                        "TSV file must contain at least one ID column: nucleotide_id, protein_id, or uniprot_id"
                    )

                # Check for local file columns
                has_gbf_col = "gbf_path" in fieldnames
                has_gff_faa_cols = "gff_path" in fieldnames and "faa_path" in fieldnames

                found_valid_row = False
                invalid_ids = []

                for row_num, row in enumerate(reader, start=2):  # start=2 because header is row 1
                    # Check if row has local files
                    has_local_gbf = has_gbf_col and row.get("gbf_path", "").strip()
                    has_local_gff_faa = (
                        has_gff_faa_cols
                        and row.get("gff_path", "").strip()
                        and row.get("faa_path", "").strip()
                    )
                    has_local_files = has_local_gbf or has_local_gff_faa

                    # Get the ID value from any ID column
                    id_value = None
                    for col in id_columns:
                        if col in row and row[col] and row[col].strip():
                            id_value = row[col].strip()
                            found_valid_row = True
                            break

                    # If no local files, validate the ID format
                    if id_value and not has_local_files:
                        result = categorize_id(id_value)
                        if result["type"] == "unmatched":
                            invalid_ids.append(f"Row {row_num}: '{id_value}'")

                if not found_valid_row:
                    raise click.BadParameter(
                        "The TSV file must contain at least one valid row with an ID value."
                    )

                if invalid_ids:
                    raise click.BadParameter(
                        "Invalid ID format(s) without local files - must be valid NCBI or UniProt ID:\n  "
                        + "\n  ".join(invalid_ids[:5])
                        + (
                            f"\n  ... and {len(invalid_ids) - 5} more"
                            if len(invalid_ids) > 5
                            else ""
                        )
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
    df = pl.read_csv(filename, separator="\t", infer_schema_length=0)
    df = df.with_row_index("og_index").with_columns(pl.col("og_index").cast(pl.Utf8))

    # Identify user-provided extra columns (not in reserved set)
    user_extra_columns = [col for col in df.columns if col not in RESERVED_COLUMNS]

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
        "failed_reason",
        "input_type",
        "premade",
        "start",
        "end",
        "strand",
    ]
    for col in expected_columns:
        if col not in df.columns:
            df = df.with_columns(pl.lit(None).alias(col))

    # Determine input_type if not explicitly set
    # Priority: nucleotide_id > protein_id > uniprot_id
    df = df.with_columns(
        pl.when(pl.col("input_type").is_not_null() & (pl.col("input_type") != ""))
        .then(pl.col("input_type"))
        .when(pl.col("nucleotide_id").is_not_null() & (pl.col("nucleotide_id") != ""))
        .then(pl.lit("nucleotide"))
        .when(pl.col("protein_id").is_not_null() & (pl.col("protein_id") != ""))
        .then(pl.lit("protein"))
        .when(pl.col("uniprot_id").is_not_null() & (pl.col("uniprot_id") != ""))
        .then(pl.lit("protein"))
        .otherwise(pl.lit(None))
        .alias("input_type")
    )

    # Store user extra columns as metadata for downstream propagation
    if user_extra_columns:
        # Add a hidden attribute to track extra columns (will be used by write_data)
        df = df.with_columns(pl.lit(",".join(user_extra_columns)).alias("_user_extra_cols"))

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
            record["input_type"] = "nucleotide"

            # Check if this is a NucID:start-end format
            # categorize_id puts the "start-end" part in protein_id when it sees ":"
            potential_range = category.get("protein_id")
            if potential_range and "-" in potential_range:
                parts = potential_range.split("-")
                if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
                    record["start"] = int(parts[0])
                    record["end"] = int(parts[1])
                    if record["start"] > record["end"]:
                        record["start"], record["end"] = record["end"], record["start"]
                        record["strand"] = "-"
                    else:
                        record["strand"] = "+"
                    record["protein_id"] = None  # Clear the misinterpreted protein_id
                else:
                    record["protein_id"] = potential_range
            else:
                record["protein_id"] = potential_range

        elif category["type"] == "uniprot":
            record["uniprot_id"] = category["id"]
            record["input_type"] = "protein"
        elif category["type"] == "uniparc":
            record["uniprot_id"] = category["id"]  # stored here, resolved later
            record["input_type"] = "uniparc"
        elif category["type"] == "unmatched":
            record["failed"] = True
            record["failed_reason"] = "not valid ID"
        elif category.get("id"):
            # Fallback for unknown types
            record["protein_id"] = category["id"]
            record["input_type"] = "protein"
        else:
            record["failed"] = True
            record["failed_reason"] = "not valid ID"

        data.append(record)
    df = pl.DataFrame(data, infer_schema_length=len(data))
    return df


def uniparc2ncbi(df: pl.DataFrame) -> pl.DataFrame:
    """
    Resolve UniParc IDs (UPI...) directly to NCBI protein IDs.

    Queries the UniProt REST API (``/uniparc/{upi}``) and extracts protein IDs
    from active cross-references:

    1. Prefers **RefSeq** (e.g. ``YP_232970``)
    2. Falls back to **EMBL** (e.g. ``AAO89367``) if no active RefSeq entry

    Sets ``protein_id`` directly — no intermediate UniProt KB lookup needed.
    """
    if "input_type" not in df.columns or "uniprot_id" not in df.columns:
        return df

    mask = df["input_type"] == "uniparc"
    if mask.sum() == 0:
        return df

    upi_ids = df.filter(mask)["uniprot_id"].drop_nulls().unique().to_list()
    if not upi_ids:
        return df

    import asyncio

    import httpx

    info(f"🔎 Resolving {len(upi_ids)} UniParc ID(s) via UniProt REST API (parallel)...")

    upi_to_protein: dict[str, str] = {}
    base_url = "https://rest.uniprot.org/uniparc"
    MAX_CONCURRENT = 10  # be nice to the API

    def _pick_protein_id(cross_refs: list[dict]) -> str | None:
        """Pick the best protein ID from cross-references: RefSeq > EMBL.

        Appends the version number (e.g. ``WP_292220451.1``) when available,
        since downstream IPG lookups require versioned accessions.
        """

        def _versioned(xref: dict) -> str:
            pid = xref["id"]
            ver = xref.get("version")
            if ver is not None and "." not in pid:
                return f"{pid}.{ver}"
            return pid

        for xref in cross_refs:
            if xref.get("database") == "RefSeq" and xref.get("active") is True:
                return _versioned(xref)
        for xref in cross_refs:
            if xref.get("database") == "EMBL" and xref.get("active") is True:
                return _versioned(xref)
        return None

    async def _fetch_one(
        client: httpx.AsyncClient,
        sem: asyncio.Semaphore,
        upi: str,
    ) -> tuple[str, str | None]:
        async with sem:
            try:
                resp = await client.get(f"{base_url}/{upi}", params={"format": "json"})
                resp.raise_for_status()
                entry = resp.json()
                xrefs = entry.get("uniParcCrossReferences", [])
                return upi, _pick_protein_id(xrefs)
            except Exception as e:
                warn(f"⚠️  Failed to resolve {upi}: {e}")
                return upi, None

    async def _resolve_all() -> dict[str, str]:
        sem = asyncio.Semaphore(MAX_CONCURRENT)
        async with httpx.AsyncClient(timeout=30) as client:
            tasks = [_fetch_one(client, sem, upi) for upi in upi_ids]
            results = await asyncio.gather(*tasks)
        return {upi: pid for upi, pid in results if pid is not None}

    try:
        loop = asyncio.get_running_loop()
    except RuntimeError:
        loop = None

    if loop and loop.is_running():
        import nest_asyncio

        nest_asyncio.apply()

    upi_to_protein = asyncio.run(_resolve_all())

    if not upi_to_protein:
        warn("Could not resolve any UniParc IDs. Marking them as failed.")
        return df.with_columns(
            pl.when(mask).then(pl.lit(True)).otherwise(pl.col("failed")).alias("failed"),
            pl.when(mask)
            .then(pl.lit("UniParc ID could not be resolved."))
            .otherwise(pl.col("failed_reason"))
            .alias("failed_reason"),
        )

    # Build mapping frame
    map_df = pl.DataFrame(
        {"upi": list(upi_to_protein.keys()), "resolved_id": list(upi_to_protein.values())}
    )
    df = df.join(map_df, left_on="uniprot_id", right_on="upi", how="left")

    # Set protein_id directly and update input_type
    df = df.with_columns(
        pl.when(mask & pl.col("resolved_id").is_not_null())
        .then(pl.col("resolved_id"))
        .otherwise(pl.col("protein_id"))
        .alias("protein_id"),
        pl.when(mask & pl.col("resolved_id").is_not_null())
        .then(pl.lit("protein"))
        .otherwise(pl.col("input_type"))
        .alias("input_type"),
    ).drop("resolved_id")

    # Mark unresolved UniParc IDs as failed
    still_uniparc = df["input_type"] == "uniparc"
    if still_uniparc.sum() > 0:
        df = df.with_columns(
            pl.when(still_uniparc).then(pl.lit(True)).otherwise(pl.col("failed")).alias("failed"),
            pl.when(still_uniparc)
            .then(pl.lit("UniParc ID could not be resolved to RefSeq or EMBL."))
            .otherwise(pl.col("failed_reason"))
            .alias("failed_reason"),
        )

    resolved = len(upi_to_protein)
    info(f"✔️  Resolved {resolved}/{len(upi_ids)} UniParc IDs to protein IDs.")
    return df


def uniprot2ncbi(df: pl.DataFrame) -> pl.DataFrame:
    """
    Map UniProt accessions to NCBI protein IDs using a local parquet database.

    Uses ``idmapping_selected.parquet`` (columns: UniprotKB-AC, RefSeq, EMBL-CDS).
    Prefers the RefSeq mapping; falls back to EMBL-CDS when RefSeq is null.

    - Only attempts mapping for rows with `uniprot_id` present, `protein_id` null/empty,
      and no local files (`gff_path`/`faa_path` missing).
    - On failure (no match found), sets the `failed` column with a descriptive message.
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

    to_map = df.filter(mask)["uniprot_id"].drop_nulls().unique().to_list()
    if not to_map:
        return df

    from hoodini.download.idmapping import get_idmapping_path

    idmap_path = get_idmapping_path()
    if not idmap_path.exists():
        warn(
            "ID-mapping database not found. "
            "Run 'hoodini download databases' first. Marking UniProt entries as failed."
        )
        return df.with_columns(
            pl.when(mask).then(pl.lit(True)).otherwise(pl.col("failed")).alias("failed"),
            pl.when(mask)
            .then(pl.lit("ID-mapping database not available."))
            .otherwise(pl.col("failed_reason"))
            .alias("failed_reason"),
        )

    # Use DuckDB for memory-efficient semi-join against the parquet
    import duckdb

    con = duckdb.connect(":memory:")
    con.execute('SET memory_limit = "2GB"')

    # Register the lookup IDs as a temp table
    con.execute("CREATE TEMP TABLE lookup_ids (uniprot_ac VARCHAR)")
    con.executemany("INSERT INTO lookup_ids VALUES (?)", [(uid,) for uid in to_map])

    idmap = con.execute(f"""
        SELECT
            m."UniprotKB-AC",
            CASE
                WHEN m."RefSeq" IS NOT NULL AND TRIM(m."RefSeq") != ''
                    THEN m."RefSeq"
                WHEN m."EMBL-CDS" IS NOT NULL AND TRIM(m."EMBL-CDS") != ''
                    THEN m."EMBL-CDS"
                ELSE NULL
            END AS mapped_id
        FROM read_parquet('{str(idmap_path)}') m
        SEMI JOIN lookup_ids l ON m."UniprotKB-AC" = l.uniprot_ac
        WHERE mapped_id IS NOT NULL
        """).pl()

    con.close()

    # Keep first match per UniProt accession
    idmap = idmap.unique(subset=["UniprotKB-AC"], keep="first")

    # Join mapping onto the records
    df = df.join(
        idmap.select("UniprotKB-AC", "mapped_id"),
        left_on="uniprot_id",
        right_on="UniprotKB-AC",
        how="left",
    )

    # Fill in protein_id from mapped_id where it was missing
    df = df.with_columns(
        pl.when(mask & pl.col("mapped_id").is_not_null())
        .then(pl.col("mapped_id"))
        .otherwise(pl.col("protein_id"))
        .alias("protein_id")
    ).drop("mapped_id")

    # Mark rows that had a uniprot_id but could not be mapped
    failed_mask = mask & df["protein_id"].cast(pl.Utf8, strict=False).is_null()
    df = df.with_columns(
        pl.when(failed_mask).then(pl.lit(True)).otherwise(pl.col("failed")).alias("failed"),
        pl.when(failed_mask)
        .then(pl.lit("No associated NCBI ID found for the UniProt entry."))
        .otherwise(pl.col("failed_reason"))
        .alias("failed_reason"),
    )

    return df
