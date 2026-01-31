import re
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import polars as pl
import requests

from hoodini.utils.logging_utils import error, info, warn

# Regex pattern for RFAM IDs: RF followed by exactly 5 digits
RFAM_PATTERN = re.compile(r"^RF\d{5}$")


def is_rfam_id(value: str) -> bool:
    """Check if a string is a valid RFAM ID (RF followed by 5 digits)."""
    return bool(RFAM_PATTERN.match(value))


def validate_ncrna_input(ncrna_input: str) -> None:
    """
    Validate --ncrna input early (before pipeline runs).

    Raises:
        ValueError: If RFAM ID format is invalid
        FileNotFoundError: If CM file path doesn't exist
    """
    parts = [p.strip() for p in ncrna_input.split(",")]

    # If all parts are valid RFAM IDs, we're good
    if all(is_rfam_id(p) for p in parts):
        return

    # Check for invalid RFAM ID format
    invalid_rfam = [p for p in parts if p.upper().startswith("RF") and not is_rfam_id(p)]
    if invalid_rfam:
        raise ValueError(
            f"Invalid RFAM ID format: {', '.join(invalid_rfam)}. "
            "RFAM IDs must be 'RF' followed by exactly 5 digits (e.g., RF00001, RF02348)"
        )

    # Treat as file path - check existence
    path = Path(ncrna_input)
    if not path.exists():
        raise FileNotFoundError(f"ncRNA CM file not found: {path}")


def parse_ncrna_input(ncrna_input: str) -> tuple[bool, list[str] | Path]:
    """
    Parse the --ncrna input to determine if it's a path or RFAM IDs.

    Returns:
        Tuple of (is_rfam_ids, value) where:
        - is_rfam_ids: True if input contains RFAM IDs, False if it's a path
        - value: list of RFAM IDs or Path to CM file
    """
    # Check if it's a comma-separated list of RFAM IDs
    parts = [p.strip() for p in ncrna_input.split(",")]

    # If all parts are valid RFAM IDs, treat as RFAM input
    if all(is_rfam_id(p) for p in parts):
        return True, parts

    # Check if any part looks like an invalid RFAM ID (starts with RF but wrong format)
    invalid_rfam = [p for p in parts if p.upper().startswith("RF") and not is_rfam_id(p)]
    if invalid_rfam:
        error(f"Invalid RFAM ID format: {', '.join(invalid_rfam)}")
        error("RFAM IDs must be 'RF' followed by exactly 5 digits (e.g., RF00001, RF02348)")
        raise ValueError(
            f"Invalid RFAM ID format: {', '.join(invalid_rfam)}. Expected format: RF##### (5 digits)"
        )

    # Otherwise treat as a file path
    path = Path(ncrna_input)
    if not path.exists():
        error(f"ncRNA CM file not found: {path}")
        raise FileNotFoundError(f"CM file not found: {path}")
    return False, path


def download_rfam_cm(rfam_id: str) -> tuple[str, str | None]:
    """
    Download a CM model from RFAM.

    Args:
        rfam_id: RFAM family ID (e.g., RF00001)

    Returns:
        Tuple of (rfam_id, cm_content) or (rfam_id, None) on error
    """
    url = f"https://rfam.org/family/{rfam_id}/cm"
    try:
        response = requests.get(url, timeout=30)
        response.raise_for_status()
        return rfam_id, response.text
    except requests.RequestException as e:
        warn(f"Failed to download CM for {rfam_id}: {e}")
        return rfam_id, None


def download_rfam_cms(rfam_ids: list[str], num_threads: int = 4) -> str:
    """
    Download multiple RFAM CM models in parallel and concatenate them.

    Args:
        rfam_ids: List of RFAM family IDs
        num_threads: Number of parallel downloads

    Returns:
        Concatenated CM content as string
    """
    info(f"📥\tDownloading {len(rfam_ids)} CM models from RFAM...")

    cms = []
    failed = []

    with ThreadPoolExecutor(max_workers=min(num_threads, len(rfam_ids))) as executor:
        futures = {executor.submit(download_rfam_cm, rfid): rfid for rfid in rfam_ids}

        for future in as_completed(futures):
            rfam_id, cm_content = future.result()
            if cm_content:
                cms.append(cm_content)
                info(f"   ✓ Downloaded {rfam_id}")
            else:
                failed.append(rfam_id)

    if failed:
        warn(f"Failed to download: {', '.join(failed)}")

    if not cms:
        error("No CM models could be downloaded")
        raise RuntimeError("Failed to download any CM models from RFAM")

    info(f"   Downloaded {len(cms)}/{len(rfam_ids)} CM models")
    return "\n".join(cms)


def run_ncrna(all_neigh, den_data, output, num_threads, valid_unique_ids, ncrna_input: str):
    """
    Run Infernal for ncRNA annotation.

    Args:
        ncrna_input: Either a path to a CM file or comma-separated RFAM IDs (e.g., RF00001,RF00002)
    """
    info("🔬\tRunning Infernal for ncRNA annotation...")
    output = Path(output)
    ncrna_dir = output / "ncrna"
    ncrna_dir.mkdir(parents=True, exist_ok=True)

    # Parse input to determine if it's a path or RFAM IDs
    is_rfam, parsed_value = parse_ncrna_input(ncrna_input)

    temp_cm_file = None
    if is_rfam:
        # Download CMs from RFAM and create temporary concatenated file
        cm_content = download_rfam_cms(parsed_value, num_threads)
        temp_cm_file = tempfile.NamedTemporaryFile(  # noqa: SIM115
            mode="w", suffix=".cm", delete=False, dir=ncrna_dir
        )
        temp_cm_file.write(cm_content)
        temp_cm_file.close()
        cm_path = Path(temp_cm_file.name)
        info(f"   Created temporary CM file: {cm_path.name}")
    else:
        cm_path = parsed_value

    stockholm_file = ncrna_dir / "results.sto"
    tblout_file = ncrna_dir / "results.txt"

    try:
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
            str(cm_path),
            str(output / "neighborhood" / "neighborhoods.fasta"),
        ]
        subprocess.run(command, check=True, stdout=subprocess.DEVNULL)
    finally:
        # Clean up temporary file if created
        if temp_cm_file is not None:
            try:
                Path(temp_cm_file.name).unlink()
                info("   Cleaned up temporary CM file")
            except OSError:
                pass

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
