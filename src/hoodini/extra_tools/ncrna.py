import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import polars as pl
import pyhmmer
import pyinfernal
import requests

from hoodini.utils.logging_utils import error, info, warn

# Regex pattern for RFAM IDs: RF followed by exactly 5 digits
RFAM_PATTERN = re.compile(r"^RF\d{5}$")

# --- Consensus secondary-structure derivation (.cm parsing) -------------------
#
# pyinfernal does not expose consensus secondary structure (WUSS/dot-bracket)
# anywhere in its Python API, so it is reimplemented here directly from the
# plain-text Infernal ``.cm`` format. This reimplements the node/state
# traversal Infernal itself uses internally (``CreateEmitMap()`` in
# ``display.c``) to compute, for each consensus column of a CM, whether it is
# unpaired or base-paired with another column, with no dependency on the
# ``cmsearch``/``cmemit`` binaries.
#
# Node numbers (``nd``) are assigned by Infernal at model-construction time
# and are NOT necessarily encountered in file/print order: only the *main
# chain* (non-BIF nodes) is guaranteed contiguous, where a node's child is
# always ``nd + 1``. Node headers are therefore parsed by their explicit
# index, not by encounter order, and bifurcation (``BIF``) children are
# resolved via the ``cfirst``/``cnum`` state pointers on their bifurcation
# state, exactly as Infernal does.
#
# Validated against ``cmsearch``/``cmalign``-derived ``SS_cons`` ground truth
# for five structurally diverse Rfam models (5S rRNA, tRNA, tmRNA, and the
# much larger, heavily-bifurcated SSU/LSU rRNA models), with byte-for-byte
# matches in all cases.

_CM_NODE_RE = re.compile(r"\[\s*(\w+)\s+(\d+)\s*\]")
_CM_STATE_TYPES = {"S", "B", "D", "MP", "ML", "MR", "IL", "IR", "E"}


def parse_cm_consensus_structures(cm_path: str) -> dict[str, str]:
    """Parse a (possibly multi-model) ``.cm`` file into consensus dot-bracket structures.

    Args:
        cm_path: Path to a plain-text Infernal ``.cm`` file, containing one or
            more concatenated covariance models.

    Returns:
        Mapping of CM ``NAME`` to its consensus structure string, using a
        simplified Vienna-style dot-bracket notation (``(``, ``)``, ``.``) of
        length ``CLEN``, indexed the same way as ``cm_from``/``cm_to`` reported
        by ``pyinfernal`` hits (1-based, inclusive).
    """
    structures = {}
    with open(cm_path) as f:
        lines = f.readlines()

    # Split into per-model blocks (each ends with a line containing only "//").
    blocks = []
    start = 0
    for i, line in enumerate(lines):
        if line.strip() == "//":
            blocks.append(lines[start : i + 1])
            start = i + 1

    for block in blocks:
        # Rfam distributes each CM bundled with a companion HMMER3 filter
        # profile (used internally by cmsearch), concatenated in the same
        # file and also delimited by "//". Multiple models concatenated
        # together may also leave blank separator lines between blocks.
        # Skip anything that isn't a CM block.
        first_line = next((line for line in block if line.strip()), "")
        if not first_line.startswith("INFERNAL"):
            continue
        name, structure = _parse_single_cm(block)
        if name is not None:
            structures[name] = structure
    return structures


def _parse_single_cm(lines: list[str]) -> tuple[str | None, str]:
    name = None
    clen = None
    nnodes = None
    ndtype: list[str | None] = []
    nodemap: list[int | None] = []
    cfirst: dict[int, int] = {}
    cnum: dict[int, int] = {}
    ndidx_of_state: dict[int, int] = {}

    current_nd = None
    for line in lines:
        if line.startswith("NAME"):
            name = line.split(None, 1)[1].strip()
            continue
        if line.startswith("CLEN"):
            clen = int(line.split()[1])
            continue
        if line.startswith("NODES"):
            nnodes = int(line.split()[1])
            ndtype = [None] * nnodes
            nodemap = [None] * nnodes
            continue

        m = _CM_NODE_RE.search(line)
        if m and line.lstrip().startswith("["):
            ntype, nidx = m.groups()
            nidx = int(nidx)
            ndtype[nidx] = ntype
            current_nd = nidx
            continue

        parts = line.split()
        if not parts or current_nd is None or parts[0] not in _CM_STATE_TYPES:
            continue
        v = int(parts[1])
        if nodemap[current_nd] is None:
            nodemap[current_nd] = v
        ndidx_of_state[v] = current_nd
        cfirst[v] = int(parts[4])
        cnum[v] = int(parts[5])

    if clen is None or nnodes is None:
        return name, ""

    # Reimplementation of Infernal's CreateEmitMap() (display.c): an iterative
    # pre-order traversal assigning each MATP/MATL/MATR node its consensus
    # column position(s), recursing into BIF children via their state pointers.
    lpos: list[int | None] = [None] * nnodes
    rpos: list[int | None] = [None] * nnodes
    cpos = 0
    stack = [(0, 0)]
    while stack:
        nd, on_right = stack.pop()
        if on_right:
            rpos[nd] = cpos + 1
            if ndtype[nd] in ("MATP", "MATR"):
                cpos += 1
        else:
            if ndtype[nd] in ("MATP", "MATL"):
                cpos += 1
            lpos[nd] = cpos
            if ndtype[nd] == "BIF":
                v0 = nodemap[nd]
                right_child_nd = ndidx_of_state[cnum[v0]]
                left_child_nd = ndidx_of_state[cfirst[v0]]
                stack.append((nd, 1))
                stack.append((right_child_nd, 0))
                stack.append((left_child_nd, 0))
            else:
                stack.append((nd, 1))
                if ndtype[nd] != "END":
                    stack.append((nd + 1, 0))

    structure = ["."] * (clen + 2)  # 1-indexed, +1 slack like Infernal (0..clen+1)
    for nd in range(nnodes):
        if ndtype[nd] == "MATP":
            structure[lpos[nd]] = "("
            structure[rpos[nd]] = ")"
        elif ndtype[nd] == "MATL":
            structure[lpos[nd]] = "."
        elif ndtype[nd] == "MATR":
            structure[rpos[nd]] = "."

    return name, "".join(structure[1 : clen + 1])


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
    Run Infernal (via pyinfernal) for ncRNA annotation.

    Args:
        ncrna_input: Either a path to a CM file or comma-separated RFAM IDs (e.g., RF00001,RF00002)
    """
    info("🔬\tRunning Infernal (pyinfernal) for ncRNA annotation...")
    output = Path(output)
    ncrna_dir = output / "ncrna"
    ncrna_dir.mkdir(parents=True, exist_ok=True)

    # Parse input to determine if it's a path or RFAM IDs
    is_rfam, parsed_value = parse_ncrna_input(ncrna_input)

    if is_rfam:
        # Download CMs from RFAM and write a concatenated CM file
        cm_content = download_rfam_cms(parsed_value, num_threads)
        cm_path = ncrna_dir / "downloaded_models.cm"
        cm_path.write_text(cm_content)
        info(f"   Wrote downloaded CM file: {cm_path.name}")
    else:
        cm_path = parsed_value

    fasta_path = output / "neighborhood" / "neighborhoods.fasta"

    # Consensus secondary structure (per CM model), derived purely from the
    # .cm text file without needing the `cmsearch`/`cmemit` binaries, since
    # pyinfernal does not expose it directly.
    consensus_structures = parse_cm_consensus_structures(str(cm_path))

    alphabet = pyhmmer.easel.Alphabet.rna()
    with pyhmmer.easel.SequenceFile(
        fasta_path, format="fasta", digital=True, alphabet=alphabet
    ) as seq_file:
        sequences = list(seq_file)

    if not sequences:
        warn(f"No sequences found in {fasta_path}")
        empty_df = pl.DataFrame()
        empty_df.write_csv(ncrna_dir / "ncrna_results.tsv", separator="\t", include_header=False)
        return empty_df

    # Pre-fetch targets into a DigitalSequenceBlock: as of pyinfernal 0.1.x,
    # `cmsearch()` only auto-computes the pipeline's Z parameter (total
    # database length, needed for E-value calibration) when passed a
    # DigitalSequenceBlock directly; building it ourselves works around that.
    targets = pyhmmer.easel.DigitalSequenceBlock(alphabet, sequences)

    with pyinfernal.cm.CMFile(str(cm_path), alphabet=alphabet) as cm_file:
        cms = list(cm_file)

    if not cms:
        error(f"No covariance models could be read from {cm_path}")
        raise RuntimeError(f"Failed to load any CM models from {cm_path}")

    rows = []
    for hits in pyinfernal.cmsearch(cms, targets, cpus=num_threads, E=1e-5, incE=1e-5):
        query_name = hits.query.name
        query_accession = hits.query.accession
        structure = consensus_structures.get(query_name, "")
        for hit in hits:
            aln = hit.alignment
            sequence = aln.target_sequence.replace("-", "")
            hit_structure = _slice_hit_structure(
                structure, aln.cm_from, aln.cm_sequence, aln.target_sequence
            )
            rows.append(
                {
                    "nucid": hit.name,
                    "--": query_accession,
                    "nc_feature": query_name,
                    "cm": "cm",
                    "mdlfrom": aln.cm_from,
                    "mdlto": aln.cm_to,
                    "seqfrom": aln.target_from,
                    "seqto": aln.target_to,
                    "strand_ncrna": hit.strand,
                    "score": hit.score,
                    "E-value": hit.evalue,
                    "sequence": sequence,
                    "structure": hit_structure,
                }
            )

    if not rows:
        warn("No ncRNA found by Infernal (pyinfernal)")
        empty_df = pl.DataFrame()
        empty_df.write_csv(ncrna_dir / "ncrna_results.tsv", separator="\t", include_header=False)
        return empty_df

    cmdf = pl.DataFrame(rows)
    cmdf = cmdf.with_columns(
        [
            pl.col("seqfrom").cast(pl.Int64),
            pl.col("seqto").cast(pl.Int64),
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


def _slice_hit_structure(
    consensus_structure: str, cm_from: int, cm_sequence: str, target_sequence: str
) -> str:
    """Slice a per-model consensus structure down to a single hit's alignment.

    Walks the CM/target alignment columns in lockstep with the model's
    consensus-column pointer (starting at ``cm_from``): match columns
    (uppercase in ``cm_sequence``) consume one structure character and
    advance the pointer; insert columns (lowercase) are always unpaired in
    the output; columns deleted in the target (``-``) consume the pointer
    but contribute nothing to the output, keeping the result aligned 1:1
    with the gap-free hit sequence.
    """
    pointer = cm_from
    out = []
    for cm_char, target_char in zip(cm_sequence, target_sequence):
        is_match = cm_char.isupper()
        if target_char == "-":
            if is_match:
                pointer += 1
            continue
        if is_match:
            out.append(consensus_structure[pointer - 1])
            pointer += 1
        else:
            out.append(".")
    return "".join(out)
