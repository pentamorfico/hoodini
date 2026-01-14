"""Utilities to seed the pipeline from a single protein ID or FASTA string."""

from __future__ import annotations

import re
import subprocess
from pathlib import Path

import requests

try:
    from Bio import SearchIO
    from Bio.Blast import NCBIWWW
except Exception:  # biopython optional; fallback handled below
    NCBIWWW = None
    SearchIO = None

from hoodini.utils.logging_utils import error, info, warn

UNIPROT_RE = re.compile(r"^[A-NR-Z][0-9][A-Z0-9]{3}[0-9](?:-[0-9]+)?$")


def _looks_like_fasta(text: str) -> bool:
    return text.strip().startswith(">") or ("\n" in text and len(text.strip()) > 0)


def _fetch_fasta_for_id(prot_id: str) -> str:
    """Fetch protein FASTA from NCBI (efetch) or UniProt."""
    prot_id = prot_id.strip()
    if UNIPROT_RE.match(prot_id):
        url = f"https://rest.uniprot.org/uniprotkb/{prot_id}.fasta"
        try:
            resp = requests.get(url, timeout=15)
            if resp.status_code == 200 and resp.text.strip().startswith(">"):
                return resp.text
            warn(f"UniProt fetch failed ({resp.status_code}) for {prot_id}")
        except Exception as e:
            warn(f"UniProt fetch error for {prot_id}: {e}")

    # Fallback to NCBI efetch
    cmd = [
        "efetch",
        "-db",
        "protein",
        "-id",
        prot_id,
        "-format",
        "fasta",
    ]
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=30)
        if result.stdout and result.stdout.strip().startswith(">"):
            return result.stdout
        error(f"efetch returned no FASTA for {prot_id}")
        return ""
    except Exception as e:
        error(f"efetch failed for {prot_id}: {e}")
        return ""


def _run_remote_blast(
    fasta_text: str,
    evalue: float,
    max_targets: int,
    db: str = "nr_cluster_seq",
) -> list[str]:
    def _clean_blast_id(val: str) -> str:
        val = val.strip()
        if "|" in val:
            tokens = [t for t in val.split("|") if t]
            if len(tokens) >= 2:
                return tokens[1]
            if tokens:
                return tokens[0]
        return val

    if NCBIWWW is None or SearchIO is None:
        error("Biopython with SearchIO is required for remote BLAST; please install biopython.")
        return []

    def _qblast_once(database: str) -> list[str]:
        info(f"🔍  Running remote BLAST (qblast) on {database}...")
        handle = NCBIWWW.qblast(
            program="blastp",
            database=database,
            sequence=fasta_text,
            expect=evalue,
            hitlist_size=max_targets,
            format_type="XML",
        )
        try:
            record = SearchIO.read(handle, "blast-xml")
        except Exception as e_parse:  # noqa: BLE001
            warn(f"Failed to parse BLAST XML: {e_parse}")
            return []
        return [_clean_blast_id(hit.id) for hit in record.hits]

    try:
        hits = _qblast_once(db)
        if hits:
            return hits
        if db != "nr":
            warn(f"No hits from {db}; retrying on nr")
            hits = _qblast_once("nr")
        return hits
    except Exception as e_qblast:  # noqa: BLE001
        error(f"Remote BLAST (qblast) failed on {db}: {e_qblast}")
        return []


def prepare_single_query_input(
    query: str,
    output_dir: Path,
    evalue: float = 1e-5,
    max_targets: int = 100,
    db: str = "nr_cluster_seq",
) -> Path | None:
    """
    Given a query (protein ID or FASTA string), run remote BLAST and
    emit a single-column input list file with the hit IDs.
    """
    query = query.strip()
    output_dir.mkdir(parents=True, exist_ok=True)
    if _looks_like_fasta(query) or re.fullmatch(r"[A-Z*]+", query.replace("\n", ""), re.I):
        info("⚙️  Using provided FASTA/sequence as query.")
        fasta_txt = query if query.startswith(">") else f">query\n{query}\n"
    else:
        info(f"⚙️  Fetching FASTA for query ID: {query}")
        fasta_txt = _fetch_fasta_for_id(query)
        if not fasta_txt.strip():
            error(f"Could not fetch FASTA for query '{query}'.")
            return None

    hits = _run_remote_blast(fasta_txt, evalue=evalue, max_targets=max_targets, db=db)
    if not hits:
        error("Remote BLAST returned no hits; aborting.")
        return None

    unique_hits = []
    seen = set()
    for h in hits:
        if h not in seen:
            unique_hits.append(h)
            seen.add(h)

    # Include the query ID itself if not already present
    if not _looks_like_fasta(query) and query not in seen:
        unique_hits.insert(0, query)

    input_list_path = output_dir / "input_from_blast.txt"
    input_list_path.write_text("\n".join(unique_hits), encoding="utf-8")
    info(f"✔️  Seeded {len(unique_hits)} protein IDs from remote BLAST.")
    return input_list_path
