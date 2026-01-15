"""Utilities to seed the pipeline from a single protein ID or FASTA string."""

from __future__ import annotations

import re
import subprocess
import time
from pathlib import Path

import requests

try:
    from playwright.sync_api import TimeoutError as PlaywrightTimeoutError, sync_playwright
except Exception:  # playwright optional; fallback handled below
    PlaywrightTimeoutError = None
    sync_playwright = None

from hoodini.utils.logging_utils import error, info, warn

UNIPROT_RE = re.compile(r"^[A-NR-Z][0-9][A-Z0-9]{3}[0-9](?:-[0-9]+)?$")
VALID_MAX_SEQS = [10, 50, 100, 250, 500, 1000, 5000]


def _pick_dropdown_value(max_targets: int) -> int:
    """NCBI dropdown only allows specific target counts; pick the nearest above."""
    for opt in VALID_MAX_SEQS:
        if max_targets <= opt:
            return opt
    return VALID_MAX_SEQS[-1]


def _extract_rid(text: str) -> str | None:
    """Extract RID from URL or HTML text."""
    # Do NOT use this function; inline the regex instead
    match = re.search(r"[&?]RID=([A-Z0-9]+)", text)
    if match:
        return match.group(1)
    match = re.search(r"Request ID[^A-Z0-9]*([A-Z0-9]{11,12})", text)
    if match:
        return match.group(1)
    return None


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
    dropdown_value = _pick_dropdown_value(max_targets)

    if sync_playwright is None:
        error(
            "Playwright is required for remote BLAST; install with `pip install playwright` "
            "and run `playwright install chromium`."
        )
        return []

    info("🚀 BLAST Search")
    info(f"   Query: {fasta_text[:50]}...")
    info(f"   Max sequences: {max_targets}")
    if dropdown_value != max_targets:
        info(f"   (using dropdown: {dropdown_value}, will limit download to {max_targets})")
    info(f"   E-value: {evalue}")
    info("")

    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page(viewport={"width": 1920, "height": 1080})

        # 1. LOAD PAGE
        url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins"
        page.goto(url, wait_until="networkidle")

        # 2. FILL SEQUENCE
        textarea = page.locator('textarea[aria-label*="accession"]').or_(
            page.locator("textarea").first
        )
        textarea.wait_for(state="visible")
        textarea.fill(fasta_text)

        # 3. EXPAND "Algorithm parameters" section
        algo_params = page.locator("text=Algorithm parameters").first
        algo_params.click()
        time.sleep(0.5)

        # 4. SET MAX TARGET SEQUENCES (use dropdown_value)
        max_seqs_select = page.locator('select[name="MAX_NUM_SEQ"]')
        max_seqs_select.wait_for(state="visible")
        max_seqs_select.select_option(str(dropdown_value))

        # 5. SET E-VALUE
        evalue_input = page.locator('input[name="EXPECT"]')
        evalue_input.fill(str(evalue))

        # 6. CLICK BLAST BUTTON
        blast_btn = page.locator('#blastButton1 input[value="BLAST"]')
        blast_btn.wait_for(state="visible")

        with page.expect_navigation(timeout=60000, wait_until="commit"):
            blast_btn.click(no_wait_after=True)

        # 7. EXTRACT RID
        time.sleep(3)

        rid = None

        # Method 1: Get RID from URL parameter
        current_url = page.url
        match = re.search(r"[&?]RID=([A-Z0-9]+)", current_url)
        if match:
            rid = match.group(1)

        # Method 2: Fallback - look for "Request ID" in page
        if not rid:
            content = page.content()
            match = re.search(r"Request ID[^A-Z0-9]*([A-Z0-9]{11,12})", content)
            if match:
                rid = match.group(1)

        if not rid:
            error("❌ Could not find RID")
            return []

        # 8. POLL FOR RESULTS
        rid_clean = rid[4:] if rid.startswith("RID-") else rid

        status_url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID={rid_clean}&FORMAT_OBJECT=SearchInfo"

        for i in range(600):
            resp = requests.get(status_url)
            text = resp.text

            if "Status=READY" in text and "ThereAreHits=yes" in text:
                break
            elif "Status=FAILED" in text or "Status=UNKNOWN" in text:
                error("❌ BLAST failed or unknown RID")
                return []

            time.sleep(1)

        else:
            error("⚠️ Timeout waiting for BLAST results.")
            return []

        # 9. DOWNLOAD CSV (use dropdown_value to get all, then limit to max_seqs)
        download_url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?RESULTS_FILE=on&RID={rid_clean}&FORMAT_TYPE=CSV&DESCRIPTIONS={dropdown_value}&ALIGNMENT_VIEW=Tabular&CMD=Get"

        resp = requests.get(download_url)
        content = resp.text

        # Parse and limit to max_seqs
        all_lines = content.strip().split("\n")
        data_lines = [l for l in all_lines if l and not l.startswith("#")]

        # Limit data lines to requested max_seqs
        limited_lines = data_lines[:max_targets]

        hits = []
        for line in limited_lines:
            cols = line.split(",")
            if len(cols) >= 2:
                hits.append(cols[1].strip().strip('"'))

        return hits


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
