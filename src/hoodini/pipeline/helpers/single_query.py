"""Utilities to seed the pipeline from a single protein ID or FASTA string."""

from __future__ import annotations

import os
import re
import subprocess
import time
from pathlib import Path

import requests
from playwright.sync_api import sync_playwright

from hoodini.utils.browser_setup import ensure_playwright_firefox
from hoodini.utils.logging_utils import error, info, warn
from hoodini.utils.runtime_env import apply_ld_library_path

UNIPROT_RE = re.compile(r"^[A-NR-Z][0-9][A-Z0-9]{3}[0-9](?:-[0-9]+)?$")
VALID_MAX_SEQS_BLASTP = [10, 50, 100, 250, 500, 1000, 5000]
VALID_MAX_SEQS_PSIBLAST = [10, 50, 100, 250, 500, 1000, 5000, 10000, 20000]
PSI_BLAST_THRESHOLD = 5000


def _pick_dropdown_value(max_targets: int, use_psiblast: bool = False) -> int:
    """NCBI dropdown only allows specific target counts; pick the nearest above."""
    valid = VALID_MAX_SEQS_PSIBLAST if use_psiblast else VALID_MAX_SEQS_BLASTP
    for opt in valid:
        if max_targets <= opt:
            return opt
    return valid[-1]


def _extract_rid(text: str) -> str | None:
    """Extract RID from URL or HTML text."""
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
    """Run remote BLAST via NCBI using Playwright Firefox.

    Automatically switches to PSI-BLAST when ``max_targets`` exceeds 5000,
    since NCBI blastp caps at 5000 while PSI-BLAST supports up to 20000.
    """
    use_psiblast = max_targets > PSI_BLAST_THRESHOLD
    dropdown_value = _pick_dropdown_value(max_targets, use_psiblast=use_psiblast)
    program_label = "PSI-BLAST" if use_psiblast else "BLASTp"

    if not ensure_playwright_firefox():
        error("❌ Could not install Playwright Firefox")
        return []

    info(f"🚀 {program_label} Search")
    info(f"   Query: {fasta_text[:50]}...")
    info(f"   Max sequences: {max_targets}")
    if dropdown_value != max_targets:
        info(f"   (using dropdown: {dropdown_value}, will limit download to {max_targets})")
    info(f"   E-value: {evalue}")
    if use_psiblast:
        info(f"   Using PSI-BLAST (max_targets > {PSI_BLAST_THRESHOLD})")
    info("")

    original_ld_path = os.environ.get("LD_LIBRARY_PATH")
    apply_ld_library_path()

    try:
        return _playwright_blast(fasta_text, evalue, max_targets, dropdown_value, use_psiblast)
    finally:
        if original_ld_path is None:
            os.environ.pop("LD_LIBRARY_PATH", None)
        else:
            os.environ["LD_LIBRARY_PATH"] = original_ld_path


def _playwright_blast(
    fasta_text: str,
    evalue: float,
    max_targets: int,
    dropdown_value: int,
    use_psiblast: bool,
) -> list[str]:
    """Run the actual Playwright browser session for BLAST."""
    with sync_playwright() as p:
        browser = p.firefox.launch(headless=True)
        page = browser.new_page(viewport={"width": 1920, "height": 1080})

        url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins"
        page.goto(url, wait_until="networkidle")

        if use_psiblast:
            page.locator("text=PSI-BLAST (Position-Specific Iterated BLAST)").click()
            time.sleep(0.5)

        textarea = page.locator('textarea[aria-label*="accession"]').or_(
            page.locator("textarea").first
        )
        textarea.wait_for(state="visible")
        textarea.fill(fasta_text)

        algo_params = page.locator("text=Algorithm parameters").first
        algo_params.click()
        time.sleep(0.5)

        max_seqs_select = page.locator('select[name="MAX_NUM_SEQ"]')
        max_seqs_select.wait_for(state="visible")
        max_seqs_select.select_option(str(dropdown_value))

        evalue_input = page.locator('input[name="EXPECT"]')
        evalue_input.fill(str(evalue))

        blast_btn = page.locator('#blastButton1 input[value="BLAST"]')
        blast_btn.wait_for(state="visible")

        with page.expect_navigation(timeout=60000, wait_until="commit"):
            blast_btn.click(no_wait_after=True)

        rid = None
        for attempt in range(10):
            time.sleep(3)

            current_url = page.url
            match = re.search(r"[&?]RID=([A-Z0-9]+)", current_url)
            if match:
                rid = match.group(1)
                break

            content = page.content()
            match = re.search(r"Request ID[^A-Z0-9]*([A-Z0-9]{11,12})", content)
            if match:
                rid = match.group(1)
                break

        if not rid:
            error("❌ Could not find RID")
            info(f"debug: Current URL: {page.url}")
            info(f"debug: Page Title: {page.title()}")
            content_snippet = page.content()[:1000].replace("\n", " ")
            info(f"debug: Page Content Start: {content_snippet}...")
            return []

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

        download_url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?RESULTS_FILE=on&RID={rid_clean}&FORMAT_TYPE=CSV&DESCRIPTIONS={dropdown_value}&ALIGNMENT_VIEW=Tabular&CMD=Get"

        resp = requests.get(download_url)
        content = resp.text

        all_lines = content.strip().split("\n")
        data_lines = [l for l in all_lines if l and not l.startswith("#")]
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

    if not _looks_like_fasta(query) and query not in seen:
        unique_hits.insert(0, query)

    input_list_path = output_dir / "input_from_blast.txt"
    input_list_path.write_text("\n".join(unique_hits), encoding="utf-8")
    info(f"✔️  Seeded {len(unique_hits)} protein IDs from remote BLAST.")
    return input_list_path
