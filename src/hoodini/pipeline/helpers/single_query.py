"""Utilities to seed the pipeline from a single protein ID or FASTA string."""

from __future__ import annotations

import asyncio
import re
import subprocess
from pathlib import Path

import requests
from playwright.async_api import async_playwright

from hoodini.utils.browser_setup import ensure_playwright_firefox
from hoodini.utils.logging_utils import error, info, warn
from hoodini.utils.runtime_env import apply_ld_library_path

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




async def _run_remote_blast(
    fasta_text: str,
    evalue: float,
    max_targets: int,
    db: str = "nr_cluster_seq",
) -> list[str]:
    """Run remote BLAST via NCBI using Playwright Firefox."""
    dropdown_value = _pick_dropdown_value(max_targets)

    # Ensure Firefox is installed
    if not ensure_playwright_firefox():
        error("❌ Could not install Playwright Firefox")
        return []

    info("🚀 BLAST Search (Firefox + Playwright)")
    info(f"   Max sequences: {max_targets}")
    if dropdown_value != max_targets:
        info(f"   (using dropdown: {dropdown_value}, will limit download to {max_targets})")
    info(f"   E-value: {evalue}")
    info("")

    # Set up environment for Playwright Firefox in conda/pixi/mamba
    apply_ld_library_path()

    async with async_playwright() as p:
        browser = await p.firefox.launch(headless=True)
        page = await browser.new_page()

        try:
            # 1) Open BLAST proteins page
            await page.goto("https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins")

            # 2) Fill textarea with FASTA
            await page.fill('textarea[aria-label*="accession" i]', fasta_text)

            # 3) Try to expand Algorithm parameters
            try:
                await page.click("summary:has-text('Algorithm parameters')")
                await page.wait_for_timeout(300)
            except Exception:
                pass

            # 4) Set MAX_NUM_SEQ dropdown
            try:
                await page.select_option('select[name="MAX_NUM_SEQ"]', value=str(dropdown_value))
            except Exception as e:
                warn(f"Could not set MAX_NUM_SEQ: {e}")

            # 5) Set e-value
            try:
                await page.fill('input[name="EXPECT"]', str(evalue))
            except Exception as e:
                warn(f"Could not set EXPECT: {e}")

            # 6) Click BLAST button
            await page.click('input[value="BLAST"]')

            # 7) Extract RID from URL or HTML
            rid = None
            for _ in range(120):  # 60 seconds
                cur_url = page.url
                m = re.search(r"[&?]RID=([A-Z0-9]+)", cur_url)
                if m:
                    rid = m.group(1)
                    break

                content = await page.content()
                m = re.search(r"Request ID[^A-Z0-9]*([A-Z0-9]{11,12})", content)
                if m:
                    rid = m.group(1)
                    break

                await page.wait_for_timeout(500)

            if not rid:
                error("❌ Could not find RID.")
                return []

            rid_clean = rid[4:] if rid.startswith("RID-") else rid

            # 8) Poll status via requests
            status_url = (
                "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
                f"?CMD=Get&RID={rid_clean}&FORMAT_OBJECT=SearchInfo"
            )

            for _ in range(600):
                resp = requests.get(status_url, timeout=30)
                text = resp.text

                if "Status=READY" in text and "ThereAreHits=yes" in text:
                    break
                if "Status=READY" in text and "ThereAreHits=no" in text:
                    warn("⚠️ BLAST finished but ThereAreHits=no")
                    return []
                if "Status=FAILED" in text or "Status=UNKNOWN" in text:
                    error("❌ BLAST failed or unknown RID")
                    return []

                await page.wait_for_timeout(1000)
            else:
                error("⚠️ Timeout waiting for BLAST results.")
                return []

            # 9) Download CSV
            download_url = (
                "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
                f"?RESULTS_FILE=on&RID={rid_clean}"
                f"&FORMAT_TYPE=CSV"
                f"&DESCRIPTIONS={dropdown_value}"
                f"&ALIGNMENT_VIEW=Tabular"
                f"&CMD=Get"
            )

            resp = requests.get(download_url, timeout=60)
            content = resp.text.strip()
            if not content:
                error("❌ Empty BLAST CSV download.")
                return []

            lines = [l for l in content.split("\n") if l and not l.startswith("#")]
            lines = lines[:max_targets]

            hits: list[str] = []
            for line in lines:
                cols = line.split(",")
                if len(cols) >= 2:
                    hits.append(cols[1].strip().strip('"'))

            return hits

        finally:
            await browser.close()




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

    # Run async BLAST in synchronous context
    hits = asyncio.run(_run_remote_blast(fasta_txt, evalue=evalue, max_targets=max_targets, db=db))
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
