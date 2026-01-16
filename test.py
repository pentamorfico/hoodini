"""
Remote BLAST seeding via NCBI web UI (Firefox + Selenium) + download via requests.

Requires:
  - firefox
  - geckodriver
  - selenium
  - requests
"""

from __future__ import annotations

import re
import subprocess
import time
from pathlib import Path

import requests

from hoodini.utils.logging_utils import error, info, warn

UNIPROT_RE = re.compile(r"^[A-NR-Z][0-9][A-Z0-9]{3}[0-9](?:-[0-9]+)?$")
VALID_MAX_SEQS = [10, 50, 100, 250, 500, 1000, 5000]


def _pick_dropdown_value(max_targets: int) -> int:
    for opt in VALID_MAX_SEQS:
        if max_targets <= opt:
            return opt
    return VALID_MAX_SEQS[-1]


def _looks_like_fasta(text: str) -> bool:
    return text.strip().startswith(">") or ("\n" in text and len(text.strip()) > 0)


def _fetch_fasta_for_id(prot_id: str) -> str:
    prot_id = prot_id.strip()

    # UniProt first
    if UNIPROT_RE.match(prot_id):
        url = f"https://rest.uniprot.org/uniprotkb/{prot_id}.fasta"
        try:
            resp = requests.get(url, timeout=15)
            if resp.status_code == 200 and resp.text.strip().startswith(">"):
                return resp.text
            warn(f"UniProt fetch failed ({resp.status_code}) for {prot_id}")
        except Exception as e:
            warn(f"UniProt fetch error for {prot_id}: {e}")

    # Fallback to NCBI efetch (entrez-direct)
    cmd = ["efetch", "-db", "protein", "-id", prot_id, "-format", "fasta"]
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=30)
        if result.stdout and result.stdout.strip().startswith(">"):
            return result.stdout
        error(f"efetch returned no FASTA for {prot_id}")
        return ""
    except Exception as e:
        error(f"efetch failed for {prot_id}: {e}")
        return ""


def _run_remote_blast_firefox_selenium(
    fasta_text: str,
    evalue: float,
    max_targets: int,
) -> list[str]:
    try:
        from selenium import webdriver
        from selenium.webdriver.common.by import By
        from selenium.webdriver.support.ui import WebDriverWait, Select
        from selenium.webdriver.support import expected_conditions as EC
    except Exception:
        error("Selenium not available. Install with: mamba install -c conda-forge selenium")
        return []

    dropdown_value = _pick_dropdown_value(max_targets)

    info("🚀 BLAST Search (Firefox + Selenium)")
    info(f"   Max sequences: {max_targets}")
    if dropdown_value != max_targets:
        info(f"   (using dropdown: {dropdown_value}, will limit download to {max_targets})")
    info(f"   E-value: {evalue}")
    info("")

    options = webdriver.FirefoxOptions()
    options.add_argument("-headless")

    driver = None
    try:
        driver = webdriver.Firefox(options=options)
        wait = WebDriverWait(driver, 60)

        # 1) Open BLAST proteins page
        driver.get("https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins")

        # 2) Fill textarea
        textarea = None
        for by, sel in [
            (By.CSS_SELECTOR, 'textarea[aria-label*="accession" i]'),
            (By.CSS_SELECTOR, "textarea"),
        ]:
            try:
                textarea = wait.until(EC.presence_of_element_located((by, sel)))
                if textarea:
                    break
            except Exception:
                pass

        if not textarea:
            error("❌ Could not find BLAST query textarea.")
            return []

        textarea.clear()
        textarea.send_keys(fasta_text)

        # 3) Expand Algorithm parameters (robust)
        try:
            algo = wait.until(
                EC.element_to_be_clickable(
                    (
                        By.XPATH,
                        "//summary[contains(., 'Algorithm parameters')]"
                        " | //a[contains(., 'Algorithm parameters')]"
                        " | //button[contains(., 'Algorithm parameters')]",
                    )
                )
            )
            driver.execute_script("arguments[0].click();", algo)
            time.sleep(0.3)
        except Exception:
            pass

        # 4) Set MAX_NUM_SEQ dropdown
        try:
            max_sel = wait.until(EC.presence_of_element_located((By.NAME, "MAX_NUM_SEQ")))
            driver.execute_script("arguments[0].scrollIntoView({block:'center'});", max_sel)
            time.sleep(0.2)
            Select(max_sel).select_by_value(str(dropdown_value))
        except Exception as e:
            warn(f"Could not set MAX_NUM_SEQ: {e}")

        # 5) Set e-value (robust: scroll + JS fallback)
        try:
            ev = wait.until(EC.presence_of_element_located((By.NAME, "EXPECT")))
            driver.execute_script("arguments[0].scrollIntoView({block:'center'});", ev)
            time.sleep(0.2)

            try:
                ev.clear()
                ev.send_keys(str(evalue))
            except Exception:
                driver.execute_script("arguments[0].value = arguments[1];", ev, str(evalue))
                driver.execute_script(
                    "arguments[0].dispatchEvent(new Event('change', {bubbles:true}));", ev
                )
        except Exception as e:
            warn(f"Could not set EXPECT: {e}")

        # 6) Click BLAST
        blast_btn = None
        for by, sel in [
            (By.CSS_SELECTOR, '#blastButton1 input[value="BLAST"]'),
            (By.CSS_SELECTOR, 'input[value="BLAST"]'),
        ]:
            try:
                blast_btn = wait.until(EC.element_to_be_clickable((by, sel)))
                if blast_btn:
                    break
            except Exception:
                pass

        if not blast_btn:
            error("❌ Could not find BLAST button.")
            return []

        driver.execute_script("arguments[0].scrollIntoView({block:'center'});", blast_btn)
        time.sleep(0.2)
        driver.execute_script("arguments[0].click();", blast_btn)

        # 7) Extract RID from URL or HTML
        rid = None
        t0 = time.time()
        while time.time() - t0 < 60:
            cur_url = driver.current_url
            m = re.search(r"[&?]RID=([A-Z0-9]+)", cur_url)
            if m:
                rid = m.group(1)
                break

            html = driver.page_source
            m = re.search(r"Request ID[^A-Z0-9]*([A-Z0-9]{11,12})", html)
            if m:
                rid = m.group(1)
                break

            time.sleep(0.5)

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

            time.sleep(1)
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
        if driver is not None:
            try:
                driver.quit()
            except Exception:
                pass


def prepare_single_query_input(
    query: str,
    output_dir: Path,
    evalue: float = 1e-5,
    max_targets: int = 100,
) -> Path | None:
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

    hits = _run_remote_blast_firefox_selenium(
        fasta_txt,
        evalue=evalue,
        max_targets=max_targets,
    )
    if not hits:
        error("Remote BLAST returned no hits; aborting.")
        return None

    unique_hits: list[str] = []
    seen: set[str] = set()
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


if __name__ == "__main__":
    out = prepare_single_query_input(
        query="WP_455324848.1",
        output_dir=Path("blast_seed_test"),
        evalue=1e-5,
        max_targets=100,
    )
    print(out)
