"""Utilities to seed the pipeline from a single protein ID or FASTA string."""

from __future__ import annotations

import contextlib
import json
import re
import time
from pathlib import Path

import requests

from hoodini.utils.browser_setup import ensure_lightpanda
from hoodini.utils.logging_utils import error, info, warn
from hoodini.utils.mcp_browser import MCPBrowserError, MCPBrowserSession, find_lightpanda_binary
from hoodini.utils.ncbi_eutils import efetch

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


def _looks_like_fasta(text: str) -> bool:
    return text.strip().startswith(">") or ("\n" in text and len(text.strip()) > 0)


def _fetch_fasta_for_id(prot_id: str) -> str:
    """Fetch protein FASTA from NCBI (efetch.fcgi) or UniProt."""
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

    try:
        text = efetch(db="protein", ids=prot_id, rettype="fasta", retmode="text", timeout=30)
        if text and text.strip().startswith(">"):
            return text
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
    """Run remote BLAST via NCBI using lightpanda (headless browser over MCP stdio, no port).

    Automatically switches to PSI-BLAST when ``max_targets`` exceeds 5000,
    since NCBI blastp caps at 5000 while PSI-BLAST supports up to 20000.
    """
    use_psiblast = max_targets > PSI_BLAST_THRESHOLD
    dropdown_value = _pick_dropdown_value(max_targets, use_psiblast=use_psiblast)
    program_label = "PSI-BLAST" if use_psiblast else "BLASTp"

    if not ensure_lightpanda():
        error("❌ Could not start lightpanda")
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

    hits = _browser_blast(fasta_text, evalue, max_targets, dropdown_value, use_psiblast)
    if not hits and find_lightpanda_binary():
        warn("browser driver returned no hits — falling back to lightpanda")
        hits = _browser_blast(
            fasta_text, evalue, max_targets, dropdown_value, use_psiblast, driver="lightpanda"
        )
    return hits


def _build_blast_submit_js(
    fasta_text: str, evalue: float, dropdown_value: int, use_psiblast: bool
) -> str:
    """Build the in-page JS that fills the BLAST form fields and submits it.

    Only the fields hoodini cares about are touched (QUERY, MAX_NUM_SEQ,
    EXPECT and, when needed, the PSI-BLAST radio); everything else keeps the
    page's own defaults. Submission uses ``form.submit()`` directly rather
    than clicking the BLAST button, since that button has ``type="button"``
    and is wired up via jQuery delegate handlers that don't reliably trigger
    from synthetic events. The result is JSON-stringified so it survives the
    round-trip through the browser tool as text.
    """
    seq = re.sub(r"^>.*\n?", "", fasta_text).replace("\n", "").strip()
    psi_js = (
        """
      var psiRadio = document.querySelector('input[name="BLAST_PROGRAMS"][value="psiBlast"]');
      if (psiRadio) {
        psiRadio.checked = true;
        psiRadio.dispatchEvent(new Event('change', {bubbles: true}));
      }
        """
        if use_psiblast
        else ""
    )
    return f"""
    (() => {{
      var q = document.querySelector('textarea[name="QUERY"]') || document.querySelector('textarea');
      if (!q) return JSON.stringify({{ok: false, reason: 'query textarea not found'}});
      q.focus();
      q.value = {json.dumps(seq)};
      q.dispatchEvent(new Event('input', {{bubbles: true}}));
      q.dispatchEvent(new Event('change', {{bubbles: true}}));
      {psi_js}
      var maxSeqs = document.querySelector('select[name="MAX_NUM_SEQ"]');
      if (maxSeqs) {{
        var match = Array.from(maxSeqs.options).find(o => String(o.value) === String({dropdown_value}));
        if (!match) {{
          match = document.createElement('option');
          match.value = String({dropdown_value});
          match.text = String({dropdown_value});
          maxSeqs.add(match);
        }}
        maxSeqs.value = match.value;
        maxSeqs.dispatchEvent(new Event('change', {{bubbles: true}}));
      }}
      var expect = document.querySelector('input[name="EXPECT"]');
      if (expect) {{
        expect.value = {json.dumps(str(evalue))};
        expect.dispatchEvent(new Event('input', {{bubbles: true}}));
        expect.dispatchEvent(new Event('change', {{bubbles: true}}));
      }}
      var form = q.form;
      if (!form) return JSON.stringify({{ok: false, reason: 'form not found'}});
      form.submit();
      return JSON.stringify({{ok: true}});
    }})()
    """


def _browser_blast(
    fasta_text: str,
    evalue: float,
    max_targets: int,
    dropdown_value: int,
    use_psiblast: bool,
    driver: str | None = None,
) -> list[str]:
    """Run the browser session for BLAST (obscura or lightpanda over MCP stdio, no port)."""
    mcp = MCPBrowserSession(driver=driver)
    try:
        # Navigate via JS instead of the goto tool: NCBI's page pulls slow
        # third-party resources that can exceed lightpanda's navigation
        # timeout, but the form is usable long before the page fully settles.
        mcp.goto("about:blank")
        mcp.evaluate("location.href = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE=Proteins'")
        if not mcp.wait_for_script(
            "document.querySelector('textarea[name=\"QUERY\"]') !== null",
            timeout=120,
        ):
            error("❌ BLAST form did not load")
            return []

        submit_js = _build_blast_submit_js(fasta_text, evalue, dropdown_value, use_psiblast)

        rid = None
        if mcp.driver == "obscura":
            # obscura: native form.submit() and CDP clicks are unreliable in
            # this engine; submit the serialized form via a same-origin fetch
            # POST from the page (exactly the fields the form would send) and
            # read the RID straight from the response's hidden input.
            fill_only = submit_js.replace("form.submit();", "")
            try:
                mcp.evaluate(fill_only)
            except MCPBrowserError:
                error("❌ Could not fill the BLAST form")
                return []
            post_js = """fetch(document.querySelector('textarea[name="QUERY"]').form.action || location.href, {
                method: 'POST',
                headers: {'Content-Type': 'application/x-www-form-urlencoded'},
                body: new URLSearchParams(new FormData(document.querySelector('textarea[name="QUERY"]').form)).toString()
              }).then(r => r.text()).then(t => {
                const m = t.match(/name="RID"[^>]*value="([^"]+)"/);
                window.__rid = m ? m[1] : 'POSTED-NO-RID';
              }).catch(e => { window.__rid = 'FETCH-ERROR: ' + e; });"""
            with contextlib.suppress(MCPBrowserError):
                mcp.evaluate(post_js, timeout=90)  # the fetch keeps running in the page; the poll reads it
            for _attempt in range(40):
                time.sleep(2)
                try:
                    val = mcp.evaluate("String(window.__rid)")
                except MCPBrowserError:
                    continue
                if val and val not in ("null", "undefined"):
                    if val.startswith("FETCH-ERROR") or val == "POSTED-NO-RID":
                        error(f"❌ BLAST submission failed: {val}")
                        return []
                    rid = val
                    break
            if not rid:
                error("❌ Could not find RID (fetch POST did not complete)")
                return []
        else:
            # lightpanda: native form.submit() works; the RID is in the
            # status page's hidden input
            submitted = None
            for _attempt in range(5):
                try:
                    submitted = mcp.evaluate(submit_js)
                except MCPBrowserError:
                    submitted = None
                    time.sleep(1.5)
                    continue
                if isinstance(submitted, dict) and submitted.get("ok"):
                    break
                time.sleep(1.5)
            if not (isinstance(submitted, dict) and submitted.get("ok")):
                error(f"❌ Could not submit BLAST form: {submitted}")
                return []
            for _attempt in range(10):
                time.sleep(2)
                try:
                    rid = mcp.evaluate(
                        "document.querySelector('input[name=\"RID\"]')"
                        " ? document.querySelector('input[name=\"RID\"]').value : null"
                    )
                except MCPBrowserError:
                    continue  # evaluate during navigation can fail transiently
                if rid:
                    info(f"✓ RID: {rid}")
                    break

        if not rid:
            error("❌ Could not find RID")
            href = mcp.evaluate("String(location.href)")
            info(f"debug: Current URL: {href}")
            return []
    finally:
        mcp.close()

    status_url = (
        f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID={rid}&FORMAT_OBJECT=SearchInfo"
    )

    for _i in range(600):
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

    download_url = f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?RESULTS_FILE=on&RID={rid}&FORMAT_TYPE=CSV&DESCRIPTIONS={dropdown_value}&ALIGNMENT_VIEW=Tabular&CMD=Get"

    resp = requests.get(download_url)
    content = resp.text

    all_lines = content.strip().split("\n")
    data_lines = [line for line in all_lines if line and not line.startswith("#")]
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
