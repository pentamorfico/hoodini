"""Thin wrapper around NCBI E-utilities' efetch, replacing the entrez-direct
``efetch`` CLI binary with a direct HTTP call.

``efetch`` (from Bioconda's entrez-direct package) is itself just a bash
script that shells out to ``nquire``, which does:

    curl -X POST https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi \
        -d "db=<db>&id=<ids>&rettype=<rettype>&retmode=<retmode>" \
        -H "api-key: <NCBI_API_KEY>"   # only if set

i.e. a POST to the public efetch.fcgi endpoint with the id list in the POST
body (never a GET query string, which matters once you're batching hundreds
of accessions per request) and the API key as an HTTP header rather than a
query parameter. This module does exactly that directly, dropping the need
for the entrez-direct conda package as a runtime dependency.
"""

from __future__ import annotations

import os
import time

import requests

EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def efetch(
    db: str,
    ids: list[str] | str,
    rettype: str,
    retmode: str = "text",
    api_key: str | None = None,
    tool: str = "hoodini",
    timeout: int = 90,
) -> str:
    """Fetch records from NCBI via E-utilities' efetch.fcgi.

    Mirrors ``efetch -db <db> -id <ids> -format <rettype> -mode <retmode>``.
    Raises ``requests.HTTPError`` on non-2xx responses; callers are expected
    to handle retries themselves (NCBI's efetch endpoint occasionally 500s
    under load, same as the CLI binary did).
    """
    id_str = ",".join(ids) if isinstance(ids, list) else ids
    data = {
        "db": db,
        "id": id_str,
        "rettype": rettype,
        "retmode": retmode,
        "tool": tool,
    }
    key = api_key or os.environ.get("NCBI_API_KEY")
    headers = {"api-key": key} if key else {}
    resp = requests.post(EFETCH_URL, data=data, headers=headers, timeout=timeout)
    resp.raise_for_status()
    return resp.text


def efetch_with_retries(
    db: str,
    ids: list[str] | str,
    rettype: str,
    retmode: str = "text",
    api_key: str | None = None,
    tool: str = "hoodini",
    timeout: int = 90,
    max_retries: int = 3,
    backoff_seconds: float = 5.0,
) -> str:
    """``efetch`` with simple retry-on-failure, matching the old CLI wrapper's
    behavior of retrying transient 500s/timeouts/empty responses."""
    last_error: Exception | None = None
    for attempt in range(max_retries):
        try:
            text = efetch(
                db,
                ids,
                rettype,
                retmode=retmode,
                api_key=api_key,
                tool=tool,
                timeout=timeout,
            )
            if text and text.strip():
                return text
            last_error = RuntimeError("efetch returned an empty response")
        except (requests.RequestException, RuntimeError) as e:
            last_error = e
        if attempt < max_retries - 1:
            time.sleep(backoff_seconds)
    raise RuntimeError(f"efetch failed after {max_retries} attempts: {last_error}")
