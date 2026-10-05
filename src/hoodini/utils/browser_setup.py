"""Browser setup utilities for ensuring a headless browser is available.

Hoodini drives obscura (conda-forge) or lightpanda (bioconda / pip wheel)
over MCP stdio — see ``mcp_browser``. No CDP server and no network port are
involved.
"""

from __future__ import annotations

from hoodini.utils.logging_utils import error
from hoodini.utils.mcp_browser import find_browser_binary


def ensure_lightpanda() -> bool:
    """Ensure an obscura or lightpanda browser binary is available."""
    if find_browser_binary():
        return True
    error(
        "✗ no headless browser found. Install one via: mamba install -c conda-forge obscura "
        "or: pip install lightpanda (bundles the browser), or: mamba install -c bioconda lightpanda"
    )
    return False


__all__ = ["ensure_lightpanda"]
