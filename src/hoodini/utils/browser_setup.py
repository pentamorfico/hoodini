"""Browser setup utilities for ensuring lightpanda is available.

This module only checks that the ``lightpanda`` headless browser binary is
installed. Browser sessions run over MCP stdio (see ``mcp_browser``), so no
CDP server and no network port are involved.
"""

from __future__ import annotations

from hoodini.utils.logging_utils import error
from hoodini.utils.mcp_browser import find_lightpanda_binary


def ensure_lightpanda() -> bool:
    """Ensure a lightpanda binary is available (MCP stdio mode needs no server)."""
    if find_lightpanda_binary():
        return True
    error(
        "✗ lightpanda binary not found. Install it via: pip install lightpanda "
        "(bundles the browser) or: mamba install -c bioconda lightpanda"
    )
    return False


__all__ = ["ensure_lightpanda"]
