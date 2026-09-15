"""Browser setup utilities for ensuring lightpanda is available.

This module ensures the ``lightpanda`` headless browser binary is installed
and its CDP server is running before use. Unlike Playwright/Firefox, lightpanda
is a single static binary (from bioconda) with no GTK/X11 dependencies, so
"ensuring" it is available just means checking the binary exists and starting
its ``serve`` process if needed.
"""

from hoodini.utils.cdp_browser import (
    find_lightpanda_binary,
    lightpanda_server_alive,
    start_lightpanda_server,
)
from hoodini.utils.logging_utils import error, info


def ensure_lightpanda() -> bool:
    """
    Ensure the lightpanda CDP server is available, starting it if needed.

    Returns:
        True if the lightpanda server is running (or was successfully
        started), False otherwise.
    """
    if lightpanda_server_alive():
        return True

    lp_bin = find_lightpanda_binary()
    if not lp_bin:
        error(
            "✗ lightpanda binary not found. Install it via: " "mamba install -c bioconda lightpanda"
        )
        return False

    try:
        info("Starting lightpanda CDP server...")
        start_lightpanda_server(lp_bin)
        info("✓ lightpanda server is running")
        return True
    except Exception as e:
        error(f"✗ Failed to start lightpanda server: {e}")
        return False


__all__ = ["ensure_lightpanda"]
