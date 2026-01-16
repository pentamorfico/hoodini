"""
Browser setup utilities for installing and verifying Playwright dependencies.

This module ensures Playwright Firefox is available before use.
"""

import asyncio
import subprocess
import sys


def ensure_playwright_firefox() -> bool:
    """
    Ensure Playwright Firefox is installed.

    Automatically installs Firefox if not already present.
    This should be called once per environment during initial setup.

    Returns:
        True if Firefox is available or was successfully installed, False otherwise.
    """
    try:
        # Try to import playwright first
        from playwright.async_api import async_playwright

        # Check if Firefox binary exists by attempting to launch it
        async def check_firefox():
            try:
                async with async_playwright() as p:
                    browser = await p.firefox.launch(headless=True)
                    await browser.close()
                    return True
            except Exception:
                return False

        if asyncio.run(check_firefox()):
            return True
    except Exception:
        pass

    # Firefox not available, attempt to install
    print("Installing Playwright Firefox browser...")
    try:
        subprocess.run(
            [sys.executable, "-m", "playwright", "install", "firefox"],
            check=True,
            capture_output=True,
            text=True,
        )
        print("✓ Playwright Firefox installed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"✗ Failed to install Playwright Firefox:\n{e.stderr}")
        return False


__all__ = ["ensure_playwright_firefox"]
