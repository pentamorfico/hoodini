"""
Playwright utilities for browser automation with proper environment setup.

This module provides async helpers for Playwright Firefox that handle
LD_LIBRARY_PATH setup for conda/pixi/mamba environments automatically.
"""

from collections.abc import AsyncGenerator
from contextlib import asynccontextmanager

from playwright.async_api import Page, async_playwright

from .runtime_env import apply_ld_library_path


@asynccontextmanager
async def get_browser_context() -> AsyncGenerator[Page, None]:
    """
    Context manager for async Playwright browser with proper environment setup.
    
    Automatically:
    - Sets up LD_LIBRARY_PATH for conda/pixi/mamba environments
    - Launches Firefox browser with appropriate arguments
    - Provides a page context for automation
    - Cleans up browser on exit
    
    Yields:
        playwright.async_api.Page: Browser page for automation
        
    Example:
        >>> async with get_browser_context() as page:
        ...     await page.goto("https://example.com")
        ...     await page.fill("#search", "query")
    """
    # Apply environment setup before launching browser
    apply_ld_library_path()
    
    async with async_playwright() as p:
        browser = await p.firefox.launch(headless=True)
        page = await browser.new_page()
        try:
            yield page
        finally:
            await browser.close()


async def open_blast_page(url: str) -> tuple[Page, str | None]:
    """
    Open NCBI BLAST page and return page object for automation.
    
    This is a utility function for testing the Playwright setup.
    
    Args:
        url: The NCBI BLAST URL to open
        
    Returns:
        Tuple of (page object, error message if any)
        
    Example:
        >>> page, error = await open_blast_page("https://blast.ncbi.nlm.nih.gov/...")
        >>> if not error:
        ...     # Automate BLAST workflow
        ...     pass
    """
    try:
        apply_ld_library_path()
        async with async_playwright() as p:
            browser = await p.firefox.launch(headless=True)
            page = await browser.new_page()
            await page.goto(url)
            return page, None
    except Exception as e:
        return None, str(e)


__all__ = [
    "get_browser_context",
    "open_blast_page",
]
