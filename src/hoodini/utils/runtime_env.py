"""
Runtime environment utilities for managing LD_LIBRARY_PATH across conda/pixi/mamba.

This module handles automatic detection and setup of library paths needed by
Playwright Firefox in conda/mamba/pixi environments without requiring sudo.
"""

import os
from pathlib import Path


def find_candidate_lib_dirs() -> list[str]:
    """
    Find candidate library directories from conda/mamba/pixi environments.

    Searches for:
    - $CONDA_PREFIX/lib (conda/mamba active environment)
    - ~/.pixi/envs/*/lib (pixi environments)
    - ~/.mamba/envs/*/lib (mamba environments)

    Returns:
        List of absolute paths to lib directories, deduplicated and ordered.
    """
    candidates: list[str] = []

    # Check CONDA_PREFIX (active conda/mamba environment)
    if conda_prefix := os.environ.get("CONDA_PREFIX"):
        lib_dir = os.path.join(conda_prefix, "lib")
        if os.path.isdir(lib_dir):
            candidates.append(lib_dir)

    # Check pixi environments
    pixi_envs_dir = Path.home() / ".pixi" / "envs"
    if pixi_envs_dir.exists():
        for env_dir in pixi_envs_dir.iterdir():
            if env_dir.is_dir():
                lib_dir = env_dir / "lib"
                if lib_dir.exists():
                    candidates.append(str(lib_dir))

    # Check mamba environments
    mamba_envs_dir = Path.home() / ".mamba" / "envs"
    if mamba_envs_dir.exists():
        for env_dir in mamba_envs_dir.iterdir():
            if env_dir.is_dir():
                lib_dir = env_dir / "lib"
                if lib_dir.exists():
                    candidates.append(str(lib_dir))

    # Deduplicate while preserving order
    seen = set()
    unique_candidates = []
    for path in candidates:
        if path not in seen:
            seen.add(path)
            unique_candidates.append(path)

    return unique_candidates


def verify_gtk_availability() -> str | None:
    """
    Verify that GTK3 libraries are available in the current environment.

    Returns:
        Path to the lib directory containing libgtk-3.so.0, or None if not found.
    """
    for lib_dir in find_candidate_lib_dirs():
        gtk_lib = Path(lib_dir) / "libgtk-3.so.0"
        if gtk_lib.exists():
            return lib_dir
    return None


def apply_ld_library_path() -> None:
    """
    Apply LD_LIBRARY_PATH modifications for Playwright Firefox in conda/pixi/mamba.

    Prepends candidate library directories to LD_LIBRARY_PATH so that Playwright
    Firefox can find GTK3 and other required libraries.

    This should be called before launching Playwright browser instances.
    """
    candidates = find_candidate_lib_dirs()
    if not candidates:
        # No conda/pixi/mamba environments found, nothing to do
        return

    current_ld_path = os.environ.get("LD_LIBRARY_PATH", "")

    # Build new LD_LIBRARY_PATH: candidates first, then existing path
    ld_path_parts = candidates.copy()
    if current_ld_path:
        ld_path_parts.append(current_ld_path)

    new_ld_path = ":".join(ld_path_parts)
    os.environ["LD_LIBRARY_PATH"] = new_ld_path


__all__ = [
    "find_candidate_lib_dirs",
    "verify_gtk_availability",
    "apply_ld_library_path",
]
