"""Download and manage the UniProt ID-mapping parquet database."""

from __future__ import annotations

from importlib.resources import files
from pathlib import Path

from hoodini.utils.logging_utils import info, warn

REMOTE_IDMAPPING_URL = "https://storage.hoodini.bio/idmapping_selected.parquet"
DATA_DIR = files("hoodini").joinpath("data")
IDMAPPING_PARQUET = DATA_DIR.joinpath("idmapping_selected.parquet")


def idmapping_db_exists() -> bool:
    """Return True if the local idmapping parquet file exists."""
    return Path(IDMAPPING_PARQUET).exists()


def get_idmapping_path() -> Path:
    """Return the path to the local idmapping parquet."""
    return Path(IDMAPPING_PARQUET)


def download_idmapping(dest: Path | None = None, num_threads: int = 0) -> bool:
    """Download idmapping_selected.parquet from remote storage.

    Parameters
    ----------
    dest : Path, optional
        Destination path.  Defaults to the package data directory.
    num_threads : int
        Thread count forwarded to aria2c (0 = auto).

    Returns
    -------
    bool
        True on success.
    """
    from hoodini.download.databases import _download_url

    if dest is None:
        dest = get_idmapping_path()

    dest.parent.mkdir(parents=True, exist_ok=True)
    info(f"Downloading UniProt ID-mapping database → {dest}")
    info("   This may take a few minutes...")

    ok = _download_url(REMOTE_IDMAPPING_URL, dest, num_threads=num_threads)
    if ok and dest.exists():
        info(f"✔️  Downloaded ID-mapping database to {dest}")
        return True

    warn(
        "⚠️  Failed to download idmapping_selected.parquet. "
        "Run 'hoodini download idmapping' manually to retry."
    )
    return False
