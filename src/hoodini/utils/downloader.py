import contextlib
import re
import time
from pathlib import Path
from urllib.parse import unquote, urlparse

import bytehaul
import requests
from rich.progress import BarColumn, Progress, TextColumn, TimeRemainingColumn, TransferSpeedColumn

_TERMINAL_STATES = {"completed", "failed", "cancelled", "paused"}


def _resolve_out_names(urls, out_names):
    """Determine the output filename for each URL (explicit, Content-Disposition, or URL path)."""
    resolved = []
    for idx, url in enumerate(urls):
        out_name = None
        if out_names and idx < len(out_names) and out_names[idx]:
            out_name = out_names[idx]
        if not out_name:
            parsed = urlparse(url)
            out_name = unquote(Path(parsed.path).name)
            if not out_name:
                try:
                    r = requests.head(url, allow_redirects=True, timeout=5)
                    cd = r.headers.get("content-disposition")
                    if cd:
                        m = re.search(r"filename\*?=(?:UTF-8'')?\"?([^\";]+)\"?", cd)
                        if m:
                            out_name = m.group(1)
                except Exception:
                    out_name = ""
        if not out_name:
            out_name = f"downloaded_file_{idx}"
        resolved.append(out_name)
    return resolved


def _fmt_bytes(n: float) -> str:
    for unit in ["B", "KiB", "MiB", "GiB", "TiB", "PiB"]:
        if n < 1024 or unit == "PiB":
            return f"{n:.2f} {unit}"
        n /= 1024


def download_urls(
    urls,
    dest_dir,
    connections=16,
    show_progress=True,
    out_names=None,
    num_threads: int = 0,
):
    """
    Download URLs to dest_dir concurrently using bytehaul, with Rich progress.

    Returns a list of downloaded file paths.
    """
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)

    out_name_list = _resolve_out_names(urls, out_names)
    max_conn = num_threads or connections or 16

    downloader = bytehaul.Downloader()
    tasks = {
        out_name: downloader.download(
            url,
            output_dir=str(dest_dir),
            output_path=out_name,
            max_connections=max_conn,
        )
        for url, out_name in zip(urls, out_name_list)
    }

    try:
        if show_progress:
            with Progress(
                TextColumn("[bold blue]{task.description}"),
                BarColumn(),
                TextColumn("[progress.percentage]{task.percentage:>5.1f}%"),
                TextColumn("• {task.fields[bytes_text]}"),
                TransferSpeedColumn(),
                TimeRemainingColumn(),
                refresh_per_second=10,
                transient=True,
            ) as progress:
                progress_task_ids = {
                    out_name: progress.add_task(out_name, total=None, bytes_text="…")
                    for out_name in out_name_list
                }

                pending = set(out_name_list)
                while pending:
                    for out_name in list(pending):
                        snap = tasks[out_name].progress()
                        progress_task = progress_task_ids[out_name]
                        if snap.total_size:
                            progress.update(
                                progress_task,
                                total=snap.total_size,
                                completed=snap.downloaded,
                                bytes_text=(
                                    f"{_fmt_bytes(snap.downloaded)} / {_fmt_bytes(snap.total_size)}"
                                ),
                            )
                        else:
                            progress.update(
                                progress_task,
                                completed=snap.downloaded,
                                bytes_text=_fmt_bytes(snap.downloaded),
                            )
                        if snap.state in _TERMINAL_STATES:
                            pending.discard(out_name)
                    if pending:
                        time.sleep(0.1)
        for task in tasks.values():
            task.wait()
    except BaseException:
        for task in tasks.values():
            with contextlib.suppress(Exception):
                task.cancel()
        raise

    results = []
    for out_name in out_name_list:
        candidate = dest_dir / out_name
        if candidate.exists():
            results.append(str(candidate))
    return results
