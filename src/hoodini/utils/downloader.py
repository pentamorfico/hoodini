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
    *,
    connections=4,
    show_progress=True,
    out_names=None,
    num_threads: int = 0,
    max_retries=8,
    retry_base_delay=2.0,
    retry_max_delay=60.0,
    max_retry_elapsed=300.0,
    max_concurrent_files=2,
):
    """
    Download URLs to dest_dir using bytehaul, with Rich progress.

    Defaults are deliberately conservative: some servers (e.g. NCBI's FTP/HTTPS
    mirrors) throttle or return HTTP 503 once too many concurrent connections
    per IP are open. Downloading many URLs at once, each with its own pool of
    connections, easily stacks up dozens of simultaneous connections, so at
    most `max_concurrent_files` files are ever in flight together, and each
    keeps a modest number of connections. The retry knobs give transient
    rate-limiting time to clear.

    Returns a list of downloaded file paths.
    """
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)

    out_name_list = _resolve_out_names(urls, out_names)
    max_conn = num_threads or connections or 4
    batch_size = max(1, max_concurrent_files)

    downloader = bytehaul.Downloader()

    def _start(url, out_name):
        return downloader.download(
            url,
            output_dir=str(dest_dir),
            output_path=out_name,
            max_connections=max_conn,
            max_retries=max_retries,
            retry_base_delay=retry_base_delay,
            retry_max_delay=retry_max_delay,
            max_retry_elapsed=max_retry_elapsed,
        )

    items = list(zip(urls, out_name_list))
    tasks: dict[str, object] = {}

    def _run_batch(batch, progress=None, progress_task_ids=None):
        batch_tasks = {out_name: _start(url, out_name) for url, out_name in batch}
        tasks.update(batch_tasks)
        pending = set(batch_tasks)
        while pending:
            for out_name in list(pending):
                snap = batch_tasks[out_name].progress()
                if progress is not None:
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
        for task in batch_tasks.values():
            task.wait()

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
                    out_name: progress.add_task(out_name, total=None, bytes_text="queued")
                    for _, out_name in items
                }
                for start in range(0, len(items), batch_size):
                    _run_batch(items[start : start + batch_size], progress, progress_task_ids)
        else:
            for start in range(0, len(items), batch_size):
                _run_batch(items[start : start + batch_size])
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
