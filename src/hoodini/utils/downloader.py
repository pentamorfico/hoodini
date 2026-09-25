import contextlib
import re
import time
from pathlib import Path
from urllib.parse import unquote, urlparse

import bytehaul
import requests
from rich.progress import BarColumn, Progress, TextColumn, TimeRemainingColumn, TransferSpeedColumn

from hoodini.utils.logging_utils import info, warn

_TERMINAL_STATES = {"completed", "failed", "cancelled", "paused"}
_SUCCESS_STATES = {"completed"}
_CHUNK_SIZE = 1 << 20


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


def _is_osf_url(url: str) -> bool:
    """OSF/GCS serves large files fine but bytehaul's multi-connection strategy
    deterministically fails there with HTTP 400, so OSF is downloaded with a
    single connection."""
    host = (urlparse(url).netloc or "").split(":")[0]
    return host == "osf.io" or host.endswith(".osf.io")


def _single_stream_download(
    url,
    dest_path: Path,
    *,
    max_retries=8,
    retry_base_delay=2.0,
    retry_max_delay=60.0,
    max_retry_elapsed=300.0,
    progress=None,
    progress_task=None,
):
    """Download a URL with a single requests stream, with exponential backoff.

    Fallback for servers that reject the multi-connection/range strategy used
    by bytehaul (e.g. OSF files larger than ~500 MB fail there with HTTP 400).
    """
    start = time.monotonic()
    delay = retry_base_delay
    last_error = None
    tmp_path = dest_path.with_suffix(dest_path.suffix + ".part")
    for attempt in range(1, max_retries + 1):
        try:
            with requests.get(url, stream=True, timeout=(30, 120), allow_redirects=True) as r:
                r.raise_for_status()
                if progress is not None and progress_task is not None:
                    total = r.headers.get("content-length")
                    if total:
                        progress.update(progress_task, total=int(total), completed=0)
                with open(tmp_path, "wb") as f:
                    for chunk in r.iter_content(chunk_size=_CHUNK_SIZE):
                        if chunk:
                            f.write(chunk)
                            if progress is not None and progress_task is not None:
                                progress.update(progress_task, completed=f.tell())
            tmp_path.replace(dest_path)
            return
        except Exception as e:
            last_error = e
            with contextlib.suppress(Exception):
                tmp_path.unlink(missing_ok=True)
            if attempt >= max_retries or (time.monotonic() - start) > max_retry_elapsed:
                break
            time.sleep(delay)
            delay = min(delay * 2, retry_max_delay)
    raise RuntimeError(f"single-stream download failed: {last_error}")


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
    Download URLs to dest_dir with bytehaul, with Rich progress.

    OSF URLs are downloaded with a single connection because bytehaul's
    multi-connection strategy deterministically fails there with HTTP 400 on
    files larger than ~500 MB. Any bytehaul failure falls back to a
    single-connection requests stream, and the (possibly corrupt,
    preallocated) artifact bytehaul leaves behind is removed first.

    Defaults are deliberately conservative: some servers (e.g. NCBI's FTP/HTTPS
    mirrors) throttle or return HTTP 503 once too many concurrent connections
    per IP are open. Downloading many URLs at once, each with its own pool of
    connections, easily stacks up dozens of simultaneous connections, so at
    most `max_concurrent_files` files are ever in flight together, and each
    keeps a modest number of connections. The retry knobs give transient
    rate-limiting time to clear.

    Raises RuntimeError if any URL could not be downloaded. Otherwise returns
    a list of downloaded file paths.
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
            max_connections=1 if _is_osf_url(url) else max_conn,
            max_retries=max_retries,
            retry_base_delay=retry_base_delay,
            retry_max_delay=retry_max_delay,
            max_retry_elapsed=max_retry_elapsed,
        )

    items = list(zip(urls, out_name_list))
    tasks: dict[str, object] = {}
    all_failures: list[str] = []

    def _cleanup_artifacts(dest: Path):
        """Remove the (possibly corrupt) output bytehaul left behind."""
        for leftover in (dest, dest.parent / f"{dest.name}.bytehaul"):
            with contextlib.suppress(Exception):
                leftover.unlink(missing_ok=True)

    def _run_batch(batch, progress=None, progress_task_ids=None):
        batch_tasks = {out_name: _start(url, out_name) for url, out_name in batch}
        tasks.update(batch_tasks)
        pending = set(batch_tasks)
        final_states: dict[str, str] = {}
        final_sizes: dict[str, int] = {}
        while pending:
            for out_name in list(pending):
                snap = batch_tasks[out_name].progress()
                final_states[out_name] = snap.state
                final_sizes[out_name] = snap.total_size or 0
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

        failures = []
        url_by_name = {out_name: url for url, out_name in batch}
        for out_name, task in batch_tasks.items():
            dest = dest_dir / out_name
            error = None
            try:
                # wait() consumes the task: bytehaul raises on ANY further
                # access (even progress()), so the state must be read from
                # final_states, captured while polling above.
                task.wait()
            except Exception as e:
                error = str(e) or type(e).__name__
            else:
                state = final_states.get(out_name, "")
                if state not in _SUCCESS_STATES:
                    error = f"bytehaul ended in state '{state}'"
                elif not dest.exists() or dest.stat().st_size == 0:
                    error = "bytehaul reported success but the output file is missing or empty"
                elif final_sizes.get(out_name) and dest.stat().st_size != final_sizes[out_name]:
                    error = (
                        f"size mismatch: got {dest.stat().st_size}, "
                        f"expected {final_sizes[out_name]}"
                    )

            if error is None:
                continue

            # bytehaul preallocates the output file, so a failed download can
            # leave a corrupt artifact with the full expected size behind;
            # remove it before retrying so it is not mistaken for a valid file.
            _cleanup_artifacts(dest)
            warn(
                f"bytehaul failed for {out_name} ({error}); "
                "retrying with a single-connection download..."
            )
            try:
                _single_stream_download(
                    url_by_name[out_name],
                    dest,
                    max_retries=max_retries,
                    retry_base_delay=retry_base_delay,
                    retry_max_delay=retry_max_delay,
                    max_retry_elapsed=max_retry_elapsed,
                )
                info(
                    f"[green]✔[/green] {out_name} downloaded successfully (single-connection fallback)"
                )
            except Exception as fallback_error:
                _cleanup_artifacts(dest)
                failures.append(f"{out_name}: {error}; fallback: {fallback_error}")

        if failures:
            raise RuntimeError("; ".join(failures))

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
                    try:
                        _run_batch(items[start : start + batch_size], progress, progress_task_ids)
                    except RuntimeError as e:
                        all_failures.append(str(e))
        else:
            for start in range(0, len(items), batch_size):
                try:
                    _run_batch(items[start : start + batch_size])
                except RuntimeError as e:
                    all_failures.append(str(e))
    except BaseException:
        for task in tasks.values():
            with contextlib.suppress(Exception):
                task.cancel()
        raise

    if all_failures:
        raise RuntimeError("; ".join(all_failures))

    results = []
    for out_name in out_name_list:
        candidate = dest_dir / out_name
        if candidate.exists():
            results.append(str(candidate))
    return results
