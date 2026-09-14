import concurrent.futures
import contextlib
import re
import subprocess
from pathlib import Path

import requests
from rich.progress import BarColumn, Progress, TextColumn, TimeRemainingColumn, TransferSpeedColumn


def download_with_aria2c(
    urls,
    dest_dir,
    connections=16,
    split=16,
    show_progress=True,
    show_aria2c_output=False,
    out_names=None,
    num_threads: int = 0,
):
    """
    Download URLs to dest_dir using a single aria2c subprocess with Rich progress.
    Returns a list of downloaded file paths.
    """
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    results = []

    PERCENT_RE = re.compile(r"(\d+(?:\.\d+)?)%")
    SIZE_RE = re.compile(r"([\d.]+)\s*([KMGTP]?i?B)/([\d.]+)\s*([KMGTP]?i?B)")

    UNIT = {
        "B": 1,
        "KB": 1000,
        "MB": 1000**2,
        "GB": 1000**3,
        "TB": 1000**4,
        "PB": 1000**5,
        "KiB": 1024,
        "MiB": 1024**2,
        "GiB": 1024**3,
        "TiB": 1024**4,
        "PiB": 1024**5,
    }

    def to_bytes(num_str: str, unit_str: str) -> float:
        return float(num_str) * UNIT.get(unit_str.strip(), 1)

    def fmt_bytes(n: float) -> str:
        for unit in ["B", "KiB", "MiB", "GiB", "TiB", "PiB"]:
            if n < 1024 or unit == "PiB":
                return f"{n:.2f} {unit}"
            n /= 1024

    from urllib.parse import unquote, urlparse

    out_name_list = []
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
        out_name_list.append(out_name)

    max_conn = str(num_threads or 16)
    active_processes = set()

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
                task_ids = {
                    out_name: progress.add_task(out_name, total=None, bytes_text="…")
                    for out_name in out_name_list
                }

                def download_one(url, out_name):
                    cmd = [
                        "aria2c",
                        "--summary-interval=1",
                        "--enable-color=false",
                        "--max-connection-per-server",
                        max_conn,
                        "--split",
                        max_conn,
                        "-k",
                        "1M",
                        "-d",
                        str(dest_dir),
                        "--out",
                        out_name,
                        url,
                    ]
                    proc = subprocess.Popen(
                        cmd,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.STDOUT,
                        bufsize=0,
                    )
                    active_processes.add(proc)
                    buf = b""
                    total_bytes = None
                    task = task_ids[out_name]

                    while True:
                        chunk = proc.stdout.read(1024)
                        if not chunk:
                            if proc.poll() is not None:
                                break
                            continue

                        buf += chunk
                        parts = re.split(rb"[\r\n]", buf)
                        for part in parts[:-1]:
                            line = part.decode("utf-8", "ignore")
                            downloaded, total, pct = None, None, None
                            match = SIZE_RE.search(line)
                            if match:
                                d_num, d_unit, t_num, t_unit = match.groups()
                                downloaded = to_bytes(d_num, d_unit)
                                total = to_bytes(t_num, t_unit)
                            match = PERCENT_RE.search(line)
                            if match:
                                pct = float(match.group(1))

                            if show_aria2c_output and (
                                downloaded is not None or total is not None or pct is not None
                            ):
                                progress.console.print(line)

                            if total is not None and downloaded is not None:
                                total_bytes = int(total)
                                progress.update(
                                    task,
                                    total=total_bytes,
                                    completed=int(downloaded),
                                    bytes_text=f"{fmt_bytes(downloaded)} / {fmt_bytes(total)}",
                                )
                            elif pct is not None and total_bytes is None:
                                progress.update(
                                    task, total=100.0, completed=pct, bytes_text=f"{pct:.1f}%"
                                )
                        buf = parts[-1]

                    code = proc.wait()
                    active_processes.discard(proc)
                    if code != 0:
                        raise RuntimeError(f"aria2c exited with non-zero status: {code}")
                    if total_bytes:
                        progress.update(
                            task,
                            total=total_bytes,
                            completed=total_bytes,
                            bytes_text=f"{fmt_bytes(total_bytes)} / {fmt_bytes(total_bytes)}",
                        )
                    else:
                        progress.update(task, total=100.0, completed=100.0, bytes_text="100%")

                with concurrent.futures.ThreadPoolExecutor(max_workers=len(urls)) as executor:
                    futures = [
                        executor.submit(download_one, url, out_name)
                        for url, out_name in zip(urls, out_name_list)
                    ]
                    for future in futures:
                        future.result()
        else:
            input_lines = []
            for url, out_name in zip(urls, out_name_list):
                input_lines.append(f"{url}\n  out={out_name}")
            from tempfile import NamedTemporaryFile

            with NamedTemporaryFile("w", delete=False) as f:
                for line in input_lines:
                    f.write(line + "\n")
                input_file = f.name
            try:
                subprocess.run(
                    [
                        "aria2c",
                        "--summary-interval=1",
                        "--enable-color=false",
                        "--max-connection-per-server",
                        max_conn,
                        "--split",
                        max_conn,
                        "-k",
                        "1M",
                        "-d",
                        str(dest_dir),
                        "-i",
                        input_file,
                    ],
                    check=True,
                )
            finally:
                with contextlib.suppress(Exception):
                    Path(input_file).unlink()
    finally:
        for proc in active_processes:
            if proc.poll() is None:
                proc.terminate()

    for out_name in out_name_list:
        candidate = dest_dir / out_name
        if candidate.exists():
            results.append(str(candidate))
    return results
