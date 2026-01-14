import argparse
import shutil
from importlib.resources import files
from pathlib import Path

import requests
from rich.table import Table

from hoodini.utils.downloader import download_with_aria2c
from hoodini.utils.logging_utils import console, info, logger, stage_header, warn

PROJ = "3uz2j"
PROV = "osfstorage"
ROOT_URL = f"https://api.osf.io/v2/nodes/{PROJ}/files/{PROV}/"

DATA_DIR = files("hoodini").joinpath("data", "metacerberus")


def fetch_all_items(url):
    """Yield all items from a paginated OSF API endpoint."""
    while url:
        resp = requests.get(url)
        resp.raise_for_status()
        payload = resp.json()
        yield from payload["data"]
        url = payload["links"].get("next")


def list_db_files():
    db_id = None
    for item in fetch_all_items(ROOT_URL):
        attr = item["attributes"]
        if attr["kind"] == "folder" and attr["name"] == "db":
            db_id = item["id"]
            break
    if not db_id:
        raise RuntimeError("Couldn't find a folder named 'db' in root osfstorage")
    db_url = f"{ROOT_URL}{db_id}/"
    files = []
    for node in fetch_all_items(db_url):
        attr = node["attributes"]
        if attr["kind"] == "file":
            files.append(
                {
                    "name": attr["name"],
                    "download": node["links"]["download"],
                    "size": attr.get("size", None),
                }
            )
    return files


def get_db_groups(files):
    """Return a dict: group -> list of file dicts."""
    groups = {}
    for f in files:
        name = f["name"]
        if name.endswith(".hmm.gz") or name.endswith(".tsv"):
            group = name.split(".")[0].split("_")[0].lower()
            groups.setdefault(group, []).append(f)
    return groups


def check_downloaded(groups):
    """Return dict: group -> list of (file, present:bool)"""
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    status = {}
    for group, file_list in groups.items():
        status[group] = []
        for f in file_list:
            dest = DATA_DIR / f["name"]
            status[group].append((f, dest.exists()))
    return status


def download_files(files, force=False):
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    urls = [f["download"] for f in files]
    dests = [DATA_DIR / f["name"] for f in files]
    out_names = [f["name"] for f in files]
    logger.info(f"Downloading {len(urls)} metacerberus files to {DATA_DIR}")
    result_files = download_with_aria2c(urls, DATA_DIR, show_progress=True, out_names=out_names)

    for rf, dest in zip(result_files, dests):
        rpath = Path(rf)
        if rpath.is_file():
            if rpath != dest:
                shutil.move(str(rpath), str(dest))
        elif rpath.is_dir():
            candidate = rpath / dest.name
            if candidate.exists() and candidate.is_file():
                shutil.move(str(candidate), str(dest))
            else:
                files_inside = [p for p in rpath.iterdir() if p.is_file()]
                if files_inside:
                    shutil.move(str(files_inside[0]), str(dest))
                else:
                    raise RuntimeError(
                        f"No downloaded file found inside directory returned by downloader: {rpath}"
                    )
        else:
            raise RuntimeError(f"Downloaded path not found: {rpath}")
    logger.info("Metacerberus downloads complete")


def main(selected=None, force=False):
    stage_header("MetaCerberus Databases", "🧬")
    files = list_db_files()
    groups = get_db_groups(files)
    status = check_downloaded(groups)
    if selected is None or selected == "all":
        table = Table(title="MetaCerberus Databases", show_lines=True)
        table.add_column("Database", style="bold cyan")
        table.add_column("HMM file", style="green")
        table.add_column("TSV file", style="magenta")
        for group, file_statuses in sorted(status.items()):
            hmm = tsv = "[dim]-[/dim]"
            for f, present in file_statuses:
                if f["name"].endswith(".hmm.gz"):
                    if present:
                        hmm = f"[green]✔ {f['name']}[/green]"
                    else:
                        hmm = f"[red]✗ {f['name']}[/red]"
                elif f["name"].endswith(".tsv"):
                    if present:
                        tsv = f"[green]✔ {f['name']}[/green]"
                    else:
                        tsv = f"[red]✗ {f['name']}[/red]"
            table.add_row(group, hmm, tsv)
        console.print(table)
        return
    wanted = [s.strip().lower() for s in selected.split(",") if s.strip()]
    to_download = [f for g in wanted for f in groups.get(g, [])]
    if not to_download:
        warn(f"No files found for: {', '.join(wanted)}")
        return
    if not force:
        to_download = [f for f in to_download if not (DATA_DIR / f["name"]).exists()]
    if not to_download:
        info("All requested MetaCerberus files are present!")
        return
    download_files(to_download, force=force)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "selected",
        nargs="?",
        default="all",
        help="Which DB(s) to download: all, pfam, phrogs, etc. Comma-separated.",
    )
    parser.add_argument("--force", action="store_true", help="Overwrite existing files.")
    args = parser.parse_args()
    main(args.selected, force=args.force)
