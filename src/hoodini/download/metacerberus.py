import argparse
import gzip
import shutil
import tempfile
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

# Pfam direct download URLs
PFAM_HMM_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.1/Pfam-A.hmm.gz"
PFAM_DAT_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam38.1/Pfam-A.hmm.dat.gz"


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


def parse_pfam_dat_to_tsv(dat_gz_path: Path, tsv_path: Path):
    """
    Parse Pfam-A.hmm.dat.gz Stockholm format and convert to TSV.
    
    Input format (Stockholm):
        # STOCKHOLM 1.0
        #=GF ID   2-Hacid_dh_C
        #=GF AC   PF02826.26
        #=GF DE   D-isomer specific 2-hydroxyacid dehydrogenase, NAD binding domain
        #=GF CL   CL0063
        //
    
    Output format (TSV):
        ID	Function	Clan	Accession
        2-Hacid_dh_C	D-isomer specific 2-hydroxyacid dehydrogenase, NAD binding domain	CL0063	PF02826
    
    Note: ID column uses the Gene name (from #=GF ID) because HMMER uses NAME field
    as domain_id, which corresponds to the gene name, not the PF accession.
    """
    logger.info(f"Parsing Pfam dat file: {dat_gz_path}")
    
    entries = []
    current = {}
    
    with gzip.open(dat_gz_path, "rt", encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#=GF ID"):
                # Gene name - this will be our ID for joining
                current["ID"] = line.split(None, 2)[2].strip()
            elif line.startswith("#=GF AC"):
                # Accession - strip version number (PF02826.26 -> PF02826)
                acc = line.split(None, 2)[2].strip()
                current["Accession"] = acc.split(".")[0]
            elif line.startswith("#=GF DE"):
                # Description/Function
                current["Function"] = line.split(None, 2)[2].strip()
            elif line.startswith("#=GF CL"):
                # Clan
                current["Clan"] = line.split(None, 2)[2].strip()
            elif line == "//":
                # End of entry
                if "ID" in current:
                    entries.append({
                        "ID": current.get("ID", ""),
                        "Function": current.get("Function", ""),
                        "Clan": current.get("Clan", ""),
                        "Accession": current.get("Accession", ""),
                    })
                current = {}
    
    # Write TSV
    with open(tsv_path, "w", encoding="utf-8") as out:
        out.write("ID\tFunction\tClan\tAccession\n")
        for entry in entries:
            out.write(f"{entry['ID']}\t{entry['Function']}\t{entry['Clan']}\t{entry['Accession']}\n")
    
    logger.info(f"Wrote {len(entries)} Pfam entries to {tsv_path}")


def download_pfam_direct(force=False):
    """Download Pfam HMM and dat files directly from EBI FTP and create TSV."""
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    
    hmm_dest = DATA_DIR / "Pfam.hmm.gz"
    tsv_dest = DATA_DIR / "Pfam.tsv"
    
    # Check if already downloaded
    if not force and hmm_dest.exists() and tsv_dest.exists():
        info("Pfam files already present!")
        return
    
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        
        # Download HMM file
        if force or not hmm_dest.exists():
            logger.info(f"Downloading Pfam HMM from {PFAM_HMM_URL}")
            download_with_aria2c([PFAM_HMM_URL], tmpdir, show_progress=True)
            downloaded_hmm = tmpdir / "Pfam-A.hmm.gz"
            if downloaded_hmm.exists():
                shutil.move(str(downloaded_hmm), str(hmm_dest))
                logger.info(f"Saved HMM to {hmm_dest}")
        
        # Download dat file and convert to TSV
        if force or not tsv_dest.exists():
            logger.info(f"Downloading Pfam dat from {PFAM_DAT_URL}")
            download_with_aria2c([PFAM_DAT_URL], tmpdir, show_progress=True)
            downloaded_dat = tmpdir / "Pfam-A.hmm.dat.gz"
            if downloaded_dat.exists():
                parse_pfam_dat_to_tsv(downloaded_dat, tsv_dest)
    
    info("Pfam download complete!")


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
    
    # Remove pfam from OSF groups - we'll handle it separately
    groups.pop("pfam", None)
    
    # Add pfam with local file paths for status display
    pfam_hmm = DATA_DIR / "Pfam.hmm.gz"
    pfam_tsv = DATA_DIR / "Pfam.tsv"
    groups["pfam"] = [
        {"name": "Pfam.hmm.gz", "download": PFAM_HMM_URL, "size": None},
        {"name": "Pfam.tsv", "download": PFAM_DAT_URL, "size": None},
    ]
    
    status = check_downloaded(groups)
    
    if selected is None or selected == "all":
        table = Table(title="MetaCerberus Databases", show_lines=True)
        table.add_column("Database", style="bold cyan")
        table.add_column("HMM file", style="green")
        table.add_column("TSV file", style="magenta")
        table.add_column("Source", style="dim")
        for group, file_statuses in sorted(status.items()):
            hmm = tsv = "[dim]-[/dim]"
            source = "[dim]OSF[/dim]" if group != "pfam" else "[cyan]EBI/Pfam[/cyan]"
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
            table.add_row(group, hmm, tsv, source)
        console.print(table)
        return
    
    wanted = [s.strip().lower() for s in selected.split(",") if s.strip()]
    
    # Handle pfam separately
    if "pfam" in wanted:
        download_pfam_direct(force=force)
        wanted.remove("pfam")
    
    # Handle other databases from OSF
    if wanted:
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
