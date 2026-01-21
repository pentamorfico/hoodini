import argparse
import shutil
import subprocess
import tarfile
from importlib.resources import files
from pathlib import Path

from hoodini.utils.downloader import download_with_aria2c
from hoodini.utils.logging_utils import error, info, stage_done, stage_header, warn

EMAPPER_URL = "http://eggnog6.embl.de/download/emapperdb-5.0.2/mmseqs.tar.gz"
EGGNOG_OG = "https://storage.hoodini.bio/eggnog_og.parquet"
EGGNOG_PROTS = "https://storage.hoodini.bio/eggnog_prots.parquet"
CONTIGS_URL = "https://storage.hoodini.bio/contig_lengths.parquet"


def _padloc_db_exists() -> bool:
    """Check if padloc database is installed by looking for HMM files.
    
    Padloc stores its database relative to its binary: $(dirname $(which padloc))/../data
    """
    padloc_bin = shutil.which("padloc")
    if padloc_bin:
        # padloc stores data in ../data relative to the binary
        padloc_data = Path(padloc_bin).resolve().parent.parent / "data" / "hmm"
        if padloc_data.exists() and any(padloc_data.iterdir()):
            return True
    return False


def _defensefinder_db_exists() -> bool:
    """Check if defense-finder models are installed.
    
    Defense-finder stores models in ~/.macsyfinder/models/defense-finder-models/
    """
    models_dir = Path.home() / ".macsyfinder" / "models" / "defense-finder-models"
    return models_dir.exists() and any(models_dir.iterdir())


def _genomad_db_exists(genomad_db: Path) -> bool:
    """Check if genomad database is installed."""
    # genomad creates a nested genomad_db/genomad_db structure
    inner_db = genomad_db / "genomad_db"
    if inner_db.exists():
        # Check for key database files
        return (inner_db / "genomad_db").exists() or (inner_db / "genomad_db.index").exists()
    return False


def _run_cmd(cmd, cwd: Path | None = None):
    try:
        info(f"Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True, cwd=(str(cwd) if cwd is not None else None))
        return True
    except FileNotFoundError:
        warn(f"Command not found: {cmd[0]} — skipping.")
        return False
    except subprocess.CalledProcessError as e:
        warn(f"Command failed ({cmd[0]}): {e} — continuing.")
        return False


def _download_url(url: str, dest: Path, num_threads: int = 0):
    dest.parent.mkdir(parents=True, exist_ok=True)
    info(f"Downloading {url} -> {dest}")
    try:
        out_name = Path(dest).name
        result_files = download_with_aria2c(
            [url], dest.parent, show_progress=True, out_names=[out_name], num_threads=num_threads
        )
        if any(str(dest) == f for f in result_files):
            return True
        for f in result_files:
            pf = Path(f)
            if pf.name == out_name and pf != dest:
                shutil.move(str(pf), str(dest))
                return True
            if pf.is_dir():
                candidate = pf / out_name
                if candidate.exists() and candidate.is_file():
                    shutil.move(str(candidate), str(dest))
                    return True
                files_inside = [p for p in pf.iterdir() if p.is_file()]
                if files_inside:
                    shutil.move(str(files_inside[0]), str(dest))
                    return True
        return False
    except Exception as e:
        warn(f"Download failed: {e}")
        return False


def extract_tar(tar_path: Path, dest_dir: Path, threads: int = 0) -> bool:
    """Extract tar.gz using pigz if available (fast), else fall back to tar, else Python tarfile."""
    pigz = shutil.which("pigz")
    tar_bin = shutil.which("tar")

    if pigz and tar_bin:
        pigz_cmd = f"pigz -p {threads}" if threads > 0 else "pigz"
        cmd = [
            tar_bin,
            f"--use-compress-program={pigz_cmd}",
            "-xvf",
            str(tar_path),
            "-C",
            str(dest_dir),
        ]
        try:
            info(f"Running: {' '.join(cmd)}")
            subprocess.run(cmd, check=True)
            return True
        except subprocess.CalledProcessError as e:
            warn(f"pigz extraction failed: {e}, falling back to tarfile")

    if tar_bin:
        try:
            subprocess.run([tar_bin, "-xzf", str(tar_path), "-C", str(dest_dir)], check=True)
            return True
        except subprocess.CalledProcessError as e:
            warn(f"tar extraction failed: {e}, falling back to tarfile")

    try:
        with tarfile.open(tar_path, "r:gz") as tf:
            tf.extractall(path=dest_dir)
        return True
    except Exception as e:
        error(f"tarfile extraction failed: {e}")
        return False


def main(
    force: bool = False,
    skip_padloc: bool = False,
    skip_deffinder: bool = False,
    skip_genomad: bool = False,
    skip_emapper: bool = False,
    skip_parquet: bool = False,
    skip_contig_lengths: bool = False,
    num_threads: int = 0,
):
    stage_header("Downloading databases and support files", "⬇️")

    data_dir = files("hoodini").joinpath("data")
    emapper_dir = data_dir.joinpath("emapper")
    contig_dir = data_dir.joinpath("contig_lengths")
    genomad_db = data_dir.joinpath("genomad_db")

    if not skip_padloc:
        if force or not _padloc_db_exists():
            info("Updating padloc database...")
            _run_cmd(["padloc", "--db-update"])
        else:
            info("padloc database already exists; use --force to re-download")
    else:
        info("Skipping padloc DB update (--skip-padloc)")

    if not skip_deffinder:
        if force or not _defensefinder_db_exists():
            info("Updating defense-finder models...")
            _run_cmd(["defense-finder", "update"])
        else:
            info("defense-finder models already exist; use --force to re-download")
    else:
        info("Skipping defense-finder model install (--skip-deffinder)")

    if not skip_genomad:
        genomad_db.mkdir(parents=True, exist_ok=True)
        if force or not _genomad_db_exists(genomad_db):
            info("Downloading genomad database...")
            _run_cmd(["genomad", "download-database", str(genomad_db)])
        else:
            info("genomad database already exists; use --force to re-download")
    else:
        info("Skipping genomad database download (--skip-genomad)")

    emapper_dir.mkdir(parents=True, exist_ok=True)
    emapper_tar = emapper_dir.joinpath("mmseqs.tar.gz")
    mmseqs_folder = emapper_dir.joinpath("mmseqs")

    if skip_emapper:
        info("Skipping mmseqs/emapper DB download (--skip-emapper)")
    elif force or not mmseqs_folder.exists():
        warn("Downloading and extracting mmseqs; folder missing or --force is set")

        ok = _download_url(EMAPPER_URL, emapper_tar, num_threads=num_threads)
        if ok:
            ok_extract = extract_tar(emapper_tar, emapper_dir, threads=num_threads)
            if ok_extract:
                info(f"Extracted {emapper_tar.name} into {emapper_dir}")
            else:
                warn(f"Failed to extract {emapper_tar.name}")
            try:
                emapper_tar.unlink()
                info(f"Removed {emapper_tar} to save space")
            except Exception as e:
                warn(f"Failed to remove {emapper_tar}: {e}")
            padded_prefix = mmseqs_folder.joinpath("mmseqs.db_pad")
            if not padded_prefix.exists():
                try:
                    cmd = ["mmseqs", "makepaddedseqdb", "mmseqs.db", "mmseqs.db_pad"]
                    _run_cmd(cmd, cwd=mmseqs_folder)
                except Exception as e:
                    warn(f"Failed to create padded mmseqs DB: {e}")
        else:
            warn(f"Download failed from {EMAPPER_URL}")
    else:
        info(f"{mmseqs_folder} already exists; skipping download and extraction")

    if skip_parquet:
        info("Skipping eggNOG parquet support files (--skip-parquet)")
    else:
        for url in (EGGNOG_OG, EGGNOG_PROTS):
            dest = emapper_dir.joinpath(Path(url).name)
            if force or not dest.exists():
                _download_url(url, dest, num_threads=num_threads)
            else:
                info(f"{dest} exists; use --force to re-download")

    contig_dir.mkdir(parents=True, exist_ok=True)
    contig_dest = contig_dir.joinpath("contig_lengths.parquet")
    if skip_contig_lengths:
        info("Skipping contig_lengths.parquet download (--skip-contig-lengths)")
    elif force or not contig_dest.exists():
        _download_url(CONTIGS_URL, contig_dest, num_threads=num_threads)
    else:
        info(f"{contig_dest} exists; use --force to re-download")

    stage_done("Databases download complete")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--force", action="store_true")
    args = p.parse_args()
    main(force=args.force)
