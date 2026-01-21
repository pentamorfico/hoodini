import asyncio
import concurrent.futures
import csv
import io
import logging
from pathlib import Path

import httpx
import polars as pl
import requests
from rich.live import Live
from rich.logging import RichHandler
from rich.panel import Panel
from rich.progress import (
    BarColumn,
    DownloadColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeRemainingColumn,
    TransferSpeedColumn,
)
from rich.table import Table
from rich.text import Text

from hoodini.utils.logging_utils import console

BACDIVE_URL = "https://bacdive.dsmz.de/advsearch/csv?fg%5B0%5D%5Bgc%5D=OR&fg%5B0%5D%5Bfl%5D%5B1%5D%5Bfd%5D=INSDC+accession&fg%5B0%5D%5Bfl%5D%5B1%5D%5Bfo%5D=contains&fg%5B0%5D%5Bfl%5D%5B1%5D%5Bfv%5D=%2A&fg%5B0%5D%5Bfl%5D%5B1%5D%5Bfvd%5D=sequence_genomes-insdc_acc-7"
PHAGEDIVE_URL = "https://phagedive.dsmz.de/advsearch/csv?fg%5B0%5D%5Bgc%5D=OR&fg%5B0%5D%5Bfl%5D%5B1%5D%5Bfd%5D=Assembly+accession+number&fg%5B0%5D%5Bfl%5D%5B1%5D%5Bfo%5D=contains&fg%5B0%5D%5Bfl%5D%5B1%5D%5Bfv%5D=%2A&fg%5B0%5D%5Bfl%5D%5B1%5D%5Bfvd%5D=sequence_genome-assembly_accession_number-10"

DATA_DIR = Path(__file__).parent.parent / "data"
BACDIVE_OUT = DATA_DIR / "bacdive_download.csv"
PHAGEDIVE_OUT = DATA_DIR / "phagedive_download.csv"

logger = logging.getLogger("type_dive")
logger.setLevel(logging.INFO)
handler = RichHandler(console=console, show_time=True, show_level=True, show_path=False)
logger.addHandler(handler)


def download_csv(url, out_path, desc):
    logger.info(f"Downloading {desc} from {url}")
    try:
        r = requests.get(url, stream=True, timeout=60)
        r.raise_for_status()
        with open(out_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
        logger.info(f"Downloaded {desc} to {out_path}")
    except Exception as e:
        logger.error(f"Failed to download {desc}: {e}")
        raise


def fill_empty_with_previous(rows):
    if not rows:
        return rows
    prev = rows[0][:]
    filled = [prev]
    for row in rows[1:]:
        new_row = [cell if cell.strip() else prev[i] for i, cell in enumerate(row)]
        filled.append(new_row)
        prev = new_row
    return filled


def parse_bacdive_csv(in_path, out_path):
    logger.info(f"Parsing BacDive CSV: {in_path}")
    with open(in_path, newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        for row in reader:
            if row and row[0].strip() == "ID":
                header = row
                break
        else:
            logger.error("Could not find BacDive header row!")
            return
        data = [row for row in reader if row and any(cell.strip() for cell in row)]
    data_filled = fill_empty_with_previous(data)
    idx_id = header.index("ID")
    idx_strain = header.index("strain_number_header")
    idx_assembly = header.index("INSDC accession")
    out_header = ["bacdive_id", "collection_id", "assembly_id"]
    out_rows = []
    for row in data_filled:
        if len(row) < max(idx_id, idx_strain, idx_assembly) + 1:
            row = row + [""] * (max(idx_id, idx_strain, idx_assembly) + 1 - len(row))
        out_rows.append([row[idx_id], row[idx_strain], row[idx_assembly]])
    return out_header, out_rows


def parse_phagedive_csv(in_path, out_path):
    logger.info(f"Parsing PhageDive CSV: {in_path}")
    with open(in_path, newline="", encoding="utf-8") as f:
        reader = csv.reader(f)
        for row in reader:
            if row and row[0].strip() == "ID":
                header = row
                break
        else:
            logger.error("Could not find PhageDive header row!")
            return
        data = [row for row in reader if row and any(cell.strip() for cell in row)]
    idx_phagedive_id = header.index("ID")
    idx_coll = header.index("Collection number")
    idx_assembly = header.index("Assembly accession number")
    out_header = ["phagedive_id", "collection_id", "assembly_id"]
    out_rows = []
    for row in data:
        if len(row) < max(idx_phagedive_id, idx_coll, idx_assembly) + 1:
            row = row + [""] * (max(idx_phagedive_id, idx_coll, idx_assembly) + 1 - len(row))
        out_rows.append([row[idx_phagedive_id], row[idx_coll], row[idx_assembly]])
    return out_header, out_rows


def print_table(header, rows, title):
    table = Table(title=title, show_lines=True)
    for col in header:
        table.add_column(col, style="cyan")
    for row in rows[:10]:
        table.add_row(*row)
    console.print(table)


def main():
    DATA_DIR.mkdir(exist_ok=True)

    progress_columns = [
        TextColumn("[bold blue]{task.description}", justify="left"),
        SpinnerColumn(),
        BarColumn(bar_width=None),
        DownloadColumn(),
        TransferSpeedColumn(),
        TimeRemainingColumn(),
    ]
    progress = Progress(*progress_columns, transient=False)

    log_text = Text(justify="left")
    log_panel = Panel(log_text, title="Logs", border_style="dim")

    class RichTextHandler(logging.Handler):
        def __init__(self, rich_text):
            super().__init__()
            self.rich_text = rich_text
            self.setFormatter(
                logging.Formatter("[%(levelname)s] %(asctime)s | %(message)s", datefmt="%H:%M:%S")
            )

        def emit(self, record):
            msg = self.format(record)
            style = {
                "INFO": "white",
                "WARNING": "yellow",
                "ERROR": "red",
                "CRITICAL": "bold red",
                "DEBUG": "dim",
            }.get(record.levelname, "white")
            self.rich_text.append(msg + "\n", style=style)

    for h in logger.handlers[:]:
        logger.removeHandler(h)
    rich_handler = RichTextHandler(log_text)
    logger.addHandler(rich_handler)
    logger.propagate = False

    for noisy in ["httpx", "urllib3", "requests"]:
        logging.getLogger(noisy).setLevel(logging.WARNING)

    class ProgressAndLogs:
        def __rich_console__(self, console, options):
            yield progress
            yield log_panel

    async def get_content_length(url):
        try:
            async with httpx.AsyncClient(timeout=30) as client:
                r = await client.head(url, follow_redirects=True)
                if r.status_code == 200 and "content-length" in r.headers:
                    return int(r.headers["content-length"])
        except Exception:
            pass
        return None

    async def async_download():
        loop = asyncio.get_event_loop()

        bacdive_size = 2 * 1024 * 1024
        phagedive_size = 30 * 1024
        t1 = progress.add_task("BacDive", total=bacdive_size)
        t2 = progress.add_task("PhageDive", total=phagedive_size)

        def download_csv_in_memory(url, desc, task_id):
            logger.info(f"Downloading {desc} ...")
            try:
                with requests.get(url, stream=True, timeout=60) as r:
                    r.raise_for_status()
                    buf = io.StringIO()
                    decoder = r.iter_lines(decode_unicode=True)
                    total = 0
                    for line in decoder:
                        buf.write(line + "\n")
                        total += len(line.encode("utf-8")) + 1
                        progress.update(task_id, advance=len(line.encode("utf-8")) + 1)
                    buf.seek(0)
                logger.info(f"Downloaded {desc} in memory")
                return buf
            except Exception as e:
                logger.error(f"Failed to download {desc}: {e}")
                raise

        with concurrent.futures.ThreadPoolExecutor() as executor:
            fut1 = loop.run_in_executor(
                executor, download_csv_in_memory, BACDIVE_URL, "BacDive bacteria DB", t1
            )
            fut2 = loop.run_in_executor(
                executor, download_csv_in_memory, PHAGEDIVE_URL, "PhageDive DB", t2
            )
            bacdive_buf, phagedive_buf = await asyncio.gather(fut1, fut2)
            progress.update(t1, completed=bacdive_size)
            progress.update(t2, completed=phagedive_size)
            logger.info("Parsing BacDive CSV...")
            bacdive_buf.seek(0)
            reader = csv.reader(bacdive_buf)
            for row in reader:
                if row and row[0].strip() == "ID":
                    b_header = row
                    break
            else:
                logger.error("Could not find BacDive header row!")
                return
            b_data = [row for row in reader if row and any(cell.strip() for cell in row)]
            b_data_filled = fill_empty_with_previous(b_data)
            idx_id = b_header.index("ID")
            idx_strain = b_header.index("strain_number_header")
            idx_assembly = b_header.index("INSDC accession")
            b_rows = []
            for row in b_data_filled:
                if len(row) < max(idx_id, idx_strain, idx_assembly) + 1:
                    row = row + [""] * (max(idx_id, idx_strain, idx_assembly) + 1 - len(row))
                b_rows.append([row[idx_id], row[idx_strain], row[idx_assembly]])
            b_header_out = ["bacdive_id", "collection_id", "assembly_id"]
            df_bacdive = pl.DataFrame(
                {k: [row[i] for row in b_rows] for i, k in enumerate(b_header_out)}
            )
            df_bacdive = df_bacdive.rename({"bacdive_id": "dive_id"})
            df_bacdive = df_bacdive.with_columns([pl.lit("bacteria").alias("dive_type")])
            logger.info("Parsing PhageDive CSV...")
            phagedive_buf.seek(0)
            reader = csv.reader(phagedive_buf)
            for row in reader:
                if row and row[0].strip() == "ID":
                    p_header = row
                    break
            else:
                logger.error("Could not find PhageDive header row!")
                return
            p_data = [row for row in reader if row and any(cell.strip() for cell in row)]
            idx_phagedive_id = p_header.index("ID")
            idx_coll = p_header.index("Collection number")
            idx_assembly = p_header.index("Assembly accession number")
            p_rows = []
            for row in p_data:
                if len(row) < max(idx_phagedive_id, idx_coll, idx_assembly) + 1:
                    row = row + [""] * (
                        max(idx_phagedive_id, idx_coll, idx_assembly) + 1 - len(row)
                    )
                p_rows.append([row[idx_phagedive_id], row[idx_coll], row[idx_assembly]])
            p_header_out = ["phagedive_id", "collection_id", "assembly_id"]
            df_phagedive = pl.DataFrame(
                {k: [row[i] for row in p_rows] for i, k in enumerate(p_header_out)}
            )
            df_phagedive = df_phagedive.rename({"phagedive_id": "dive_id"})
            df_phagedive = df_phagedive.with_columns([pl.lit("phage").alias("dive_type")])
            df_combined = pl.concat([df_bacdive, df_phagedive], how="vertical_relaxed")

            # Duplicate GCA_ rows with GCF_ equivalents (RefSeq mirrors GenBank)
            gca_rows = df_combined.filter(pl.col("assembly_id").str.starts_with("GCA_"))
            gcf_rows = gca_rows.with_columns(
                pl.col("assembly_id").str.replace("^GCA_", "GCF_").alias("assembly_id")
            )
            df_combined = pl.concat([df_combined, gcf_rows], how="vertical_relaxed")
            logger.info(f"Added {gcf_rows.height} GCF_ mirror rows for GCA_ assemblies")

            df_combined.write_parquet(DATA_DIR / "dive_combined.parquet")
            logger.info(
                f"Saved combined data to {DATA_DIR / 'dive_combined.parquet'} ({df_combined.height} rows)"
            )
            logger.info("DSMZ BacDive/PhageDive download and normalization complete!")

    with Live(ProgressAndLogs(), console=console, refresh_per_second=10):
        asyncio.run(async_download())


if __name__ == "__main__":
    main()
