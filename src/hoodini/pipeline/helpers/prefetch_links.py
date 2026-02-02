#!/usr/bin/env python3
"""Generate prefetched NCBI dataset links for assemblies on the fly.

Provides:
- get_prefetched_link(asm_id, filetype) -> str
- get_prefetched_link_table(asm_id_list, kinds=...) -> pandas.DataFrame
- CLI to produce TSV: assembly_id\tfiletype\tlink

This uses the same encoding approach as the user's reference script to build
datasets API URLs and fetches the sequence-report JSON when a file URL
requires the assembly ftp path.
"""

from __future__ import annotations

import argparse
import base64
import sys
import zlib
from collections.abc import Iterable
from importlib.resources import files

import duckdb
import polars as pl
import requests

from hoodini.utils.logging_utils import warn

BASE = "https://api.ncbi.nlm.nih.gov/datasets/fetch_h"
METHOD_SEQ_B64 = "QXNzZW1ibHlEYXRhc2V0SW50ZXJuYWwuR2V0U2VxdWVuY2VSZXBvcnQ"
METHOD_FILE_B64 = "R2V0UmVtb3RlRGF0YWZpbGU"

FULL_NAME = {
    "gbff": "{acc}_{asm}_genomic.gbff.gz",
    "gff": "{acc}_{asm}_genomic.gff.gz",
    "fna": "{acc}_{asm}_genomic.fna.gz",
    "faa": "{acc}_{asm}_protein.faa.gz",
}


def varint(n: int) -> bytes:
    out = bytearray()
    while True:
        b = n & 0x7F
        n >>= 7
        out.append(b | (0x80 if n else 0))
        if not n:
            return bytes(out)


def b64url_no_pad(b: bytes) -> str:
    return base64.urlsafe_b64encode(b).decode().rstrip("=")


def encode_seqrep_tail(acc: str) -> str:
    s = acc.encode("ascii")
    msg = b"\x0a" + varint(len(s)) + s
    return b64url_no_pad(zlib.compress(msg, 9))


def encode_datafile_tail(ncbi_ftp_url: str) -> str:
    s = ncbi_ftp_url.encode("utf-8")
    inner = b"\x0a" + varint(len(s)) + s
    outer = b"\x1a" + varint(len(inner)) + inner
    msg = outer + b"\x30\x01" + b"\x38\x01"
    return b64url_no_pad(zlib.compress(msg, 9))


def make_seqrep_url(acc: str) -> str:
    return f"{BASE}/{METHOD_SEQ_B64}/{encode_seqrep_tail(acc)}"


def make_file_url(ftp_path_https: str, acc: str, asm: str, kind: str):
    suffix = ftp_path_https.split("https://ftp.ncbi.nlm.nih.gov/", 1)[1].rstrip("/")
    fname = FULL_NAME[kind].format(acc=acc, asm=asm)
    ncbi_ftp_url = f"ncbi+ftp://{suffix}/{fname}"
    tail = encode_datafile_tail(ncbi_ftp_url)
    return f"{BASE}/{METHOD_FILE_B64}/{tail}", fname


def derive_asm_from_ftp_path(ftp_path_https: str, acc: str) -> str:
    last = ftp_path_https.rstrip("/").split("/")[-1]
    return last[len(acc) + 1 :] if last.startswith(acc + "_") else last


def _find_in_json(obj, key: str):
    """Recursively find first occurrence of a key in a nested JSON-like obj."""
    if isinstance(obj, dict):
        if key in obj:
            return obj[key]
        for v in obj.values():
            res = _find_in_json(v, key)
            if res is not None:
                return res
    elif isinstance(obj, list):
        for item in obj:
            res = _find_in_json(item, key)
            if res is not None:
                return res
    return None


_LOCAL_ASM_MAP: dict | None = None


def _load_local_assembly_map() -> dict:
    """Load mapping assembly_accession -> ftp_path from packaged parquet file.

    Uses DuckDB for memory-efficient querying. Returns empty dict if the file
    is missing or cannot be read.
    """
    global _LOCAL_ASM_MAP
    if _LOCAL_ASM_MAP is not None:
        return _LOCAL_ASM_MAP
    try:
        path = files("hoodini").joinpath("data", "assembly_summary.parquet")
        if not path.exists():
            _LOCAL_ASM_MAP = {}
            return _LOCAL_ASM_MAP

        con = duckdb.connect(":memory:")
        con.execute('SET memory_limit = "4GB"')
        df = con.execute(f"""
            SELECT 
                CAST(assembly_accession AS VARCHAR) as assembly_accession,
                CAST(ftp_path AS VARCHAR) as ftp_path
            FROM read_parquet('{str(path)}')
        """).pl()
        con.close()

        _LOCAL_ASM_MAP = dict(
            zip(
                df["assembly_accession"].to_list(),
                df["ftp_path"].to_list(),
            )
        )
        return _LOCAL_ASM_MAP
    except Exception:
        _LOCAL_ASM_MAP = {}
        return _LOCAL_ASM_MAP


def get_prefetched_link(
    asm_id: str,
    filetype: str,
) -> str:
    """Return a prefetch link for a single assembly and filetype.

    filetype may be one of the keys of FULL_NAME (gbff,gff,fna,faa) or
    'sequence_report' (or 'seqrep').
    """
    filetype = filetype.lower()
    if filetype in ("sequence_report", "seqrep"):
        return make_seqrep_url(asm_id)

    if filetype not in FULL_NAME:
        raise ValueError(f"unsupported filetype: {filetype}")

    filetype = filetype.lower()
    if filetype in ("sequence_report", "seqrep"):
        return make_seqrep_url(asm_id)

    if filetype not in FULL_NAME:
        raise ValueError(f"unsupported filetype: {filetype}")

    asm_map = _load_local_assembly_map()
    if asm_id not in asm_map:
        warn(f"Assembly {asm_id} not present in assembly_summary.parquet")
        raise KeyError(f"assembly {asm_id} not present")

    ftp = asm_map.get(asm_id)
    if not ftp or str(ftp).strip().lower() == "na":
        warn(f"Assembly {asm_id} has no ftp_path in assembly_summary.parquet")
        raise ValueError(f"assembly {asm_id} missing ftp_path")

    asm = derive_asm_from_ftp_path(ftp, asm_id)
    url, _ = make_file_url(ftp, asm_id, asm, filetype)
    return url


def get_prefetched_link_table(
    asm_id_list: Iterable[str],
    kinds: list[str] | None = None,
    *,
    seqrep_only: bool = False,
) -> pl.DataFrame:
    """Return a DataFrame with columns: assembly_id, filetype, link

    If ``seqrep_only`` is True, only sequence_report links are generated and
    the local `assembly_summary.parquet` is not read (no ftp_path lookup).

    For each assembly in asm_id_list and each kind in kinds, attempt to
    produce a link. If an individual link cannot be produced it will be
    skipped (no row added)."""
    if kinds is None:
        kinds = ["gbff", "gff", "fna", "faa", "sequence_report"]

    rows = []

    if seqrep_only:
        for asm in asm_id_list:
            asm = asm.strip()
            if not asm:
                continue
            try:
                link = make_seqrep_url(asm)
            except Exception as e:
                warn(f"Skipped {asm} sequence_report: {e}")
                continue
            rows.append({"assembly_id": asm, "filetype": "sequence_report", "url": link})
        return pl.DataFrame(rows)

    requests.Session()
    for asm in asm_id_list:
        asm = asm.strip()
        if not asm:
            continue
        for k in kinds:
            try:
                link = get_prefetched_link(asm, k)
            except Exception as e:
                warn(f"Skipped {asm} {k}: {e}")
                continue
            rows.append({"assembly_id": asm, "filetype": k, "url": link})

    return pl.DataFrame(rows)


def _cli():
    p = argparse.ArgumentParser(description="Generate prefetched links for NCBI assemblies")
    p.add_argument("accessions", nargs="*", help="assembly accessions (GCF_/GCA_...)")
    p.add_argument("-i", "--input", help="file with one accession per line")
    p.add_argument(
        "-k", "--kinds", default="gbff,gff,fna,faa,sequence_report", help="comma-separated kinds"
    )
    p.add_argument("-o", "--output", help="output TSV file (defaults to stdout)")
    args = p.parse_args()

    accs = list(args.accessions or [])
    if args.input:
        with open(args.input) as fh:
            accs.extend([l.strip() for l in fh if l.strip()])

    kinds = [x.strip() for x in args.kinds.split(",") if x.strip()]
    df = get_prefetched_link_table(accs, kinds=kinds)
    if args.output:
        df.write_csv(
            args.output,
            separator="\t",
            include_header=False,
            columns=["assembly_id", "filetype", ""],
            header=True,
        )
    else:
        sys.stdout.write(
            df.write_csv(
                separator="\t",
                include_header=False,
                columns=["assembly_id", "filetype", "url"],
                header=True,
            ).rstrip()
            + "\n"
        )


if __name__ == "__main__":
    _cli()
