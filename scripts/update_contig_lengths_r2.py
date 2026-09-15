#!/usr/bin/env python
"""Merge newly downloaded contig lengths into the published R2 parquet.

This is intentionally *not* part of the hoodini package (same rationale as
update_assembly_summary_r2.py): it's a maintenance script run by hand (or via
s5cmd/r2cmd for the upload step), not something end users need.

Background:
- `hoodini download contig_lengths` only fetches assemblies that are *missing*
  from the local `src/hoodini/data/contig_lengths/*.parquet` dataset **and**
  newer than the R2 published copy's Last-Modified date (see
  `hoodini.download.contig_lengths`). It never re-downloads what's already on
  R2, so the local dataset only ever contains the "delta" since the last R2
  publish.
- The local dataset's schema has grown over time and now has 4 extra columns
  (`gcPercent`, `sequenceName`, `sortOrder`, `unlocalizedCount`) that the
  currently published R2 parquet does not have. To keep the published schema
  stable for existing consumers, this script drops those extra columns when
  merging -- it does NOT attempt to backfill them into old R2 rows.

What it does:
1. Reads the remote contig_lengths.parquet schema/rows via DuckDB httpfs.
2. Reads the local delta parts, keeping only the columns that already exist
   in the remote schema.
3. Anti-joins on assemblyAccession (local wins if somehow overlapping) and
   concatenates, writing a single new parquet.
4. Uploads the result to R2 via s5cmd (the `r2cmd` shell alias), replacing
   contig_lengths.parquet in place.

Usage:
    python scripts/update_contig_lengths_r2.py [--dry-run]

Requires:
- duckdb, s5cmd (invoked via the same AWS_PROFILE=cloudflare-r2 endpoint used
  by the `r2cmd` shell alias -- no boto3 dependency needed here).
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import duckdb

REMOTE_URL = "https://storage.hoodini.bio/contig_lengths.parquet"
BUCKET_KEY = "s3://hoodini/contig_lengths.parquet"
# Fallback endpoint (matches the `r2cmd` shell alias used for local, manual
# runs). In CI, R2_ACCOUNT_ID (a GitHub secret) takes precedence.
DEFAULT_R2_ACCOUNT_ID = "dfc06c27709c4e1fdf99d5dda4b43dc9"

LOCAL_PARTS_GLOB = str(
    Path(__file__).resolve().parent.parent
    / "src"
    / "hoodini"
    / "data"
    / "contig_lengths"
    / "*.parquet"
)

# The published schema on R2 today. We only keep these columns from the local
# delta; any newer columns added to the local fetcher are dropped here so the
# published file's schema doesn't drift underneath existing consumers.
REMOTE_COLUMNS = [
    "assemblyAccession",
    "assemblyUnit",
    "assignedMoleculeLocationType",
    "chrName",
    "gcCount",
    "genbankAccession",
    "length",
    "refseqAccession",
    "role",
]


def _r2_endpoint() -> str:
    account_id = os.environ.get("R2_ACCOUNT_ID", DEFAULT_R2_ACCOUNT_ID)
    return f"https://{account_id}.r2.cloudflarestorage.com"


def _s5cmd_env() -> dict:
    """Build the subprocess env for s5cmd.

    Two supported credential sources:
    - CI (GitHub Actions): R2_ACCESS_KEY_ID / R2_SECRET_ACCESS_KEY secrets ->
      mapped to the standard AWS_ACCESS_KEY_ID / AWS_SECRET_ACCESS_KEY vars.
    - Local/manual runs: no R2_* secrets set, so fall back to the
      `cloudflare-r2` named AWS profile (the same one the `r2cmd` shell alias
      uses), which the operator is expected to have configured already.
    """
    env = dict(os.environ)
    access_key = env.get("R2_ACCESS_KEY_ID")
    secret_key = env.get("R2_SECRET_ACCESS_KEY")
    if access_key and secret_key:
        env["AWS_ACCESS_KEY_ID"] = access_key
        env["AWS_SECRET_ACCESS_KEY"] = secret_key
        env.pop("AWS_PROFILE", None)
    else:
        env["AWS_PROFILE"] = "cloudflare-r2"
    return env


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Build the merged parquet and print stats, but skip the R2 upload.",
    )
    args = parser.parse_args()

    con = duckdb.connect(":memory:")
    con.execute("INSTALL httpfs; LOAD httpfs;")
    con.execute('SET memory_limit = "8GB"')

    # union_by_name=True: local part files have accumulated schema drift over
    # time (some lack gcPercent/sequenceName/sortOrder/unlocalizedCount), so
    # match columns by name instead of position across the glob.
    local_count = con.execute(
        f"SELECT COUNT(DISTINCT assemblyAccession) "
        f"FROM read_parquet('{LOCAL_PARTS_GLOB}', union_by_name=True)"
    ).fetchone()[0]
    remote_count = con.execute(
        f"SELECT COUNT(DISTINCT assemblyAccession) FROM read_parquet('{REMOTE_URL}')"
    ).fetchone()[0]
    print(f"Local delta assemblies: {local_count:,}")
    print(f"Remote published assemblies: {remote_count:,}")

    cols_sql = ", ".join(REMOTE_COLUMNS)

    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir) / "contig_lengths.parquet"

        # Local rows win when an assemblyAccession exists in both (shouldn't
        # normally happen given the fetcher's own missing/date filtering, but
        # be defensive): drop the remote rows for any accession also present
        # locally, then union with the (column-trimmed) local delta.
        con.execute(
            f"""
            COPY (
                SELECT {cols_sql}
                FROM read_parquet('{REMOTE_URL}')
                WHERE assemblyAccession NOT IN (
                    SELECT DISTINCT assemblyAccession
                    FROM read_parquet('{LOCAL_PARTS_GLOB}', union_by_name=True)
                )
                UNION ALL
                SELECT {cols_sql}
                FROM read_parquet('{LOCAL_PARTS_GLOB}', union_by_name=True)
            ) TO '{output_path}' (FORMAT PARQUET, COMPRESSION ZSTD)
            """
        )

        merged_count = con.execute(
            f"SELECT COUNT(*), COUNT(DISTINCT assemblyAccession) FROM read_parquet('{output_path}')"
        ).fetchone()
        size_mb = output_path.stat().st_size / (1024 * 1024)
        print(f"Merged parquet: {merged_count[0]:,} rows, {merged_count[1]:,} assemblies")
        print(f"Merged file size: {size_mb:.1f} MiB at {output_path}")

        if args.dry_run:
            print("--dry-run given; skipping upload.")
            return 0

        print(f"Uploading to {BUCKET_KEY} via s5cmd...")
        result = subprocess.run(
            [
                "s5cmd",
                "--endpoint-url",
                _r2_endpoint(),
                "cp",
                str(output_path),
                BUCKET_KEY,
            ],
            env=_s5cmd_env(),
        )
        if result.returncode != 0:
            print("Upload failed.", file=sys.stderr)
            return 1

    print(f"Published {BUCKET_KEY}.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
