#!/usr/bin/env python
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import tempfile
from glob import glob
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


def _tolerant_select_sql(
    con: duckdb.DuckDBPyConnection,
    source: str,
    columns: list[str],
    expected_types: dict[str, str] | None = None,
) -> tuple[str, list[str]]:
    """Build a SELECT list over `source` that tolerates schema drift.

    The published schema (REMOTE_COLUMNS) is the contract: columns present in
    the source pass through, columns absent are filled with a typed NULL, and
    any extra source columns are dropped. This matters because union_by_name
    only back-fills a column with NULL when at least one part has it — when NO
    part has it, the unified schema lacks the column entirely and a bare SELECT
    fails with a BinderException (NCBI's sequence_report format is actively
    migrating; gcCount vanished from the 2026-09 delta).
    """
    present = {
        row[0]: row[1]
        for row in con.execute(f"DESCRIBE SELECT * FROM {source}").fetchall()
    }
    missing = [c for c in columns if c not in present]
    kinds = expected_types or {}
    select_parts = [
        c if c in present else f"CAST(NULL AS {kinds.get(c, 'VARCHAR')}) AS {c}"
        for c in columns
    ]
    return ", ".join(select_parts), missing


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Build the merged parquet and print stats, but skip the R2 upload.",
    )
    args = parser.parse_args()

    # No local delta parts means `hoodini download contig_lengths` found
    # nothing new/missing since the last R2 publish (e.g. it was already run
    # manually, or a prior cron run already merged everything). Nothing to
    # do in that case -- exit cleanly instead of failing on an empty glob.
    if not glob(LOCAL_PARTS_GLOB):
        print("No local delta parquet parts found; nothing to merge. Skipping.")
        return 0

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

    remote_source = f"read_parquet('{REMOTE_URL}')"
    local_source = f"read_parquet('{LOCAL_PARTS_GLOB}', union_by_name=True)"

    # Types come from the currently published parquet: it is the schema
    # contract, so NULLs injected on either side match its column types.
    remote_types = {
        row[0]: row[1]
        for row in con.execute(f"DESCRIBE SELECT * FROM {remote_source}").fetchall()
    }

    remote_sql, remote_missing = _tolerant_select_sql(con, remote_source, REMOTE_COLUMNS, remote_types)
    if remote_missing:
        print(f"⚠️  Parquet publicado sin columnas {remote_missing}; rellenando con NULL.")
    local_sql, local_missing = _tolerant_select_sql(con, local_source, REMOTE_COLUMNS, remote_types)
    if local_missing:
        print(
            f"⚠️  Delta local sin columnas {local_missing} — ¿formato sequence_report de "
            f"NCBI cambiado? Rellenando con NULL (sin dato para esos assemblies)."
        )

    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir) / "contig_lengths.parquet"

        # Local rows win when an assemblyAccession exists in both (shouldn't
        # normally happen given the fetcher's own missing/date filtering, but
        # be defensive): drop the remote rows for any accession also present
        # locally, then union with the (column-trimmed) local delta. Both
        # selects tolerate missing columns on either side.
        con.execute(
            f"""
            COPY (
                SELECT {remote_sql}
                FROM {remote_source}
                WHERE assemblyAccession NOT IN (
                    SELECT DISTINCT assemblyAccession
                    FROM {local_source}
                )
                UNION ALL
                SELECT {local_sql}
                FROM {local_source}
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
