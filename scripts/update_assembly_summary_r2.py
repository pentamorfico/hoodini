#!/usr/bin/env python
from __future__ import annotations

import json
import os
import sys
import tempfile
from pathlib import Path

import requests

MANIFEST_KEY = "assembly_summary_manifest.json"
PARQUET_KEY = "assembly_summary.parquet"
PUBLIC_BASE_URL = "https://storage.hoodini.bio"


def _head_metadata(url: str) -> dict:
    # Force an uncompressed response: NCBI's server omits Content-Length when
    # negotiating gzip (requests' default), which would otherwise leave us
    # with only Last-Modified as a change signal.
    r = requests.head(
        url, timeout=30, allow_redirects=True, headers={"Accept-Encoding": "identity"}
    )
    r.raise_for_status()
    return {
        "etag": r.headers.get("ETag", ""),
        "last_modified": r.headers.get("Last-Modified", ""),
        "content_length": r.headers.get("Content-Length", ""),
    }


def _current_manifest(urls: list[str]) -> dict:
    return {url: _head_metadata(url) for url in urls}


def _published_manifest() -> dict | None:
    r = requests.get(f"{PUBLIC_BASE_URL}/{MANIFEST_KEY}", timeout=30)
    if r.status_code == 404:
        return None
    r.raise_for_status()
    return r.json()


def _r2_client():
    import boto3

    account_id = os.environ["R2_ACCOUNT_ID"]
    return boto3.client(
        "s3",
        endpoint_url=f"https://{account_id}.r2.cloudflarestorage.com",
        aws_access_key_id=os.environ["R2_ACCESS_KEY_ID"],
        aws_secret_access_key=os.environ["R2_SECRET_ACCESS_KEY"],
        region_name="auto",
    )


def main() -> int:
    from hoodini.download.assembly_summary import download_assembly_db, generate_summary_urls

    urls = generate_summary_urls(["refseq", "genbank"], include_historical=True)
    print(f"Checking {len(urls)} NCBI assembly summary URLs for changes...")
    current = _current_manifest(urls)

    published = _published_manifest()
    if published and published.get("sources") == current:
        print("No changes detected upstream; skipping rebuild.")
        return 0

    print("Change detected (or no prior build found); rebuilding merged parquet...")
    columns = [
        "assembly_accession",
        "refseq_category",
        "taxid",
        "species_taxid",
        "organism_name",
        "infraspecific_name",
        "isolate",
        "assembly_level",
        "genome_rep",
        "gbrs_paired_asm",
        "paired_asm_comp",
        "group",
        "ftp_path",
        "genome_size",
        "gc_percent",
        "replicon_count",
        "scaffold_count",
        "contig_count",
        "seq_rel_date",
    ]

    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = Path(tmpdir) / PARQUET_KEY
        download_assembly_db(
            dbs=["refseq", "genbank"],
            output_path=output_path,
            columns_to_keep=columns,
            include_historical=True,
        )

        size_mb = output_path.stat().st_size / (1024 * 1024)
        print(f"Built {output_path} ({size_mb:.1f} MiB); uploading to R2...")

        bucket = os.environ["R2_BUCKET_NAME"]
        client = _r2_client()
        client.upload_file(str(output_path), bucket, PARQUET_KEY)

        manifest = {"sources": current}
        client.put_object(
            Bucket=bucket,
            Key=MANIFEST_KEY,
            Body=json.dumps(manifest, indent=2).encode("utf-8"),
            ContentType="application/json",
        )

    print(f"Published {PARQUET_KEY} and {MANIFEST_KEY} to R2 bucket '{bucket}'.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
