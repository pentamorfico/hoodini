from pathlib import Path

import polars as pl

from hoodini.utils.logging_utils import info

hive_dir = Path("/home/klaupaucius/software/hoodini/hoodini/data/contig_lengths")
output_file = Path(
    "/home/klaupaucius/software/hoodini/hoodini/data/contig_lengths/contig_lengths.parquet"
)

lazy_df = pl.scan_parquet(str(hive_dir / "part-*.parquet"), allow_missing_columns=True)

df = lazy_df.collect()

df.write_parquet(output_file, compression="zstd")

info(f"Combined parquet written to: {output_file}")
