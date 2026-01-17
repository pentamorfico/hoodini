from pathlib import Path

import duckdb
import polars as pl

from hoodini.utils.logging_utils import info

hive_dir = Path("/home/klaupaucius/software/hoodini/hoodini/data/contig_lengths")
output_file = Path(
    "/home/klaupaucius/software/hoodini/hoodini/data/contig_lengths/contig_lengths.parquet"
)

# Use DuckDB for memory-efficient parquet combining
con = duckdb.connect(":memory:")
con.execute('SET memory_limit = "4GB"')

parquet_glob = str(hive_dir / "part-*.parquet")
df = con.execute(f"""
    SELECT * FROM read_parquet('{parquet_glob}')
""").pl()
con.close()

df.write_parquet(output_file, compression="zstd")

info(f"Combined parquet written to: {output_file}")
