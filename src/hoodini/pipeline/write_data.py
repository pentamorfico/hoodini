import base64
import gzip
from datetime import datetime
from importlib.resources import files
from pathlib import Path

import polars as pl
from jinja2 import Environment

from hoodini.utils.polars_adapters import to_polars


def _append_extra_tool_gffs(
    all_gff: pl.DataFrame,
    blast_data: pl.DataFrame | None = None,
    crispr_df: pl.DataFrame | None = None,
    ncrna_data: pl.DataFrame | None = None,
    genomad_df: pl.DataFrame | None = None,
) -> pl.DataFrame:
    """
    Internal helper to transform and concatenate extra tool outputs to GFF format.

    Parameters:
    -----------
    all_gff: Base GFF DataFrame from assembly parsing
    blast_data: Raw BLAST results
    crispr_df: Raw CCTyper CRISPR array predictions
    ncrna_data: Raw ncRNA predictions from cmsearch
    genomad_df: Raw geNomad MGE predictions

    Returns:
    --------
    pl.DataFrame: Extended GFF with all tool results concatenated
    """
    result_gff = all_gff.clone()

    # BLAST regions
    if blast_data is not None and blast_data.height > 0:
        gff_df = pl.DataFrame(
            {
                "seqid": blast_data["seqid"],
                "source": "hoodini",
                "type": "region",
                "start": blast_data["start"],
                "end": blast_data["end"],
                "score": ".",
                "strand": "+",
                "phase": ".",
                "attributes": "ID=" + blast_data["nc_feature"] + ";",
            }
        )
        result_gff = pl.concat([result_gff, gff_df], how="vertical")

    # CCTyper CRISPR arrays
    if crispr_df is not None and crispr_df.height > 0:
        gff_df = crispr_df.select(
            [
                pl.col("seqid"),
                pl.lit("hoodini").alias("source"),
                pl.lit("region").alias("type"),
                pl.col("start"),
                pl.col("end"),
                pl.lit(".").alias("score"),
                pl.lit(".").alias("strand"),
                pl.lit(".").alias("phase"),
                (pl.lit("ID=") + pl.col("nc_feature") + pl.lit(";")).alias("attributes"),
            ]
        )
        result_gff = pl.concat([result_gff, gff_df], how="vertical")

    # ncRNA predictions
    if ncrna_data is not None and ncrna_data.height > 0:
        gff_df = ncrna_data.select(
            [
                pl.col("nucid").alias("seqid"),
                pl.lit("hoodini").alias("source"),
                pl.lit("ncRNA").alias("type"),
                pl.min_horizontal([pl.col("start"), pl.col("end")]).alias("start"),
                pl.max_horizontal([pl.col("start"), pl.col("end")]).alias("end"),
                pl.lit(".").alias("score"),
                pl.col("strand_ncrna").alias("strand"),
                pl.lit(".").alias("phase"),
                (pl.lit("ID=") + pl.col("nc_feature") + pl.lit(";")).alias("attributes"),
            ]
        )
        result_gff = pl.concat([result_gff, gff_df], how="vertical")

    # geNomad MGE regions
    if genomad_df is not None and genomad_df.height > 0:
        gff_df = genomad_df.select(
            [
                pl.col("seqid"),
                pl.lit("hoodini").alias("source"),
                pl.lit("region").alias("type"),
                pl.min_horizontal([pl.col("start"), pl.col("end")]).alias("start"),
                pl.max_horizontal([pl.col("start"), pl.col("end")]).alias("end"),
                pl.lit(".").alias("score"),
                pl.lit(".Z").alias("strand"),
                pl.lit(".").alias("phase"),
                (pl.lit("ID=") + pl.col("mge_type") + pl.lit(";")).alias("attributes"),
            ]
        )
        result_gff = pl.concat([result_gff, gff_df], how="vertical")

    return result_gff


def write_viz_outputs(
    *,
    output_dir: str,
    all_gff: pl.DataFrame,
    all_neigh: pl.DataFrame,
    all_prots: pl.DataFrame,
    den_data: pl.DataFrame,
    tree_str: str,
    records: pl.DataFrame | None = None,
    valid_uids: list[str] | None = None,
    nt_links: pl.DataFrame | None = None,
    pairwise_aa: pl.DataFrame | None = None,
    domains_data: pl.DataFrame | None = None,
    ncrna_data: pl.DataFrame | None = None,
    write_domains: bool = False,
    # Raw extra tool outputs (for GFF transformation)
    blast_data: pl.DataFrame | None = None,
    crispr_df: pl.DataFrame | None = None,
    genomad_df: pl.DataFrame | None = None,
) -> Path:
    """
    Write hoodini visualization-ready files into a hoodini-viz folder.

    Expected Files:
    ---------------
    - all_gff: DataFrame with genomic features
    - all_neigh: DataFrame with neighborhood metadata
    - all_prots: DataFrame with protein sequences and annotations
    - den_data: DataFrame with taxonomy and tree metadata
    - tree_str: Newick format tree string
    - domains_data: Optional DataFrame with domain annotations
    - ncrna_data: Optional DataFrame with ncRNA metadata
    - nt_links: Optional DataFrame with nucleotide alignment links
    - pairwise_aa: Optional DataFrame with protein similarity links

    Generated Files:
    ----------------
    Parquet files (for programmatic access):
        - {output}/hoodini-viz/parquet/gff.parquet
        - {output}/hoodini-viz/parquet/hoods.parquet
        - {output}/hoodini-viz/parquet/protein_metadata.parquet
        - {output}/hoodini-viz/parquet/tree_metadata.parquet
        - {output}/hoodini-viz/parquet/nucleotide_links.parquet
        - {output}/hoodini-viz/parquet/protein_links.parquet
        - {output}/hoodini-viz/parquet/domains.parquet (if domains requested)
        - {output}/hoodini-viz/parquet/domains_metadata.parquet (if domains requested)
        - {output}/hoodini-viz/parquet/ncrna_metadata.parquet (if ncrna requested)

    TSV files (human-readable):
        - {output}/hoodini-viz/tsv/gff.gff
        - {output}/hoodini-viz/tsv/hoods.txt
        - {output}/hoodini-viz/tsv/protein_metadata.txt
        - {output}/hoodini-viz/tsv/tree_metadata.txt
        - {output}/hoodini-viz/tsv/nucleotide_links.txt
        - {output}/hoodini-viz/tsv/protein_links.txt
        - {output}/hoodini-viz/tsv/domains.txt (if domains)
        - {output}/hoodini-viz/tsv/domains_metadata.txt (if domains)
        - {output}/hoodini-viz/tsv/ncrna_metadata.txt (if ncrna)

    Other files:
        - {output}/hoodini-viz/tree.nwk: Newick phylogenetic tree
        - {output}/hoodini-viz/hoodini-viz.html: Standalone interactive viewer
            (embeds all parquet data as base64-encoded strings for portability)

    Parameters (Raw Tool Outputs):
    ------------------------------
    blast_data: Optional raw BLAST DataFrame (auto-transformed to GFF)
    crispr_df: Optional raw CCTyper CRISPR array DataFrame (auto-transformed to GFF)
    genomad_df: Optional raw geNomad MGE DataFrame (auto-transformed to GFF)
    ncrna_data: Optional ncRNA metadata (exported separately + auto-transformed to GFF)

    Returns:
    --------
    Path: The hoodini-viz output directory path
    """
    # Transform and concatenate extra tool outputs to GFF
    all_gff = _append_extra_tool_gffs(
        all_gff=all_gff,
        blast_data=blast_data,
        crispr_df=crispr_df,
        ncrna_data=ncrna_data,
        genomad_df=genomad_df,
    )

    outdir = Path(output_dir) / "hoodini-viz"
    outdir.mkdir(parents=True, exist_ok=True)
    parquet_dir = outdir / "parquet"
    tsv_dir = outdir / "tsv"
    parquet_dir.mkdir(parents=True, exist_ok=True)
    tsv_dir.mkdir(parents=True, exist_ok=True)

    all_gff = to_polars(all_gff) if all_gff is not None else pl.DataFrame()
    all_neigh = to_polars(all_neigh) if all_neigh is not None else pl.DataFrame()
    all_prots = to_polars(all_prots) if all_prots is not None else pl.DataFrame()
    den_data = to_polars(den_data) if den_data is not None else pl.DataFrame()
    nt_links = to_polars(nt_links) if nt_links is not None else None
    pairwise_aa = to_polars(pairwise_aa) if pairwise_aa is not None else None
    domains_data = to_polars(domains_data) if domains_data is not None else None
    ncrna_data = to_polars(ncrna_data) if ncrna_data is not None else None

    # Filter out failed records - only include neighborhoods with valid_uids
    # This ensures failed records don't appear in visualization (gff, hoods, prots)
    if valid_uids is not None:
        valid_uid_set = {str(uid) for uid in valid_uids}
        # Filter all_neigh by unique_id
        if all_neigh.height > 0 and "unique_id" in all_neigh.columns:
            all_neigh = all_neigh.filter(pl.col("unique_id").cast(pl.Utf8).is_in(valid_uid_set))
        # Filter all_prots by unique_id
        if all_prots.height > 0 and "unique_id" in all_prots.columns:
            all_prots = all_prots.filter(pl.col("unique_id").cast(pl.Utf8).is_in(valid_uid_set))
        # Filter all_gff by seqid (matching valid neighborhoods' seqids)
        if all_gff.height > 0 and all_neigh.height > 0 and "seqid" in all_gff.columns:
            valid_seqids = set(all_neigh["seqid"].unique().to_list())
            all_gff = all_gff.filter(pl.col("seqid").is_in(valid_seqids))

    if all_gff is not None and all_gff.height > 0:
        gff_df = all_gff.clone()

        def normalize_attributes(val):
            if val is None:
                return ""
            if isinstance(val, str):
                return val
            if isinstance(val, dict):
                # Flatten dict to "key=value" pairs; preserve ordering best-effort.
                parts = []
                for k, v in val.items():
                    if v is None:
                        continue
                    parts.append(f"{k}={v}")
                return ";".join(parts)
            # Fallback to simple string conversion (avoid double-JSON encoding)
            return str(val)

        if "attributes" in gff_df.columns:
            gff_df = gff_df.with_columns(
                pl.col("attributes").map_elements(normalize_attributes).cast(pl.Utf8, strict=False)
            )
        # GFF text under tsv/, parquet under parquet/
        gff_df.write_csv(tsv_dir / "gff.gff", include_header=False, separator="\t")
        gff_df.write_parquet(parquet_dir / "gff.parquet")
    else:
        empty = pl.DataFrame()
        empty.write_csv(tsv_dir / "gff.gff", include_header=False, separator="\t")
        empty.write_parquet(parquet_dir / "gff.parquet")

    # Base hood columns
    base_hood_cols = [
        "hood_id",
        "seqid",
        "start",
        "end",
        "align_gene",
        "align_start",
        "align_end",
        "align_strand",
    ]

    if all_neigh is not None and all_neigh.height > 0:
        neigh = all_neigh.clone()
        mapping = {
            "unique_id": "hood_id",
            "start_win": "start",
            "end_win": "end",
            "target_prot": "align_gene",
        }
        present_map = {k: v for k, v in mapping.items() if k in neigh.columns}
        if present_map:
            neigh = neigh.rename(present_map)

        # Join with records to get alignment inputs (start/end/strand) if provided
        if records is not None and records.height > 0 and "hood_id" in neigh.columns:
            align_map = {}
            if "start" in records.columns:
                align_map["start"] = "align_start"
            if "end" in records.columns:
                align_map["end"] = "align_end"
            if "strand" in records.columns:
                align_map["strand"] = "align_strand"
            if align_map:
                records_align = records.select(["unique_id"] + list(align_map.keys())).rename(
                    align_map
                )
                records_align = records_align.with_columns(pl.col("unique_id").cast(pl.Utf8))
                neigh = neigh.with_columns(pl.col("hood_id").cast(pl.Utf8))
                neigh = neigh.join(
                    records_align, left_on="hood_id", right_on="unique_id", how="left"
                )

        # Join with records to get user extra columns
        user_extra_cols = []
        if records is not None and records.height > 0:
            from hoodini.utils.validation import RESERVED_COLUMNS

            # Find user extra columns in records (not reserved, not starting with _)
            user_extra_cols = [
                c for c in records.columns if c not in RESERVED_COLUMNS and not c.startswith("_")
            ]
            if user_extra_cols and "hood_id" in neigh.columns:
                # Join on unique_id (records) = hood_id (neigh after rename)
                records_extra = records.select(["unique_id"] + user_extra_cols)
                records_extra = records_extra.with_columns(pl.col("unique_id").cast(pl.Utf8))
                neigh = neigh.with_columns(pl.col("hood_id").cast(pl.Utf8))
                neigh = neigh.join(
                    records_extra, left_on="hood_id", right_on="unique_id", how="left"
                )

        # Select base columns + user extra columns only
        cols = [c for c in base_hood_cols if c in neigh.columns]
        cols.extend([c for c in user_extra_cols if c in neigh.columns])

        if cols:
            hoods_headers = "\t".join(cols) + "\n"
            csv_data = neigh.select(cols).write_csv(separator="\t", include_header=False)
            (tsv_dir / "hoods.txt").write_text(hoods_headers + csv_data, encoding="utf-8")
            neigh.select(cols).write_parquet(parquet_dir / "hoods.parquet")
        else:
            hoods_headers = "\t".join(base_hood_cols) + "\n"
            (tsv_dir / "hoods.txt").write_text(hoods_headers, encoding="utf-8")
            pl.DataFrame().write_parquet(parquet_dir / "hoods.parquet")
    else:
        hoods_headers = "\t".join(base_hood_cols) + "\n"
        (tsv_dir / "hoods.txt").write_text(hoods_headers, encoding="utf-8")
        pl.DataFrame().write_parquet(parquet_dir / "hoods.parquet")

    base_headers = [
        "id",
        "sequence",
        "product",
        "target_prot",
        "target_nuc",
        "unique_id",
        "cluster",
    ]
    if all_prots is not None and all_prots.height > 0:
        prots = all_prots.clone()
        if "fam_cluster" in prots.columns:
            prots = prots.rename({"fam_cluster": "cluster"})

        # Drop redundant aliases if present
        drop_cols = []
        if "protein_id" in prots.columns:
            drop_cols.append("protein_id")
        if "clusterID" in prots.columns:
            drop_cols.append("clusterID")
        if drop_cols:
            prots = prots.drop(drop_cols)

        for col in ["target_prot", "target_nuc", "unique_id"]:
            if col not in prots.columns:
                prots = prots.with_columns(pl.lit("").alias(col))

        if "cluster" not in prots.columns:
            prots = prots.with_columns(pl.lit(None).alias("cluster"))
        else:
            prots = prots.with_columns(pl.col("cluster").round().cast(pl.Int64))

        base_cols = [
            "id",
            "sequence",
            "product",
            "target_prot",
            "target_nuc",
            "unique_id",
            "cluster",
        ]
        base_cols_present = [c for c in base_cols if c in prots.columns]
        extra_cols = [c for c in prots.columns if c not in base_cols]
        prots = prots.select(base_cols_present + extra_cols)

        for col in prots.columns:
            col_dtype = prots.schema.get(col)
            if col_dtype in (pl.Utf8, pl.String):
                prots = prots.with_columns(
                    pl.col(col)
                    .cast(pl.Utf8)
                    .str.replace_all("\r\n", " ")
                    .str.replace_all("\n", " ")
                    .str.replace_all("\r", " ")
                    .str.replace_all("`", "'")
                    .str.replace_all('"', "'")
                    .alias(col)
                )
            elif col_dtype in (pl.Float64, pl.Float32):
                prots = prots.with_columns(pl.col(col).round(2).alias(col))
        csv_data = prots.write_csv(separator="\t", include_header=False)
        protein_headers = "\t".join(prots.columns) + "\n"
        (tsv_dir / "protein_metadata.txt").write_text(protein_headers + csv_data, encoding="utf-8")
        prots.write_parquet(parquet_dir / "protein_metadata.parquet")
    else:
        protein_headers = "\t".join(base_headers) + "\n"
        (tsv_dir / "protein_metadata.txt").write_text(protein_headers, encoding="utf-8")
        pl.DataFrame(schema=dict.fromkeys(base_headers, pl.Utf8)).write_parquet(
            parquet_dir / "protein_metadata.parquet"
        )

    base_tree_cols = [
        "leaf_id",
        "og_index",
        "superkingdom",
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
        "start_win",
        "end_win",
        "strand_win",
        "start_target",
        "end_target",
    ]
    if den_data is not None and den_data.height > 0:
        tree_meta = den_data.clone()
        if "unique_id" in tree_meta.columns:
            tree_meta = tree_meta.rename({"unique_id": "leaf_id"})
        # Include all columns (base + any extras like dive_id, collection_id, dive_type)
        tree_headers = "\t".join(tree_meta.columns) + "\n"
        csv_data = tree_meta.write_csv(separator="\t", include_header=False)
        (tsv_dir / "tree_metadata.txt").write_text(tree_headers + csv_data, encoding="utf-8")
        tree_meta.write_parquet(parquet_dir / "tree_metadata.parquet")
    else:
        tree_headers = "\t".join(base_tree_cols) + "\n"
        (tsv_dir / "tree_metadata.txt").write_text(tree_headers, encoding="utf-8")
        pl.DataFrame().write_parquet(parquet_dir / "tree_metadata.parquet")

    (outdir / "tree.nwk").write_text(tree_str or "", encoding="utf-8")

    if write_domains and domains_data is not None and domains_data.height > 0:
        # Normalize to the exact schema expected by hoodini-viz:
        # gene_id, domainName, start, end, source, evalue, coverage
        df = domains_data.clone()

        def c(name: str):
            return pl.col(name) if name in df.columns else pl.lit(None)

        df = df.with_columns(
            [
                pl.coalesce([c("gene_id"), c("protein_id")]).alias("gene_id_norm"),
                pl.coalesce([c("domainName"), c("domain_id")]).alias("domain_id_norm"),
                pl.coalesce([c("source"), c("database")]).alias("source_norm"),
                pl.coalesce([c("evalue"), c("e_value")]).alias("evalue_norm"),
                pl.coalesce([c("coverage"), c("cov")]).alias("coverage_norm"),
            ]
        )

        # Cast numeric fields; keep nulls when casting fails so we can drop incomplete rows.
        df = df.with_columns(
            [
                c("start").cast(pl.Float64, strict=False).alias("start"),
                c("end").cast(pl.Float64, strict=False).alias("end"),
                pl.col("evalue_norm").cast(pl.Float64, strict=False),
                pl.col("coverage_norm").cast(pl.Float64, strict=False),
            ]
        )

        df_domains = (
            df.select(
                [
                    pl.col("gene_id_norm").alias("gene_id"),
                    pl.col("domain_id_norm").alias("domainName"),
                    pl.col("start"),
                    pl.col("end"),
                    pl.col("source_norm").alias("source"),
                    pl.col("evalue_norm").alias("evalue"),
                    pl.col("coverage_norm").alias("coverage"),
                ]
            )
            .drop_nulls(["gene_id", "domainName", "start", "end"])
            .with_columns(
                [
                    pl.col("start").round(2),
                    pl.col("end").round(2),
                    pl.col("evalue").map_elements(
                        lambda x: f"{x:.2e}" if x is not None else None,
                        return_dtype=pl.Utf8,
                    ),
                    pl.col("coverage").round(2),
                ]
            )
        )

        if df_domains.height == 0:
            df_domains = pl.DataFrame(
                {
                    "gene_id": [],
                    "domainName": [],
                    "start": [],
                    "end": [],
                    "source": [],
                    "evalue": [],
                    "coverage": [],
                }
            )

        df_domains.select(
            ["gene_id", "domainName", "start", "end", "source", "evalue", "coverage"]
        ).write_csv(
            tsv_dir / "domains.txt",
            separator="\t",
            include_header=False,
        )
        df_domains.write_parquet(parquet_dir / "domains.parquet")

        # Build a compact metadata table keyed by domain_id (domainName) if extra columns exist.
        metadata_exclude = {
            "gene_id",
            "protein_id",
            "domain_id",
            "domainName",
            "start",
            "end",
            "database",
            "source",
            "e_value",
            "evalue",
            "cov",
            "coverage",
            "gene_id_norm",
            "domain_id_norm",
            "source_norm",
            "evalue_norm",
            "coverage_norm",
        }
        drop_meta = {"bit_score*alignment_length", "domain_id_base", "domain_id_clean"}
        meta_candidates = [
            c for c in df.columns if c not in metadata_exclude and c not in drop_meta
        ]
        if meta_candidates:
            df_meta = df.select(
                [pl.col("domain_id_norm").alias("domain_id")] + [pl.col(c) for c in meta_candidates]
            ).unique(subset=["domain_id"])
            numeric_cols = [
                c for c, dt in df_meta.schema.items() if getattr(dt, "is_numeric", lambda: False)()
            ]
            df_meta = df_meta.with_columns([pl.col(c).round(2).alias(c) for c in numeric_cols])
            df_meta.write_csv(tsv_dir / "domains_metadata.txt", separator="\t", include_header=True)
            df_meta.write_parquet(parquet_dir / "domains_metadata.parquet")
        else:
            pl.DataFrame().write_csv(
                tsv_dir / "domains_metadata.txt", separator="\t", include_header=True
            )
            pl.DataFrame().write_parquet(parquet_dir / "domains_metadata.parquet")
    else:
        # Always emit empty parquet/text files so the front-end receives a valid data URL.
        empty_domains = pl.DataFrame(
            {
                "protein_id": [],
                "domain_id": [],
                "start": [],
                "end": [],
                "database": [],
                "e_value": [],
                "cov": [],
            }
        )
        empty_domains.write_csv(tsv_dir / "domains.txt", separator="\t", include_header=False)
        empty_domains.write_parquet(parquet_dir / "domains.parquet")

        empty_meta = pl.DataFrame()
        empty_meta.write_csv(tsv_dir / "domains_metadata.txt", separator="\t", include_header=True)
        empty_meta.write_parquet(parquet_dir / "domains_metadata.parquet")

    # Write ncRNA metadata with extended schema: seqid, start, end, type, sequence, structure, score, evalue, accession, strand
    ncrna_cols = [
        "seqid",
        "start",
        "end",
        "type",
        "sequence",
        "structure",
        "score",
        "evalue",
        "accession",
        "strand",
    ]
    if ncrna_data is not None and isinstance(ncrna_data, pl.DataFrame) and ncrna_data.height > 0:
        # Normalize/mapping from cmsearch outputs
        ncrna_meta = ncrna_data
        mappings = []
        if "nucid" in ncrna_meta.columns:
            mappings.append(pl.col("nucid").alias("seqid"))
        if "nc_feature" in ncrna_meta.columns:
            mappings.append(pl.col("nc_feature").alias("type"))
        # Query accession column (fourth in tblout) was named "--" in our parser
        if "--" in ncrna_meta.columns:
            mappings.append(pl.col("--").alias("accession"))
        # Strand from tblout
        if "strand_ncrna" in ncrna_meta.columns:
            mappings.append(pl.col("strand_ncrna").alias("strand"))
        # E-value normalization
        if "E-value" in ncrna_meta.columns:
            mappings.append(pl.col("E-value").alias("evalue"))
        if mappings:
            ncrna_meta = ncrna_meta.with_columns(mappings)

        # Cast numeric fields if present
        castings = []
        if "score" in ncrna_meta.columns:
            castings.append(pl.col("score").cast(pl.Float64, strict=False).round(2).alias("score"))
        if "evalue" in ncrna_meta.columns:
            # Keep as string if scientific notation, else cast to Float64; best-effort cast
            castings.append(pl.col("evalue").cast(pl.Float64, strict=False).alias("evalue"))
        if castings:
            ncrna_meta = ncrna_meta.with_columns(castings)

        missing_cols = [c for c in ncrna_cols if c not in ncrna_meta.columns]
        if not missing_cols:
            ncrna_meta.select(ncrna_cols).write_csv(
                tsv_dir / "ncrna_metadata.txt", separator="\t", include_header=True
            )
            ncrna_meta.select(ncrna_cols).write_parquet(parquet_dir / "ncrna_metadata.parquet")
        else:
            # Backfill missing columns to maintain schema
            filled = ncrna_meta.select(
                [
                    pl.col(c) if c in ncrna_meta.columns else pl.lit(None).alias(c)
                    for c in ncrna_cols
                ]
            )
            filled.write_csv(tsv_dir / "ncrna_metadata.txt", separator="\t", include_header=True)
            filled.write_parquet(parquet_dir / "ncrna_metadata.parquet")
    else:
        empty_ncrna = pl.DataFrame({c: [] for c in ncrna_cols})
        empty_ncrna.write_csv(tsv_dir / "ncrna_metadata.txt", separator="\t", include_header=True)
        empty_ncrna.write_parquet(parquet_dir / "ncrna_metadata.parquet")

    if nt_links is not None and isinstance(nt_links, pl.DataFrame) and nt_links.height > 0:
        cols = [
            "query",
            "query_start",
            "query_end",
            "ref",
            "ref_start",
            "ref_end",
            "ani",
        ]
        present = [c for c in cols if c in nt_links.columns]
        if len(present) == len(cols):
            nt_links.select(cols).write_csv(
                tsv_dir / "nucleotide_links.txt", separator="\t", include_header=False
            )
            nt_links.select(cols).write_parquet(parquet_dir / "nucleotide_links.parquet")
        else:
            empty_nt = pl.DataFrame({c: [] for c in cols})
            empty_nt.write_csv(
                tsv_dir / "nucleotide_links.txt", separator="\t", include_header=False
            )
            empty_nt.write_parquet(parquet_dir / "nucleotide_links.parquet")
    else:
        empty_nt = pl.DataFrame(
            {
                c: []
                for c in ["query", "query_start", "query_end", "ref", "ref_start", "ref_end", "ani"]
            }
        )
        empty_nt.write_csv(tsv_dir / "nucleotide_links.txt", separator="\t", include_header=False)
        empty_nt.write_parquet(parquet_dir / "nucleotide_links.parquet")

    if pairwise_aa is not None and isinstance(pairwise_aa, pl.DataFrame) and pairwise_aa.height > 0:
        cols = ["qseqid", "sseqid", "pident"]
        present = [c for c in cols if c in pairwise_aa.columns]
        if len(present) == len(cols):
            pairwise_aa.select(cols).write_csv(
                tsv_dir / "protein_links.txt", separator="\t", include_header=False
            )
            pairwise_aa.select(cols).write_parquet(parquet_dir / "protein_links.parquet")
        else:
            empty_prot = pl.DataFrame({c: [] for c in cols})
            empty_prot.write_csv(
                tsv_dir / "protein_links.txt", separator="\t", include_header=False
            )
            empty_prot.write_parquet(parquet_dir / "protein_links.parquet")
    else:
        empty_prot = pl.DataFrame({c: [] for c in ["qseqid", "sseqid", "pident"]})
        empty_prot.write_csv(tsv_dir / "protein_links.txt", separator="\t", include_header=False)
        empty_prot.write_parquet(parquet_dir / "protein_links.parquet")

    # Render standalone HTML by injecting base64 parquet data into the template placeholders.
    resource_template = files("hoodini").joinpath("template", "template.html")
    template_html = resource_template.read_text(encoding="utf-8")

    def b64(path: Path) -> str:
        try:
            return base64.b64encode(Path(path).read_bytes()).decode()
        except Exception:
            return ""

    parquet_map = {
        "PARQUET_GFF_B64": parquet_dir / "gff.parquet",
        "PARQUET_PROT_LINKS_B64": parquet_dir / "protein_links.parquet",
        "PARQUET_NUC_LINKS_B64": parquet_dir / "nucleotide_links.parquet",
        "PARQUET_DOMAINS_B64": parquet_dir / "domains.parquet",
        "PARQUET_HOODS_B64": parquet_dir / "hoods.parquet",
        "PARQUET_PROT_META_B64": parquet_dir / "protein_metadata.parquet",
        "PARQUET_DOM_META_B64": parquet_dir / "domains_metadata.parquet",
        "PARQUET_TREE_META_B64": parquet_dir / "tree_metadata.parquet",
        "PARQUET_NCRNA_META_B64": parquet_dir / "ncrna_metadata.parquet",
    }

    env = Environment(
        autoescape=False,
        variable_start_string="##HOODINI:",
        variable_end_string="##!",
    )
    template = env.from_string(template_html)
    rendered_html = template.render(
        NEWICK_B64=base64.b64encode(gzip.compress((tree_str or "").encode())).decode(),
        **{k: b64(v) for k, v in parquet_map.items()},
    )

    # Create HTML with output name and timestamp
    output_name = Path(output_dir).name
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M")
    html_path = outdir / f"{output_name}_{timestamp}.html"
    html_path.write_text(rendered_html, encoding="utf-8")

    # If running in a Jupyter notebook, display the HTML interactively
    try:
        from IPython import get_ipython
        from IPython.display import IFrame, display

        ipython = get_ipython()
        if ipython is not None and "IPKernelApp" in ipython.config:
            # We're in a Jupyter notebook
            display(IFrame(src=str(html_path), width="100%", height=800))
    except (ImportError, AttributeError):
        # Not in a notebook or IPython not available
        pass

    return outdir
