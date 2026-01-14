#!/usr/bin/env python3
"""IPG (Identical Protein Groups) enrichment pipeline using Polars."""

from importlib.resources import files
from pathlib import Path

import polars as pl

from hoodini.pipeline.helpers.fetch_ipg_from_accessions import fetch_ipg_from_accessions
from hoodini.pipeline.helpers.nuc2asmlen import run_nuc2asmlen
from hoodini.utils.logging_utils import info, warn
from hoodini.utils.polars_adapters import to_polars

PlDF = pl.DataFrame


def safe_collect(lf: pl.LazyFrame) -> PlDF:
    """Collect a lazy frame, preferring streaming when available."""
    try:
        return lf.collect(engine="streaming")
    except (TypeError, AttributeError):
        return lf.collect(streaming=False)
    except Exception:
        return lf.collect(streaming=False)


def run_ipg(records_df: PlDF, *, cand_mode: str) -> PlDF:
    """Polars-based IPG enrichment pipeline.

    Expected Files:
    ---------------
    - records_df: DataFrame from initialize_inputs with protein IDs
    - Remote APIs: NCBI IPG web service, assembly reports (FTP)

    Generated Files:
    ----------------
    - None (all data kept in memory/DataFrame)

    Process:
    --------
    1. Fetches IPG (Identical Protein Groups) records from NCBI for input protein IDs
    2. Maps proteins to nucleotide IDs and assembly accessions
    3. Applies candidate selection mode ('all', 'best', or 'model')
    4. Enriches records with assembly metadata (taxid, organism, strain)
    5. Computes assembly lengths and filters candidates

    Returns:
    --------
    pl.DataFrame: Enriched records with columns: og_index, protein_id, nucleotide_id,
                  assembly, taxid, organism, strain, unique_id, etc.
    """
    df = records_df.clone()

    def _as_scalar_str(val):
        if isinstance(val, list):
            return val[0] if val else None
        return val

    norm_exprs = []
    for col in ("protein_id", "nucleotide_id"):
        if col in df.columns:
            norm_exprs.append(
                pl.col(col).map_elements(_as_scalar_str, return_dtype=pl.Utf8).alias(col)
            )
    if norm_exprs:
        df = df.with_columns(norm_exprs)

    if "query_protein_id" not in df.columns:
        df = df.with_columns(pl.col("protein_id").alias("query_protein_id"))

    def _is_refseq(val):
        s = _as_scalar_str(val)
        return bool(s) and str(s).startswith(("WP_", "YP_"))

    df = df.with_columns(
        pl.col("protein_id")
        .map_elements(_is_refseq, return_dtype=pl.Boolean)
        .alias("is_refseq_query")
    )
    info("🔍  Fetching IPG data...")
    df = _fetch_ipg_data(df, cand_mode)
    _trace_ipg(df, stage="after_fetch_ipg")
    info("🔍  Fetching nucleotide data...")
    df = _fetch_nucleotide_data(df)
    _trace_ipg(df, stage="after_fetch_nuc")
    info("✅  Selecting best IPG records...")
    df = _select_best_ipg(df, cand_mode)
    _trace_ipg(df, stage="after_select_best")
    df = _finalize_ipg(df, cand_mode)
    return df


def _fill_ipg_polars(records: PlDF, ipg_df: PlDF) -> PlDF:
    """Fill IPG information into records (Polars version)."""
    if ipg_df.height == 0:
        return records

    if "failed_reason" not in records.columns:
        records = records.with_columns(pl.lit(None).alias("failed_reason"))

    if "ipg_protein_id" not in ipg_df.columns and "protein_id" in ipg_df.columns:
        ipg_df = ipg_df.with_columns(pl.col("protein_id").alias("ipg_protein_id"))

    ipg_map = ipg_df.select(["protein_id", "ipg_id"]).unique(subset=["protein_id"])
    records = records.join(ipg_map, on="protein_id", how="left", suffix="_ipg")

    if "ipg_id_ipg" in records.columns:
        records = records.with_columns(
            pl.when(pl.col("ipg_id_ipg").is_not_null())
            .then(pl.col("ipg_id_ipg"))
            .otherwise(pl.col("ipg_id"))
            .alias("ipg_id")
        ).drop("ipg_id_ipg")

    cond = (
        (pl.col("protein_id").is_not_null())
        & (pl.col("failed").is_null())
        & (~pl.col("premade"))
        & (pl.col("ipg_id").is_null())
    )
    records = records.with_columns(
        pl.when(cond).then(True).otherwise(pl.col("failed")).alias("failed"),
        pl.when(cond)
        .then(pl.lit("Unable to retrieve IPG"))
        .otherwise(pl.col("failed_reason"))
        .alias("failed_reason"),
    )

    ipg_cols = [c for c in ipg_df.columns if c not in ["protein_id", "ipg_id"]]
    if "ipg_protein_id" in ipg_df.columns and "ipg_protein_id" not in ipg_cols:
        ipg_cols.append("ipg_protein_id")
    select_cols = ["ipg_id"] + ipg_cols
    if ipg_df.height > 0 and "ipg_id" in ipg_df.columns:
        select_cols_present = [c for c in select_cols if c in ipg_df.columns]
        if select_cols_present:
            ipg_subset = ipg_df.select(select_cols_present)
            records = records.join(ipg_subset, on="ipg_id", how="left", suffix="_ipg_data")

            for col in ipg_cols:
                ipg_col = f"{col}_ipg_data"
                if ipg_col in records.columns:
                    records = records.with_columns(
                        pl.when(pl.col(col).is_null())
                        .then(pl.col(ipg_col))
                        .otherwise(pl.col(col))
                        .alias(col)
                    ).drop(ipg_col)

    return records


def _fetch_ipg_data(df: PlDF, cand_mode: str) -> PlDF:
    """Fetch IPG data for proteins and enrich with assembly/taxid info."""
    cond_mask = (
        (pl.col("protein_id").is_not_null()) & (pl.col("failed").is_null()) & (~pl.col("premade"))
    )

    proteins = df.filter(cond_mask).select("protein_id").unique().to_series().to_list()
    proteins = [p for p in proteins if p and str(p).strip()]

    if not proteins:
        info("ℹ️  No records match conditions for IPG.")
        return df

    ipg_df = fetch_ipg_from_accessions(proteins)
    if ipg_df.height == 0:
        info("ℹ️  No IPG records found.")
        return df
    info(f"Fetched {ipg_df.height} IPG records for {len(proteins)} proteins.")

    assemblies = ipg_df.select("assembly").unique().to_series().to_list()
    assemblies = [a for a in assemblies if a and str(a).strip()]

    if not assemblies:
        info("ℹ️  No assemblies found in IPG data.")
        return df

    try:
        dive_path = files("hoodini").joinpath("data", "dive_combined.parquet")
        asm_ids = pl.DataFrame({"assembly_id": assemblies})
        ts = safe_collect(
            pl.scan_parquet(dive_path).join(asm_ids.lazy(), on="assembly_id", how="inner")
        )
        if ts.height > 0:
            ipg_df = ipg_df.join(ts, left_on="assembly", right_on="assembly_id", how="left")
    except Exception as e:
        warn(f"Skipping dive_combined.parquet: {e}")

    try:
        summary_path = files("hoodini").joinpath("data", "assembly_summary.parquet")
        asm_ids2 = pl.DataFrame({"assembly_accession": assemblies})
        summary = safe_collect(
            pl.scan_parquet(summary_path).join(
                asm_ids2.lazy(), on="assembly_accession", how="inner"
            )
        )
        keep_cols = [
            "assembly_accession",
            "taxid",
            "species_taxid",
            "organism_name",
            "infraspecific_name",
            "assembly_level",
            "group",
        ]
        summary = summary.select([c for c in keep_cols if c in summary.columns])
        if summary.height > 0:
            ipg_df = ipg_df.join(
                summary, left_on="assembly", right_on="assembly_accession", how="left"
            )
    except Exception as e:
        warn(f"Skipping assembly_summary.parquet: {e}")

    rename_map = {
        "nucleotide accession": "nucleotide_id",
        "assembly": "assembly_id",
        "stop": "end",
        "protein": "protein_id",
        "id": "ipg_id",
    }
    ipg_df = ipg_df.rename({k: v for k, v in rename_map.items() if k in ipg_df.columns})

    df = _fill_ipg_polars(records=df, ipg_df=ipg_df)

    from hoodini.utils.id_parsing import is_refseq_nuccore, switch_assembly_prefix

    def fix_asm(nuc_id: str | None, asm_id: str | None) -> str | None:
        if nuc_id is None or asm_id is None:
            return asm_id
        if (
            is_refseq_nuccore(nuc_id)
            and str(asm_id).startswith("GCA_")
            or not is_refseq_nuccore(nuc_id)
            and str(asm_id).startswith("GCF_")
        ):
            return switch_assembly_prefix(asm_id)
        return asm_id

    mask = (
        (pl.col("failed").is_null())
        & (~pl.col("premade"))
        & (pl.col("nucleotide_id").is_not_null())
    )

    df = df.with_columns(
        pl.when(mask)
        .then(
            pl.concat_list(["nucleotide_id", "assembly_id"]).map_elements(
                lambda x: fix_asm(x[0], x[1]), return_dtype=pl.Utf8
            )
        )
        .otherwise(pl.col("assembly_id"))
        .alias("assembly_id")
    )

    return df


def _fetch_nucleotide_data(df: PlDF) -> PlDF:
    """Fetch nucleotide sequence lengths and assembly metadata."""
    nucs = df.select("nucleotide_id").unique().to_series().to_list()
    nucs = [n for n in nucs if n and str(n).strip()]

    if not nucs:
        info("ℹ️  No nucleotide IDs to fetch.")
        return df

    base = files("hoodini").joinpath("data")
    contig_path = base / "contig_lengths"
    summary_path = base / "assembly_summary.parquet"

    try:
        contigs_lf = pl.scan_parquet(
            contig_path,
            missing_columns="insert",
            extra_columns="ignore",
        ).filter((pl.col("genbankAccession").is_in(nucs)) | (pl.col("refseqAccession").is_in(nucs)))

        summary_scan = pl.scan_parquet(summary_path).select(
            [
                "assembly_accession",
                "taxid",
                "species_taxid",
                "organism_name",
                "infraspecific_name",
                "assembly_level",
                "group",
            ]
        )

        joined_lf = contigs_lf.join(
            summary_scan, left_on="assemblyAccession", right_on="assembly_accession", how="left"
        )

        mapping_df = safe_collect(
            joined_lf.select(
                [
                    pl.col("assemblyAccession").alias("assembly_id"),
                    pl.col("length").alias("sequence_length"),
                    pl.col("taxid"),
                    pl.col("group"),
                    pl.col("species_taxid"),
                    pl.col("organism_name"),
                    pl.col("infraspecific_name"),
                    pl.col("assembly_level"),
                    pl.col("genbankAccession"),
                    pl.col("refseqAccession"),
                ]
            )
        )

    except Exception as e:
        warn(f"Failed Polars fetch in _fetch_nucleotide_data: {e}")
        mapping_df = pl.DataFrame()

    if mapping_df.height > 0:
        mapping_df = mapping_df.with_columns(
            pl.col("assembly_id").cast(pl.Utf8),
            pl.col("refseqAccession").str.replace(r"^[A-Z]{2}_", "").alias("refseq_no_prefix"),
            pl.col("genbankAccession").cast(pl.Utf8),
        )

        gcf_map = (
            mapping_df.filter(pl.col("assembly_id").str.starts_with("GCF_"))
            .select(
                [
                    pl.col("refseqAccession").alias("nucleotide_id"),
                    pl.col("refseq_no_prefix").alias("nucleotide_id_no_prefix"),
                    "assembly_id",
                    "sequence_length",
                    "taxid",
                    "group",
                    "species_taxid",
                    "organism_name",
                    "infraspecific_name",
                    "assembly_level",
                ]
            )
            .unique(subset=["nucleotide_id", "nucleotide_id_no_prefix"], keep="first")
        )

        gca_map = (
            mapping_df.filter(pl.col("assembly_id").str.starts_with("GCA_"))
            .select(
                [
                    pl.col("genbankAccession").alias("nucleotide_id"),
                    "assembly_id",
                    "sequence_length",
                    "taxid",
                    "group",
                    "species_taxid",
                    "organism_name",
                    "infraspecific_name",
                    "assembly_level",
                ]
            )
            .unique(subset=["nucleotide_id"], keep="first")
        )

        for c in [
            "sequence_length",
            "assembly_id",
            "taxid",
            "group",
            "species_taxid",
            "organism_name",
            "infraspecific_name",
            "assembly_level",
        ]:
            if c not in df.columns:
                df = df.with_columns(
                    pl.lit(None).cast(pl.Utf8 if c not in ["taxid"] else pl.Int64).alias(c)
                )

        if "nucleotide_id_no_prefix" not in df.columns:
            df = df.with_columns(
                pl.col("nucleotide_id")
                .str.replace(r"^[A-Z]{2}_", "")
                .alias("nucleotide_id_no_prefix")
            )

        has_asm_mask = df["assembly_id"].is_not_null()
        df_with_asm = df.filter(has_asm_mask)
        df_no_asm = df.filter(~has_asm_mask)

        if df_with_asm.height > 0:
            asm_meta = mapping_df.select(
                [
                    pl.col("assembly_id"),
                    "sequence_length",
                    "taxid",
                    "group",
                    "species_taxid",
                    "organism_name",
                    "infraspecific_name",
                    "assembly_level",
                ]
            ).unique(subset=["assembly_id"], keep="first")
            df_with_asm = df_with_asm.join(asm_meta, on="assembly_id", how="left", suffix="_asm")
            for col in [
                "sequence_length",
                "taxid",
                "group",
                "species_taxid",
                "organism_name",
                "infraspecific_name",
                "assembly_level",
            ]:
                src = f"{col}_asm"
                if src in df_with_asm.columns:
                    df_with_asm = df_with_asm.with_columns(
                        pl.when(pl.col(col).is_null())
                        .then(pl.col(src))
                        .otherwise(pl.col(col))
                        .alias(col)
                    ).drop(src)

        if df_no_asm.height > 0:
            df_no_asm = df_no_asm.join(
                gcf_map,
                left_on="nucleotide_id",
                right_on="nucleotide_id",
                how="left",
                suffix="_gcf",
            )
            df_no_asm = df_no_asm.join(
                gcf_map,
                left_on="nucleotide_id_no_prefix",
                right_on="nucleotide_id_no_prefix",
                how="left",
                suffix="_gcfnp",
            )

            for col in [
                "assembly_id",
                "sequence_length",
                "taxid",
                "group",
                "species_taxid",
                "organism_name",
                "infraspecific_name",
                "assembly_level",
            ]:
                for suffix_col in ("_gcf", "_gcfnp"):
                    src = f"{col}{suffix_col}"
                    if src in df_no_asm.columns:
                        df_no_asm = df_no_asm.with_columns(
                            pl.when(pl.col(col).is_null())
                            .then(pl.col(src))
                            .otherwise(pl.col(col))
                            .alias(col)
                        ).drop(src)

            df_no_asm = df_no_asm.join(gca_map, on="nucleotide_id", how="left", suffix="_gca")
            for col in gca_map.columns:
                if col != "nucleotide_id" and f"{col}_gca" in df_no_asm.columns:
                    df_no_asm = df_no_asm.with_columns(
                        pl.when(pl.col(col).is_null())
                        .then(pl.col(f"{col}_gca"))
                        .otherwise(pl.col(col))
                        .alias(col)
                    ).drop(f"{col}_gca")

        if has_asm_mask.any():
            all_cols = sorted(set(df_with_asm.columns) | set(df_no_asm.columns))

            dtype_map = {col: df_with_asm[col].dtype for col in df_with_asm.columns}
            for col in df_no_asm.columns:
                dtype_map.setdefault(col, df_no_asm[col].dtype)

            def _align(df_part: pl.DataFrame, cols: list[str]) -> pl.DataFrame:
                missing = [c for c in cols if c not in df_part.columns]
                if missing:
                    df_part = df_part.with_columns(
                        [pl.lit(None).cast(dtype_map[c]).alias(c) for c in missing]
                    )
                return df_part.select(cols)

            df = pl.concat(
                [
                    _align(df_with_asm, all_cols),
                    _align(df_no_asm, all_cols),
                ],
                how="vertical",
            )

        if "source" in df.columns:
            gca_override = gca_map.with_columns(
                pl.col("nucleotide_id").str.replace(r"^[A-Z]{2}_", "").alias("nuc_gca_no_prefix")
            ).select(
                [
                    pl.col("nucleotide_id").alias("nuc_gca_full"),
                    pl.col("nuc_gca_no_prefix"),
                    pl.col("assembly_id").alias("assembly_id_gca_override"),
                    pl.col("sequence_length").alias("sequence_length_gca_override"),
                    pl.col("taxid").alias("taxid_gca_override"),
                    pl.col("group").alias("group_gca_override"),
                    pl.col("species_taxid").alias("species_taxid_gca_override"),
                    pl.col("organism_name").alias("organism_name_gca_override"),
                    pl.col("infraspecific_name").alias("infraspecific_name_gca_override"),
                    pl.col("assembly_level").alias("assembly_level_gca_override"),
                ]
            )
            df = df.join(
                gca_override,
                left_on="nucleotide_id",
                right_on="nuc_gca_full",
                how="left",
                suffix="_gcaovr",
            )
            df = df.join(
                gca_override,
                left_on="nucleotide_id_no_prefix",
                right_on="nuc_gca_no_prefix",
                how="left",
                suffix="_gcaovrnp",
            )

            override_cols = [
                "assembly_id",
                "sequence_length",
                "taxid",
                "group",
                "species_taxid",
                "organism_name",
                "infraspecific_name",
                "assembly_level",
            ]
            for base in override_cols:
                for src in (f"{base}_gca_override", f"{base}_gca_override_gcaovrnp"):
                    if src in df.columns:
                        df = df.with_columns(
                            pl.when((pl.col("source") == "INSDC") & pl.col(src).is_not_null())
                            .then(pl.col(src))
                            .otherwise(pl.col(base))
                            .alias(base)
                        ).drop(src)

    missing_mask = df.select(pl.col("sequence_length").is_null()).to_series()
    if missing_mask.any():
        to_fetch = df.filter(missing_mask).select("nucleotide_id").unique().to_series().to_list()
        to_fetch = [t for t in to_fetch if t and str(t).strip()]
        if to_fetch:
            meta = run_nuc2asmlen(to_fetch)
            meta = to_polars(meta)
            if meta.height > 0:
                meta = meta.rename(
                    {
                        "NucleotideAccession": "nucleotide_id",
                        "length": "sequence_length",
                        "AssemblyAccession": "assembly_id",
                    }
                ).unique(subset=["nucleotide_id"], keep="first")
                df = df.join(meta, on="nucleotide_id", how="left", suffix="_edirect")
                for col in ["sequence_length", "assembly_id"]:
                    if f"{col}_edirect" in df.columns:
                        df = df.with_columns(
                            pl.when(pl.col(col).is_null())
                            .then(pl.col(f"{col}_edirect"))
                            .otherwise(pl.col(col))
                            .alias(col)
                        ).drop(f"{col}_edirect")

                new_asms = (
                    df.filter(pl.col("assembly_id").is_not_null())
                    .select("assembly_id")
                    .unique()
                    .to_series()
                    .to_list()
                )
                if new_asms:
                    try:
                        ids_asm = pl.DataFrame({"assembly_accession": new_asms})
                        summary2 = safe_collect(
                            pl.scan_parquet(summary_path).join(
                                ids_asm.lazy(), on="assembly_accession", how="inner"
                            )
                        )
                        if summary2.height > 0:
                            summary2 = summary2.select(
                                ["assembly_accession", "taxid", "group"]
                            ).rename({"assembly_accession": "assembly_id"})
                            df = df.join(summary2, on="assembly_id", how="left", suffix="_backfill")
                            for col in ["taxid", "group"]:
                                if f"{col}_backfill" in df.columns:
                                    df = df.with_columns(
                                        pl.when(pl.col(col).is_null())
                                        .then(pl.col(f"{col}_backfill"))
                                        .otherwise(pl.col(col))
                                        .alias(col)
                                    ).drop(f"{col}_backfill")
                    except Exception as e:
                        warn(f"Failed backfill join: {e}")

    cond = (
        (pl.col("input_type") == "nucleotide")
        & (pl.col("nucleotide_id").is_not_null())
        & (pl.col("failed").is_null())
        & (~pl.col("premade"))
        & (pl.col("start").is_null())
        & (pl.col("end").is_null())
    )
    df = df.with_columns(
        pl.when(cond).then(pl.lit(0)).otherwise(pl.col("start")).alias("start"),
        pl.when(cond).then(pl.col("sequence_length")).otherwise(pl.col("end")).alias("end"),
    )

    return df


def _select_best_ipg(df: pl.DataFrame, cand_mode: str) -> pl.DataFrame:
    """Select best IPG record per og_index based on ranking."""
    if "ipg_id" not in df.columns:
        return df

    mask_ipg = df.select(pl.col("ipg_id").is_not_null()).to_series()
    df_no_ipg = df.filter(~mask_ipg)
    df_ipg = df.filter(mask_ipg)

    if df_ipg.height == 0:
        return df

    def _require_query_match(df_in: pl.DataFrame) -> pl.DataFrame:
        """Keep rows whose protein_id matches the original query when present for that og_index."""
        protein_col = "ipg_protein_id" if "ipg_protein_id" in df_in.columns else "protein_id"
        df_in = df_in.with_columns(
            (pl.col(protein_col) == pl.col("query_protein_id")).alias("_is_query")
        )
        flags = df_in.group_by("og_index").agg(pl.any("_is_query").alias("_has_query"))
        df_in = df_in.join(flags, on="og_index", how="left")
        df_in = df_in.filter(
            pl.when(pl.col("_has_query")).then(pl.col("_is_query")).otherwise(True)
        )
        return df_in.drop(["_is_query", "_has_query"])

    def _rank_func(asm, lvl):
        if asm and str(asm).startswith("GCF") and lvl in ["Chromosome", "Complete Genome"]:
            return 1
        elif lvl in ["Chromosome", "Complete Genome"]:
            return 2
        elif asm:
            return 3
        else:
            return 4

    if cand_mode == "any_ipg":
        df_ipg = df_ipg.unique(subset=["ipg_id", "nucleotide_id", "start", "end"], keep="first")
        df_ipg = df_ipg.with_columns(
            pl.col("assembly_id").cast(pl.Utf8, strict=False).alias("assembly_id"),
            pl.col("assembly_id")
            .cast(pl.Utf8, strict=False)
            .str.replace(r"^(GCF_|GCA_)", "")
            .alias("_asm_core"),
            pl.col("assembly_id")
            .cast(pl.Utf8, strict=False)
            .str.starts_with("GCF_")
            .alias("_is_refseq"),
        )
        has_refseq = df_ipg.group_by(["_asm_core", "start", "end"]).agg(
            pl.any("_is_refseq").alias("_has_refseq")
        )
        df_ipg = df_ipg.join(has_refseq, on=["_asm_core", "start", "end"], how="left")
        df_ipg = df_ipg.filter(
            pl.when(pl.col("_has_refseq")).then(pl.col("_is_refseq")).otherwise(True)
        )
        df_ipg = df_ipg.drop(["_asm_core", "_is_refseq", "_has_refseq"])
    elif cand_mode in ["best_ipg", "best_id"]:
        if cand_mode == "best_id":
            df_ipg = _require_query_match(df_ipg)
        df_ipg = df_ipg.with_columns(
            pl.concat_list(["assembly_id", "assembly_level"])
            .map_elements(lambda x: _rank_func(x[0], x[1]), return_dtype=pl.Int32)
            .alias("ranked")
        )
        df_ipg = df_ipg.with_columns(
            pl.when(
                pl.all_horizontal(
                    pl.col("start").is_not_null(),
                    pl.col("end").is_not_null(),
                    pl.col("sequence_length").is_not_null(),
                )
            )
            .then(
                pl.min_horizontal(
                    pl.col("start").cast(pl.Int64, strict=False),
                    (
                        pl.col("sequence_length").cast(pl.Int64, strict=False)
                        - pl.col("end").cast(pl.Int64, strict=False)
                    ),
                )
            )
            .otherwise(None)
            .alias("_edge_buffer")
        )
        df_ipg = df_ipg.with_columns(
            pl.col("_edge_buffer").fill_null(-1).cast(pl.Int64, strict=False).alias("_edge_buffer")
        )
        df_ipg = (
            df_ipg.sort(["og_index", "ranked", "_edge_buffer"], descending=[False, False, True])
            .unique(subset=["og_index"], keep="first")
            .drop(["ranked", "_edge_buffer"])
        )
    elif cand_mode == "one_id":
        df_ipg = df_ipg.unique(subset=["og_index"], keep="first")
    elif cand_mode == "same_id":
        df_ipg = _require_query_match(df_ipg)
        df_ipg = df_ipg.unique(subset=["og_index", "nucleotide_id", "start", "end"], keep="first")

    return pl.concat([df_no_ipg, df_ipg], how="vertical")


def _finalize_ipg(df: pl.DataFrame, cand_mode: str) -> pl.DataFrame:
    """Mark records as failed if they don't meet requirements."""
    for col in [
        "group",
        "gff_path",
        "faa_path",
        "assembly_id",
        "failed",
        "premade",
        "failed_reason",
    ]:
        if col not in df.columns:
            df = df.with_columns(pl.lit(None).alias(col))

    cond1 = (
        (pl.col("failed").is_null())
        & (
            ~(
                (pl.col("gff_path").is_not_null() & pl.col("faa_path").is_not_null())
                | pl.col("assembly_id").is_null()
            )
        )
        & (~pl.col("premade"))
    )
    to_fail1 = cond1 & pl.col("assembly_id").is_null()
    df = df.with_columns(
        pl.when(to_fail1).then(True).otherwise(pl.col("failed")).alias("failed"),
        pl.when(to_fail1)
        .then(pl.lit("Unable to retrieve IPG/Nuccore data"))
        .otherwise(pl.col("failed_reason"))
        .alias("failed_reason"),
    )

    valid = {"bacteria", "viral", "archaea", "metagenomes"}
    cond2 = (
        (pl.col("failed").is_null())
        & (
            ~(
                (pl.col("gff_path").is_not_null() & pl.col("faa_path").is_not_null())
                | pl.col("assembly_id").is_null()
            )
        )
        & (~pl.col("premade"))
    )
    to_invalid = cond2 & (~pl.col("group").is_in(valid))
    df = df.with_columns(
        pl.when(to_invalid).then(True).otherwise(pl.col("failed")).alias("failed"),
        pl.when(to_invalid)
        .then(pl.lit("Invalid superkingdom"))
        .otherwise(pl.col("failed_reason"))
        .alias("failed_reason"),
    )

    return df


def _trace_ipg(df: pl.DataFrame, stage: str) -> None:
    """Write a small debug trace of IPG enrichment when debug mode is on."""
    from hoodini.utils.logging_utils import is_debug_enabled

    if not is_debug_enabled():
        return
    cols = [
        "unique_id",
        "og_index",
        "protein_id",
        "nucleotide_id",
        "assembly_id",
        "source",
        "ipg_id",
        "start",
        "end",
        "failed",
    ]
    existing = [c for c in cols if c in df.columns]
    trace_df = df.select(existing).with_columns(pl.lit(stage).alias("stage"))
    out_path = Path("ipg_debug.csv")
    include_header = not out_path.exists()
    with out_path.open("a") as fh:
        trace_df.write_csv(fh, include_header=include_header)
