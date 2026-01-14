from __future__ import annotations

import concurrent.futures as _fut
import math
from pathlib import Path

import numpy as np
import polars as pl

from hoodini.utils.logging_utils import console

try:
    from scipy.stats import hypergeom

    _HAS_SCIPY = True
except Exception:
    _HAS_SCIPY = False


def _ensure_polars(df_or_lazy) -> pl.DataFrame:
    """Convert eager/lazy Polars to eager Polars. Reject other types to avoid pandas dependency."""
    if isinstance(df_or_lazy, pl.DataFrame):
        return df_or_lazy
    if isinstance(df_or_lazy, pl.LazyFrame):
        return df_or_lazy.collect()
    raise TypeError(f"Unsupported DataFrame type: {type(df_or_lazy)}")


def _normalize_hits_df(hits: pl.DataFrame) -> pl.DataFrame:
    """
    Normalize a BLAST outfmt 6 table to internal schema:
      keep qseqid / sseqid and pident (float). Cast types for joins.
    Other columns (length, evalue, bitscore) are preserved but unused here.
    """
    need = {"qseqid", "sseqid", "pident"}
    missing = need - set(hits.columns)
    if missing:
        raise ValueError(
            f"hits_df is missing required columns: {sorted(missing)}. " f"Available: {hits.columns}"
        )
    return hits.with_columns(
        pl.col("qseqid").cast(pl.Utf8, strict=False),
        pl.col("sseqid").cast(pl.Utf8, strict=False),
        pl.col("pident").cast(pl.Float64, strict=False),
    )


def _normalize_proteins_df(prots: pl.DataFrame, require_fam: bool = False) -> pl.DataFrame:
    """
    Normalize protein annotation table to internal schema:
      prot_id, target_nuc, optional fam_cluster.
    Accepts common aliases.
    """
    cols = set(prots.columns)

    prot_col = next(
        (c for c in ["prot_id", "protein_id", "protid", "proteinId", "proteinID"] if c in cols),
        None,
    )
    nuc_col = next(
        (
            c
            for c in [
                "target_nuc",
                "nucleotide_id",
                "seqid",
                "seq_id",
                "contig_id",
                "contig",
                "scaffold_id",
                "target_seq",
            ]
            if c in cols
        ),
        None,
    )
    fam_col = next(
        (
            c
            for c in [
                "fam_cluster",
                "fam_clustter",
                "family_cluster",
                "cluster",
                "protein_family",
                "pfam_cluster",
                "fam",
            ]
            if c in cols
        ),
        None,
    )

    miss = []
    if prot_col is None:
        miss.append("prot_id (aliases: protein_id, protid, ...)")
    if nuc_col is None:
        miss.append("target_nuc (aliases: nucleotide_id, seqid, contig_id, ...)")
    if miss:
        raise ValueError(
            f"proteins_df is missing required columns: {miss}. Available: {prots.columns}"
        )

    out = prots.rename({prot_col: "prot_id", nuc_col: "target_nuc"})
    if fam_col:
        out = out.rename({fam_col: "fam_cluster"})
    elif require_fam:
        raise ValueError(
            "proteins_df needs a family column (e.g., fam_cluster) for the vContact2-like metric."
        )
    return out


def _annotate_hits(
    hits: pl.DataFrame, prots: pl.DataFrame, pident_min: float | None, exclude_self: bool
) -> pl.DataFrame:
    """Join hits to seq IDs and apply filters in Polars."""
    qp = prots.rename({"prot_id": "qseqid", "target_nuc": "q_seq"})
    tp = prots.rename({"prot_id": "sseqid", "target_nuc": "t_seq"})

    ann = (
        hits.join(qp, on="qseqid", how="left")
        .join(tp, on="sseqid", how="left")
        .drop_nulls(["q_seq", "t_seq"])
    )
    if pident_min is not None:
        ann = ann.filter(pl.col("pident") >= float(pident_min))
    if exclude_self:
        ann = ann.filter(pl.col("q_seq") != pl.col("t_seq"))
    return ann


def _sizes_per_seq(prots: pl.DataFrame, col_name="n_prots") -> pl.DataFrame:
    return prots.group_by("target_nuc").agg(pl.col("prot_id").n_unique().alias(col_name))


def _find_prot_id_col(prots: pl.DataFrame) -> str:
    """Return the protein identifier column name present in `prots`.
    Prefer common names: 'prot_id', 'protein_id', 'id'."""
    for c in ("prot_id", "protein_id", "id", "protid"):
        if c in prots.columns:
            return c
    raise ValueError(f"cannot find protein id column in proteins table. Available: {prots.columns}")


def _compute_subset_protein_ids(
    all_prots: pl.DataFrame,
    all_neigh: pl.DataFrame,
    all_gff: pl.DataFrame,
    subset_mode: str,
    win: int | None = None,
    win_mode: str = "bp",
) -> set:
    """Compute a set of protein ids to keep according to subset_mode.

    Modes supported:
      - 'target_prot': keep unique values from column 'target_prot' in `all_prots`.
      - 'target_region': keep all proteins whose coordinates (in all_gff) fall within
         the neighborhood window (start_win..end_win) from `all_neigh`.
      - 'window': like 'target_region' but expand the neighborhood window by `win` on
         both sides before selecting proteins.
    """
    mode = (subset_mode or "").lower()
    if mode not in {"target_prot", "target_region", "window"}:
        raise ValueError(f"Unknown subset_mode={subset_mode}")

    if mode == "target_prot":
        if "target_prot" not in all_prots.columns:
            raise ValueError(
                "'target_prot' column not found in all_prots for subset_mode='target_prot'"
            )
        vals = all_prots.select("target_prot").drop_nulls().unique().to_series().to_list()
        return {v for v in vals if v is not None}

    if not {"seqid", "start_win", "end_win"}.issubset(set(all_neigh.columns)):
        raise ValueError(
            "all_neigh must contain columns: seqid, start_win, end_win for region/window subsetting"
        )
    prot_col = None
    for c in ("protein_id", "prot_id", "id"):
        if c in all_gff.columns:
            prot_col = c
            break
    if prot_col is None:
        if "attributes" in all_gff.columns:
            temp_gff = all_gff.with_columns(
                pl.col("attributes").str.extract(r"ID=([^;]+)").alias("id")
            )
            prot_col = "id"
            all_gff = temp_gff
        else:
            raise ValueError(
                "all_gff must contain a protein id column (protein_id/id) or an 'attributes' column with ID= to perform region/window subsetting"
            )

    neigh = all_neigh.select(["seqid", "start_win", "end_win"]).drop_nulls()
    gff = all_gff.select(["seqid", "start", "end", prot_col]).drop_nulls()

    if mode == "window":
        w = int(win or 0)
        wm = (win_mode or "").lower()
        gene_mode = wm in ("genes", "win_genes", "win-genes", "win_genes", "gene", "win_gene")
        if gene_mode:
            gff = gff.with_columns(
                ((pl.col("start").rank(method="ordinal").over("seqid") - 1).cast(pl.Int64)).alias(
                    "gene_idx"
                )
            )

            joined = gff.join(neigh, on="seqid", how="inner")
            contained = joined.filter(
                (pl.col("start") >= pl.col("start_win")) & (pl.col("end") <= pl.col("end_win"))
            )
            if contained.is_empty():
                return set()
            min_idx = int(contained.select(pl.col("gene_idx").min()).to_series().iloc[0])
            max_idx = int(contained.select(pl.col("gene_idx").max()).to_series().iloc[0])
            start_idx = max(0, min_idx - w)
            end_idx = max_idx + w
            sel = gff.filter((pl.col("gene_idx") >= start_idx) & (pl.col("gene_idx") <= end_idx))
        else:
            neigh = neigh.with_columns(
                [
                    (pl.col("start_win") - w).alias("start_exp"),
                    (pl.col("end_win") + w).alias("end_exp"),
                ]
            )
            neigh = neigh.with_columns(
                [
                    pl.when(pl.col("start_exp") < 0)
                    .then(0)
                    .otherwise(pl.col("start_exp"))
                    .alias("start_exp"),
                ]
            )
            joined = gff.join(neigh, on="seqid", how="inner")
            sel = joined.filter(
                (pl.col("start") >= pl.col("start_exp")) & (pl.col("end") <= pl.col("end_exp"))
            )
    else:
        joined = gff.join(neigh, on="seqid", how="inner")
        sel = joined.filter(
            (pl.col("start") >= pl.col("start_win")) & (pl.col("end") <= pl.col("end_win"))
        )

    if sel.is_empty():
        return set()
    vals = sel.select(prot_col).unique().to_series().to_list()
    return {v for v in vals if v is not None}


def _rbh_from_ann(ann: pl.DataFrame) -> pl.DataFrame:
    """
    Compute Reciprocal Best Hits fully in Polars.
    Tie-breaking: ordinal rank => first occurrence wins within each group.
    Returns columns: query_seq, target_seq, qseqid, sseqid, pident
    """
    best_q_to_target = (
        ann.with_columns(
            pl.col("pident")
            .rank(method="ordinal", descending=True)
            .over(["qseqid", "t_seq"])
            .alias("_rk")
        )
        .filter(pl.col("_rk") == 1)
        .select(
            [
                pl.col("qseqid"),
                pl.col("t_seq").alias("target_seq"),
                pl.col("sseqid"),
                pl.col("pident").alias("pident_best_qt"),
            ]
        )
    )

    best_t_to_query = (
        ann.with_columns(
            pl.col("pident")
            .rank(method="ordinal", descending=True)
            .over(["sseqid", "q_seq"])
            .alias("_rk")
        )
        .filter(pl.col("_rk") == 1)
        .select(
            [
                pl.col("sseqid"),
                pl.col("q_seq").alias("query_seq"),
                pl.col("qseqid"),
                pl.col("pident").alias("pident_best_tq"),
            ]
        )
    )

    rbh = best_q_to_target.join(
        best_t_to_query,
        left_on=["qseqid", "sseqid"],
        right_on=["qseqid", "sseqid"],
        how="inner",
        suffix="_r",
    ).select(
        [
            pl.col("query_seq"),
            pl.col("target_seq"),
            pl.col("qseqid"),
            pl.col("sseqid"),
            pl.col("pident_best_qt").alias("pident"),
        ]
    )

    return rbh


def compute_wgrr(
    hits_df: pl.DataFrame | pl.LazyFrame,
    proteins_df: pl.DataFrame | pl.LazyFrame,
    pident_min: float | None = 30.0,
    exclude_self: bool = True,
    symmetric: bool = True,
) -> pl.DataFrame:
    """
    Weighted Gene Repertoire Relatedness (wGRR) between proteomes.

    Definition (symmetric by construction):
      - Build RBHs between proteomes A and B.
      - Sum identity FRACTIONS (pident/100) across RBH pairs per sequence pair.
      - Normalize by min(nA, nB).
      - Clip to [0, 1].

    Returns Polars DataFrame with columns: ["qseqid", "sseqid", "wGRR_sym", "AAI"]
    """
    hits_raw = _ensure_polars(hits_df)
    prots_raw = _ensure_polars(proteins_df)

    hits = _normalize_hits_df(hits_raw)
    prots = _normalize_proteins_df(prots_raw)

    sizes = _sizes_per_seq(prots, col_name="n_prots")

    ann = _annotate_hits(hits, prots, pident_min=pident_min, exclude_self=exclude_self)
    rbh = _rbh_from_ann(ann)
    if rbh.is_empty():
        return pl.DataFrame(columns=["qseqid", "sseqid", "wGRR_sym", "AAI"])

    w = (
        rbh.group_by(["query_seq", "target_seq"])
        .agg(pl.col("pident").sum().alias("sum_pident_pct"))
        .join(sizes.rename({"target_nuc": "query_seq"}), on="query_seq", how="left")
        .join(sizes.rename({"target_nuc": "target_seq"}), on="target_seq", how="left", suffix="_B")
        .with_columns(
            [
                (pl.col("sum_pident_pct") / 100.0).alias("sum_pident_frac"),
                pl.min_horizontal(pl.col(["n_prots", "n_prots_B"])).alias("n_min"),
            ]
        )
        .with_columns(
            (
                pl.when(pl.col("n_min") > 0)
                .then(pl.col("sum_pident_frac") / pl.col("n_min"))
                .otherwise(0.0)
            ).alias("wGRR_sym")
        )
        .with_columns(pl.col("wGRR_sym").clip(0.0, 1.0))
        .select(["query_seq", "target_seq", "wGRR_sym"])
    )
    w = w.with_columns(pl.col("wGRR_sym").alias("AAI"))
    return w.rename({"query_seq": "qseqid", "target_seq": "sseqid"}).select(
        ["qseqid", "sseqid", "wGRR_sym", "AAI"]
    )


def compute_aai_rbh(
    hits_df: pl.DataFrame | pl.LazyFrame,
    proteins_df: pl.DataFrame | pl.LazyFrame,
    pident_min: float | None = 30.0,
    exclude_self: bool = True,
) -> pl.DataFrame:
    """
    Average Amino-acid Identity using Reciprocal Best Hits (RBH).
    Returns Polars DataFrame.
    """
    hits_raw = _ensure_polars(hits_df)
    prots_raw = _ensure_polars(proteins_df)

    hits = _normalize_hits_df(hits_raw)
    prots = _normalize_proteins_df(prots_raw)

    ann = _annotate_hits(hits, prots, pident_min=pident_min, exclude_self=exclude_self)
    rbh = _rbh_from_ann(ann)
    if rbh.is_empty():
        return pl.DataFrame(columns=["qseqid", "sseqid", "AAI", "n_RBH", "RBH_frac_min"])

    aai = rbh.group_by(["query_seq", "target_seq"]).agg(
        [
            pl.len().alias("n_RBH"),
            pl.col("pident").mean().alias("AAI"),
        ]
    )

    sizes = _sizes_per_seq(prots, col_name="n_prots")
    aai = (
        aai.join(sizes.rename({"target_nuc": "query_seq"}), on="query_seq", how="left")
        .join(sizes.rename({"target_nuc": "target_seq"}), on="target_seq", how="left", suffix="_B")
        .with_columns(
            (pl.col("n_RBH") / pl.min_horizontal(pl.col(["n_prots", "n_prots_B"]))).alias(
                "RBH_frac_min"
            )
        )
        .select(["query_seq", "target_seq", "AAI", "n_RBH", "RBH_frac_min"])
    )
    return aai.rename({"query_seq": "qseqid", "target_seq": "sseqid"}).select(
        ["qseqid", "sseqid", "AAI", "n_RBH", "RBH_frac_min"]
    )


def compute_vcontact2_hypergeom(
    proteins_df: pl.DataFrame | pl.LazyFrame,
    max_df_frac: float = 0.2,
    min_shared: int = 2,
    multiple_test: str = "BH",
) -> pl.DataFrame:
    """
    vContact2-like hypergeometric similarity on presence/absence of protein families.
    Vectorized version (no Python loops). Returns Polars DataFrame.
    """
    prots_raw = _ensure_polars(proteins_df)

    fam_like = [
        c
        for c in [
            "fam_cluster",
            "fam_clustter",
            "family_cluster",
            "cluster",
            "protein_family",
            "pfam_cluster",
            "fam",
        ]
        if c in prots_raw.columns
    ]
    if not fam_like:
        return pl.DataFrame(
            columns=["query_seq", "target_seq", "k", "K_A", "K_B", "M", "pval", "p_adj", "score"]
        )

    prots = _normalize_proteins_df(prots_raw, require_fam=True)

    df_pa = prots.select(["target_nuc", "fam_cluster"]).drop_nulls().unique()
    if df_pa.is_empty():
        return pl.DataFrame(
            columns=["query_seq", "target_seq", "k", "K_A", "K_B", "M", "pval", "p_adj", "score"]
        )

    fam_counts = df_pa.group_by("fam_cluster").agg(pl.len().alias("df"))
    Nseqs = df_pa["target_nuc"].n_unique()
    keep_thresh = int(math.ceil(max_df_frac * max(1, Nseqs)))
    keep_fams = fam_counts.filter(pl.col("df") <= keep_thresh).select("fam_cluster")
    if keep_fams.is_empty():
        return pl.DataFrame(
            columns=["query_seq", "target_seq", "k", "K_A", "K_B", "M", "pval", "p_adj", "score"]
        )

    df_pa = df_pa.join(keep_fams, on="fam_cluster", how="inner")

    sizes = df_pa.group_by("target_nuc").agg(pl.col("fam_cluster").n_unique().alias("K"))

    M = int(df_pa["fam_cluster"].n_unique())

    pairs = (
        df_pa.join(df_pa, on="fam_cluster", how="inner", suffix="_r")
        .filter(pl.col("target_nuc") < pl.col("target_nuc_r"))
        .group_by(["target_nuc", "target_nuc_r"])
        .agg(pl.len().alias("k"))
        .rename({"target_nuc": "query_seq", "target_nuc_r": "target_seq"})
        .filter(pl.col("k") >= int(min_shared))
    )
    if pairs.is_empty():
        return pl.DataFrame(
            columns=["query_seq", "target_seq", "k", "K_A", "K_B", "M", "pval", "p_adj", "score"]
        )

    pairs = (
        pairs.join(
            sizes.rename({"target_nuc": "query_seq", "K": "K_A"}), on="query_seq", how="left"
        )
        .join(sizes.rename({"target_nuc": "target_seq", "K": "K_B"}), on="target_seq", how="left")
        .with_columns(pl.lit(M).alias("M"))
    )

    if pairs.is_empty():
        return pl.DataFrame(
            columns=[
                "query_seq",
                "target_seq",
                "k",
                "K_A",
                "K_B",
                "M",
                "pval",
                "p_adj",
                "score",
                "AAI",
            ]
        )

    try:
        k_vals = pairs["k"].to_numpy()
        Ka = pairs["K_A"].to_numpy()
        Kb = pairs["K_B"].to_numpy()
        Mvals = pairs["M"].to_numpy()

        if _HAS_SCIPY:
            pvals = hypergeom(Mvals, Ka, Kb).sf(k_vals - 1)
        else:
            from math import comb

            def hypergeom_sf(k, Mv, K, n):
                denom = comb(int(Mv), int(n))
                top = 0
                for i in range(int(k), int(min(K, n)) + 1):
                    top += comb(int(K), int(i)) * comb(int(Mv - K), int(n - i))
                return top / denom if denom != 0 else 1.0

            pvals = np.array(
                [
                    hypergeom_sf(k, Mv, Ka_i, Kb_i)
                    for k, Ka_i, Kb_i, Mv in zip(k_vals, Ka, Kb, Mvals)
                ]
            )
    except Exception:
        pvals = np.ones(pairs.height, dtype=float)

    pvals = np.clip(pvals, 1e-300, 1.0)
    mtests = len(pvals)

    if multiple_test.lower() == "bonferroni":
        padj = np.minimum(1.0, pvals * mtests)
    else:
        order = np.argsort(pvals)
        ranks = np.empty_like(order)
        ranks[order] = np.arange(1, mtests + 1)
        padj = np.minimum(1.0, pvals * mtests / ranks)

    score = -np.log10(np.maximum(1e-300, padj))

    pairs = pairs.with_columns(
        pl.Series("pval", pvals),
        pl.Series("p_adj", padj),
        pl.Series("score", score),
        pl.Series("AAI", score),
    )
    return pairs.rename({"query_seq": "qseqid", "target_seq": "sseqid"})


def run_proteome_similarity(
    all_prots: pl.DataFrame | pl.DataFrame | pl.LazyFrame,
    pairwise_aa: pl.DataFrame | pl.DataFrame | pl.LazyFrame,
    all_neigh: pl.DataFrame | pl.DataFrame | pl.LazyFrame | None = None,
    all_gff: pl.DataFrame | pl.DataFrame | pl.LazyFrame | None = None,
    outdir: str | None = None,
    pident_min: float = 30.0,
    mode: str = "all",
    subset_mode: str | None = None,
    win: int | None = None,
    win_mode: str = "bp",
    parallel: bool = False,
    num_threads: int | None = None,
) -> pl.DataFrame | tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """
    Compute wGRR, AAI (RBH), and/or vContact2-like hypergeometric scores.

    Parameters
    ----------
    all_prots : Polars DataFrame
    pairwise_aa : Polars DataFrame (BLAST outfmt 6 or normalized)
    all_neigh : unused (kept for API compatibility)
    outdir : optional path to write TSVs
    pident_min : identity threshold for hits (default 30.0)
    mode : {'wgrr','aai','hyper','all'}
    parallel : if True and mode='all', run metrics concurrently (3 threads)
    num_threads : optional override for parallel thread pool size (defaults to 3 when parallel=True)

    Returns
    -------
    Polars DataFrame or tuple of Polars DataFrame
    """
    mode = (mode or "all").lower()
    valid = {"wgrr", "aai", "hyper", "all"}
    if mode not in valid:
        raise ValueError(f"Invalid mode='{mode}'. Choose from {sorted(valid)}")

    all_prots_pl = _ensure_polars(all_prots)
    pairwise_aa_pl = _ensure_polars(pairwise_aa)
    all_neigh_pl = None
    all_gff_pl = None
    if all_neigh is not None:
        all_neigh_pl = _ensure_polars(all_neigh)
    if all_gff is not None:
        all_gff_pl = _ensure_polars(all_gff)

    try:
        if "unique_id" in set(all_prots_pl.columns):
            console.log(
                "run_proteome_similarity: using 'unique_id' from proteins table to compute neighborhood-wise similarities"
            )
            all_prots_pl = all_prots_pl.with_columns(
                pl.col("unique_id").cast(pl.Utf8).alias("target_nuc")
            )
        elif (
            all_neigh_pl is not None
            and not all_neigh_pl.is_empty()
            and "target_nuc" in set(all_prots_pl.columns)
        ):
            if {"temp_seqid", "seqid", "unique_id"}.intersection(all_neigh_pl.columns):
                an_pl = all_neigh_pl.select(
                    [c for c in ["temp_seqid", "seqid", "unique_id"] if c in all_neigh_pl.columns]
                )
                mapping = {}

                if "temp_seqid" in an_pl.columns:
                    temp_map = (
                        an_pl.select(
                            [
                                pl.col("temp_seqid").cast(pl.Utf8).alias("k"),
                                pl.col("unique_id").cast(pl.Utf8).alias("v"),
                            ]
                        )
                        .drop_nulls()
                        .unique(subset=["k", "v"])
                    )
                    mapping.update(dict(zip(temp_map["k"].to_list(), temp_map["v"].to_list())))

                if "seqid" in an_pl.columns:
                    seq_group = (
                        an_pl.filter(
                            pl.col("seqid").is_not_null() & pl.col("unique_id").is_not_null()
                        )
                        .group_by("seqid")
                        .agg(
                            [
                                pl.col("unique_id").n_unique().alias("n"),
                                pl.col("unique_id").first().cast(pl.Utf8).alias("uid_first"),
                            ]
                        )
                        .filter(pl.col("n") == 1)
                        .select(
                            [
                                pl.col("seqid").cast(pl.Utf8).alias("k"),
                                pl.col("uid_first").alias("v"),
                            ]
                        )
                    )
                    mapping.update(dict(zip(seq_group["k"].to_list(), seq_group["v"].to_list())))

                if mapping:
                    console.log(
                        f"run_proteome_similarity: mapping target_nuc -> unique_id for {len(mapping)} seq keys (unambiguous mapping)"
                    )
                    all_prots_pl = (
                        all_prots_pl.with_columns(
                            pl.col("target_nuc")
                            .cast(pl.Utf8)
                            .replace(mapping, default=None)
                            .alias("target_nuc_mapped")
                        )
                        .with_columns(
                            pl.when(pl.col("target_nuc_mapped").is_not_null())
                            .then(pl.col("target_nuc_mapped"))
                            .otherwise(pl.col("target_nuc"))
                            .alias("target_nuc")
                        )
                        .drop("target_nuc_mapped")
                    )
                    console.log(
                        f"run_proteome_similarity: after mapping, unique target_nuc values={all_prots_pl.select('target_nuc').n_unique()}"
                    )
    except Exception as e:
        console.log(f"run_proteome_similarity: neighborhood remapping skipped due to error: {e}")

    if subset_mode:
        if all_neigh_pl is None or all_gff_pl is None:
            raise ValueError("subset_mode requires both all_neigh and all_gff to be provided")
        keep_ids = _compute_subset_protein_ids(
            all_prots_pl, all_neigh_pl, all_gff_pl, subset_mode, win=win, win_mode=win_mode
        )
        if keep_ids:
            pid_col = _find_prot_id_col(all_prots_pl)
            all_prots_pl = all_prots_pl.filter(pl.col(pid_col).is_in(list(keep_ids)))
        else:
            all_prots_pl = all_prots_pl.filter(pl.lit(False))
        if not pairwise_aa_pl.is_empty():
            if {"qseqid", "sseqid"}.issubset(set(pairwise_aa_pl.columns)):
                pairwise_aa_pl = pairwise_aa_pl.filter(
                    (pl.col("qseqid").is_in(list(keep_ids)))
                    & (pl.col("sseqid").is_in(list(keep_ids)))
                )
            else:
                raise ValueError("pairwise_aa must contain qseqid/sseqid to apply subsetting")

    wgrr_df = aai_df = vcon_df = None

    if mode == "all" and parallel:
        max_workers = max(1, int(num_threads)) if num_threads is not None else 3
        with _fut.ThreadPoolExecutor(max_workers=max_workers) as ex:
            fut_w = ex.submit(compute_wgrr, pairwise_aa_pl, all_prots_pl, pident_min, True, True)
            fut_a = ex.submit(compute_aai_rbh, pairwise_aa_pl, all_prots_pl, pident_min, True)
            fut_h = ex.submit(compute_vcontact2_hypergeom, all_prots_pl, 0.2, 2, "BH")
            wgrr_df = fut_w.result()
            aai_df = fut_a.result()
            vcon_df = fut_h.result()
    else:
        if mode in {"wgrr", "all"}:
            wgrr_df = compute_wgrr(
                pairwise_aa_pl, all_prots_pl, pident_min=pident_min, symmetric=True
            )
        if mode in {"aai", "all"}:
            aai_df = compute_aai_rbh(pairwise_aa_pl, all_prots_pl, pident_min=pident_min)
        if mode in {"hyper", "all"}:
            vcon_df = compute_vcontact2_hypergeom(all_prots_pl, max_df_frac=0.2, min_shared=2)

    if outdir:
        outdir_path = Path(outdir)
        outdir_path.mkdir(parents=True, exist_ok=True)
        if wgrr_df is not None:
            wgrr_df.write_csv(outdir_path / "wgrr.tsv", separator="\t", include_header=False)
        if aai_df is not None:
            aai_df.write_csv(outdir_path / "aai.tsv", separator="\t", include_header=False)
        if vcon_df is not None:
            vcon_df.write_csv(
                outdir_path / "vcontact_hypergeom.tsv", separator="\t", include_header=False
            )

    if mode == "wgrr":
        return wgrr_df
    elif mode == "aai":
        return aai_df
    elif mode == "hyper":
        return vcon_df
    else:
        if wgrr_df is None:
            wgrr_df = pl.DataFrame(columns=["qseqid", "sseqid", "wGRR_sym", "AAI"])
        if aai_df is None:
            aai_df = pl.DataFrame(columns=["qseqid", "sseqid", "AAI", "n_RBH", "RBH_frac_min"])
        if vcon_df is None:
            vcon_df = pl.DataFrame(
                columns=[
                    "qseqid",
                    "sseqid",
                    "k",
                    "K_A",
                    "K_B",
                    "M",
                    "pval",
                    "p_adj",
                    "score",
                    "AAI",
                ]
            )
        return wgrr_df, aai_df, vcon_df
