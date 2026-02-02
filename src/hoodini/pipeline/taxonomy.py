import os
import subprocess
import warnings
from collections.abc import Iterable
from pathlib import Path

warnings.filterwarnings("ignore", category=UserWarning, module="UniProtMapper")
warnings.filterwarnings("ignore", category=UserWarning, module="numpy.core.getlimits")

import numpy as np  # noqa: E402
import polars as pl  # noqa: E402
from alphafetcher import AlphaFetcher  # noqa: E402
from ete3 import NCBITaxa  # noqa: E402
from scipy.cluster import hierarchy  # noqa: E402
from scipy.spatial.distance import pdist, squareform  # noqa: E402
from UniProtMapper import ProtMapper  # noqa: E402

from hoodini.utils.logging_utils import console, success  # noqa: E402
from hoodini.utils.seq_io import read_fasta, to_fasta  # noqa: E402


def parse_taxonomy_and_build_tree(
    records,
    all_gff,
    all_prots,
    all_neigh,
    output_dir,
    tree_mode,
    tree_file=None,
    num_threads=4,
    pairwise_aai=None,
    pairwise_ani=None,
    valid_uids=None,
    aai_mode: str | None = None,
    ani_mode: str | None = None,
    aai_subset_mode: str | None = None,
    nj_algorithm: str | None = None,
):
    """
    Parse taxonomic information and build phylogenetic tree.

    Expected Files:
    ---------------
    - records: DataFrame with taxid, organism, unique_id
    - all_prots: DataFrame with protein sequences and fam_cluster
    - all_neigh: DataFrame with neighborhood metadata
    - tree_file: Optional user-provided Newick tree file
    - {output}/target_prots.aln (if tree_mode == 'target_tree')
    - {output}/aai_matrix.tsv (if tree_mode == 'aai_tree')
    - {output}/ani_matrix.tsv (if tree_mode == 'ani_tree')

    Generated Files:
    ----------------
    - {output}/tree.nwk: Newick format phylogenetic tree

    Process:
    --------
    1. Enriches records with NCBI taxonomy (superkingdom, phylum, class, order, family, genus, species)
    2. Builds phylogenetic tree based on tree_mode:
            - 'target_tree': Uses target protein alignment
            - 'aai_tree': Uses average amino acid identity matrix
            - 'ani_tree': Uses average nucleotide identity matrix
            - 'user': Uses provided tree_file
            - 'taxonomy': Uses NCBI taxonomy hierarchy
    3. Creates dendrogram metadata (den_data) for visualization

    Returns:
    --------
    tuple: (tree_str: str, den_data: pl.DataFrame)
        - tree_str: Newick format tree string
        - den_data: DataFrame with leaf_id, taxonomy columns, and neighborhood coordinates
    """

    os.makedirs(output_dir, exist_ok=True)

    uid_map = {}
    if all_neigh is not None and len(all_neigh.columns) > 0:
        cols_available = set(all_neigh.columns)
        has_seq = "seqid" in cols_available
        has_temp = "temp_seqid" in cols_available
        if "unique_id" in cols_available and (has_seq or has_temp):
            for r in all_neigh.select(
                [c for c in ["unique_id", "seqid", "temp_seqid"] if c in cols_available]
            ).iter_rows(named=True):
                uid = r.get("unique_id")
                if uid is None:
                    continue
                for cand in [r.get("seqid"), r.get("temp_seqid")]:
                    if cand is None:
                        continue
                    s = str(cand)
                    uid_map[s] = uid
                    p = Path(s)
                    uid_map[p.name] = uid
                    uid_map[p.stem] = uid

    def _map_label(val):
        if val is None:
            return None
        s = str(val)
        return uid_map.get(s) or uid_map.get(Path(s).name) or uid_map.get(Path(s).stem) or s

    if tree_mode == "taxonomy":
        tree_str = _make_taxonomic_tree(records)
    elif tree_mode == "fast_nj":
        tree_str = _make_fast_phylo_tree(
            records, all_prots, output_dir, nj_algorithm, threads=num_threads
        )
    elif tree_mode == "aai_tree":
        chosen = pairwise_aai
        if isinstance(pairwise_aai, tuple) and len(pairwise_aai) >= 2:
            wgrr_df, aai_df = pairwise_aai[0], pairwise_aai[1]
            if aai_df is not None and not aai_df.is_empty():
                chosen = aai_df
            elif wgrr_df is not None:
                chosen = wgrr_df

        if chosen is None or (hasattr(chosen, "is_empty") and chosen.is_empty()):
            raise ValueError("pairwise_aai is empty")

        cols = set(chosen.columns)
        qcol = scol = pcol = None
        if {"qseqid", "sseqid"}.issubset(cols):
            qcol, scol = "qseqid", "sseqid"
            if "pident" in cols:
                pcol = "pident"
            elif "AAI" in cols:
                pcol = "AAI"
        elif {"qprot", "tprot"}.issubset(cols):
            qcol, scol = "qprot", "tprot"
            pcol = "pident" if "pident" in cols else ("AAI" if "AAI" in cols else None)
        elif {"A", "B"}.issubset(cols):
            qcol, scol = "A", "B"
            if "AAI" in cols:
                pcol = "AAI"
            elif "wGRR_sym" in cols:
                pcol = "wGRR_sym"

        if not all([qcol, scol, pcol]):
            raise ValueError(f"pairwise_aai is missing required columns; found {sorted(cols)}")

        chosen = chosen.with_columns(
            [
                pl.col(qcol).map_elements(_map_label).alias(qcol),
                pl.col(scol).map_elements(_map_label).alias(scol),
            ]
        )

        tree_str = aai_tree(
            chosen,
            qcol=qcol,
            scol=scol,
            pcol=pcol,
            valid_uids=valid_uids,
            algorithm=nj_algorithm,
            threads=num_threads,
            mode=aai_mode,
            subset_mode=aai_subset_mode,
        )
    elif tree_mode == "ani_tree":
        if pairwise_ani is None or (hasattr(pairwise_ani, "is_empty") and pairwise_ani.is_empty()):
            raise ValueError("pairwise_ani is empty")

        cols = set(pairwise_ani.columns)
        qcol = scol = pcol = None

        if {"qseqid", "sseqid"}.issubset(cols):
            qcol, scol = "qseqid", "sseqid"
            pcol = "pident" if "pident" in cols else ("ANI" if "ANI" in cols else None)
        elif {"Ref_name", "Query_name"}.issubset(cols):
            qcol, scol = "Query_name", "Ref_name"
            if "ANI" in cols:
                pcol = "ANI"
        elif {"A", "B"}.issubset(cols):
            qcol, scol = "A", "B"
            pcol = "ANI" if "ANI" in cols else None

        if not all([qcol, scol, pcol]):
            raise ValueError(f"pairwise_ani is missing required columns; found {sorted(cols)}")

        pairwise_ani = pairwise_ani.with_columns(
            [
                pl.col(qcol).map_elements(_map_label).alias(qcol),
                pl.col(scol).map_elements(_map_label).alias(scol),
            ]
        )

        tree_str = ani_tree(
            pairwise_ani,
            qcol=qcol,
            scol=scol,
            pcol=pcol,
            valid_uids=valid_uids,
            algorithm=nj_algorithm,
            threads=num_threads,
            mode=ani_mode,
        )
    elif tree_mode == "fast_ml":
        tree_str = _make_tree(records, all_prots, output_dir, num_threads)
    elif tree_mode == "use_input_tree":
        with open(tree_file) as f:
            tree_str = f.read()
    elif tree_mode == "foldmason_tree":
        tree_str = _make_foldmason_tree(records, all_prots, output_dir, num_threads)
    elif tree_mode == "neigh_similarity_tree":
        tree_str = _make_neigh_similarity_tree(all_prots, all_neigh)
    elif tree_mode == "neigh_phylo_tree":
        tree_str = _make_neigh_phylo_tree(records, all_prots, all_neigh, all_gff)
    else:
        raise ValueError(f"Unsupported tree mode: {tree_mode}")

    with open(f"{output_dir}/tree.nwk", "w") as out:
        out.write(tree_str)

    den_data = _build_leaf_metadata(records, all_neigh)
    success("Tree saved as Newick")
    return tree_str, den_data


def _build_leaf_metadata(records: pl.DataFrame, all_neigh: pl.DataFrame) -> pl.DataFrame:
    """Annotate leaf metadata with taxonomy using Polars only."""
    from hoodini.utils.logging_utils import info

    taxcols = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]

    records = records.with_columns(
        pl.col("taxid").fill_null("32644").cast(pl.Utf8).alias("taxid"),
        pl.col("unique_id").cast(pl.Utf8),
    )
    all_neigh = all_neigh.with_columns(pl.col("unique_id").cast(pl.Utf8))

    info("Loading NCBI taxonomy database (may take a few minutes on first run)...")
    ncbi = NCBITaxa()

    # Get unique taxids as integers (filter out non-numeric)
    unique_taxids = records.select("taxid").unique().to_series().to_list()
    taxids_int = []
    taxid_str_to_int = {}
    for t in unique_taxids:
        try:
            tid = int(str(t))
            taxids_int.append(tid)
            taxid_str_to_int[str(t)] = tid
        except (ValueError, TypeError):
            pass
    info(f"Looking up taxonomy for {len(taxids_int)} unique taxids...")

    # Batch lookup all lineages at once (much faster than individual calls)
    lineage_map = ncbi.get_lineage_translator(taxids_int) if taxids_int else {}

    # Collect all taxids from lineages for batch rank/name lookup
    all_lineage_taxids = set()
    for lineage in lineage_map.values():
        if lineage:
            all_lineage_taxids.update(lineage)

    # Batch lookup ranks and names for ALL taxids at once
    all_ranks = ncbi.get_rank(list(all_lineage_taxids)) if all_lineage_taxids else {}
    all_names = ncbi.get_taxid_translator(list(all_lineage_taxids)) if all_lineage_taxids else {}

    # Build taxonomy rows using cached data
    tax_rows = []
    for taxid in unique_taxids:
        row = {"taxid": str(taxid)}
        tid = taxid_str_to_int.get(str(taxid))
        lineage = lineage_map.get(tid, []) if tid else []
        for ltid in lineage:
            rank_name = all_ranks.get(ltid, "")
            if rank_name in taxcols:
                row[rank_name] = all_names.get(ltid, "")
        for c in taxcols:
            row.setdefault(c, None)
        tax_rows.append(row)

    taxdf = (
        pl.DataFrame(tax_rows, schema={"taxid": pl.Utf8, **dict.fromkeys(taxcols, pl.Utf8)})
        if tax_rows
        else pl.DataFrame([{"taxid": "32644", **dict.fromkeys(taxcols)}])
    )

    records = records.join(taxdf, on="taxid", how="left")

    missing_cols = [c for c in taxcols if c not in records.columns]
    if missing_cols:
        records = records.with_columns([pl.lit(None).alias(c) for c in missing_cols])

    records = records.with_columns(
        [
            pl.when(pl.col(c).is_null()).then(pl.lit("unclassified")).otherwise(pl.col(c)).alias(c)
            for c in taxcols
        ]
    )

    # Include dive columns (BacDive/PhageDive) if present
    dive_cols = ["dive_id", "collection_id", "dive_type"]
    base_cols = ["unique_id", "og_index"] + taxcols
    extra_cols = [c for c in dive_cols if c in records.columns]
    
    # Include user-provided extra columns from inputsheet
    # These are columns that are not in the reserved set (defined in validation.py)
    from hoodini.utils.validation import RESERVED_COLUMNS
    user_extra_cols = [c for c in records.columns 
                       if c not in RESERVED_COLUMNS 
                       and c not in base_cols 
                       and c not in extra_cols
                       and not c.startswith("_")]  # Exclude internal columns

    # Filter out failed records - only include successful ones in tree metadata
    valid_records = records.filter(pl.col("failed").is_null())
    den_data = valid_records.select(base_cols + extra_cols + user_extra_cols)

    # Select neighborhood columns including seqid and target_prot (as align_gene)
    neigh_cols = ["unique_id", "start_win", "end_win", "strand_win", "start_target", "end_target"]
    if "seqid" in all_neigh.columns:
        neigh_cols.append("seqid")
    if "target_prot" in all_neigh.columns:
        neigh_cols.append("target_prot")

    neigh_data = all_neigh.select(neigh_cols)
    # Rename target_prot to align_gene for visualization consistency
    if "target_prot" in neigh_data.columns:
        neigh_data = neigh_data.rename({"target_prot": "align_gene"})

    den_data = den_data.join(neigh_data, on="unique_id", how="left")
    return den_data


def _make_tree(records, all_prots, output_dir, threads):
    from hoodini.utils.logging_utils import info, warn

    valid = records.filter(pl.col("failed").is_null())
    prots = (
        valid.select(["protein_id", "unique_id"])
        .drop_nulls(subset=["protein_id", "unique_id"])
        .unique(subset=["unique_id"])
    )
    faa = all_prots.join(prots, left_on="id", right_on="protein_id", how="inner")

    if "unique_id" not in faa.columns:
        id_to_uid = valid.select(["protein_id", "unique_id"]).drop_nulls()
        faa = faa.join(id_to_uid, left_on="id", right_on="protein_id", how="left")

    # ensure one sequence per unique_id (neighborhood) and preserve original mapping
    faa = faa.unique(subset=["unique_id"])
    info(f"Building tree from {faa.height} target protein sequences...")

    # Handle edge case: no target proteins available
    if faa.height == 0:
        # Check if there are ANY valid records for fallback
        valid_count = records.filter(pl.col("failed").is_null()).height
        if valid_count == 0:
            from hoodini.utils.logging_utils import error

            error("No valid records found. All input sequences failed during processing.")
            raise SystemExit(1)
        warn("No target proteins found for tree building. Falling back to taxonomy tree.")
        return _make_taxonomic_tree(records)

    # Handle edge case: only 1 sequence (can't build a tree)
    if faa.height == 1:
        uid = faa["unique_id"][0]
        warn(f"Only 1 target protein found ({uid}). Creating single-leaf tree.")
        return f"({uid}:0.0);"

    # Handle edge case: only 2 sequences (FAMSA works but VeryFastTree may have issues)
    if faa.height == 2:
        uids = faa["unique_id"].to_list()
        warn("Only 2 target proteins found. Creating simple 2-leaf tree.")
        return f"({uids[0]}:0.1,{uids[1]}:0.1);"

    to_fasta(faa, "unique_id", "sequence", f"{output_dir}/target_prots.fasta")

    info("Running FAMSA alignment...")
    subprocess.run(
        [
            "famsa",
            "-t",
            str(threads),
            f"{output_dir}/target_prots.fasta",
            f"{output_dir}/target_prots.aln",
            "-remove-rare-columns",
            "0.10",
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    info("Running VeryFastTree...")
    result = subprocess.run(
        ["VeryFastTree", "-threads", str(threads), f"{output_dir}/target_prots.aln"],
        capture_output=True,
        text=True,
    )
    return result.stdout


def _make_fast_phylo_tree(records, all_prots, output_dir, nj_algorithm, threads: int = 1):
    from hoodini.pipeline.helpers.decenttree_builder import run_decenttree_from_matrix

    prots = (
        records.filter(pl.col("failed").is_null())
        .select(["protein_id", "unique_id"])
        .drop_nulls(subset=["protein_id", "unique_id"])
        .unique(subset=["unique_id"])
    )
    faa = all_prots.join(prots, left_on="id", right_on="protein_id", how="inner")

    if "unique_id" not in faa.columns:
        id_to_uid = records.select(["protein_id", "unique_id"]).drop_nulls()
        faa = faa.join(id_to_uid, left_on="id", right_on="protein_id", how="left")

    # one sequence per unique_id; header uses unique_id to match original record
    faa = faa.unique(subset=["unique_id"])
    to_fasta(faa, "unique_id", "sequence", f"{output_dir}/target_prots.fasta")
    subprocess.run(
        [
            "famsa",
            "-dist_export",
            "-square_matrix",
            f"{output_dir}/target_prots.fasta",
            f"{output_dir}/distance_matrix.csv",
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    dispd = pl.read_csv(f"{output_dir}/distance_matrix.csv")
    first_col = dispd.columns[0]
    dispd = dispd.rename({first_col: "id"})
    dispd = dispd.select(["id"] + [c for c in dispd.columns if c != "id"])
    newick = run_decenttree_from_matrix(
        dispd, algorithm=nj_algorithm, threads=max(1, int(threads or 1))
    )
    return newick


def _make_taxonomic_tree(records):
    from hoodini.utils.logging_utils import error, warn

    valid = records.filter(pl.col("failed").is_null()).with_row_count("_idx")
    uids = valid["unique_id"].to_list()
    n = len(uids)

    # Handle edge cases
    if n == 0:
        error("No valid records found for tree building. All input sequences failed.")
        raise SystemExit(1)
    if n == 1:
        warn(f"Only 1 valid record for taxonomy tree ({uids[0]}). Creating single-leaf tree.")
        return f"({uids[0]}:0.0);"
    if n == 2:
        warn("Only 2 valid records for taxonomy tree. Creating simple 2-leaf tree.")
        return f"({uids[0]}:0.1,{uids[1]}:0.1);"

    # Filter out records without taxid and warn
    valid_with_taxid = valid.filter(pl.col("taxid").is_not_null())
    missing_taxid_count = n - valid_with_taxid.height
    if missing_taxid_count > 0:
        warn(f"{missing_taxid_count} records have no taxid and will use a default taxid for tree building.")
    
    # Use a default taxid (1 = root) for records without taxid
    DEFAULT_TAXID = 1
    valid = valid.with_columns(
        pl.col("taxid").fill_null(DEFAULT_TAXID).alias("taxid")
    )

    taxids = valid.select("taxid").unique().to_series().to_list()
    distances = calculate_taxid_distances(taxids, update_db=False)
    mat = np.zeros((n, n), dtype=float)
    for a in range(n):
        for b in range(a, n):
            if a == b:
                mat[a, b] = 0.0
            else:
                taxid_i = int(valid[a, "taxid"])
                taxid_j = int(valid[b, "taxid"])
                d = (
                    0
                    if taxid_i == taxid_j
                    else distances.get((taxid_i, taxid_j), distances.get((taxid_j, taxid_i), 1e6))
                )
                mat[a, b] = mat[b, a] = d
    linkage = hierarchy.linkage(squareform(mat), method="single")
    return _linkage_to_newick(linkage, uids)


def _make_neigh_similarity_tree(all_prot, all_neigh=None):
    pa = all_prot.filter(pl.col("id") != pl.col("target_prot"))
    counts = (
        pa.group_by(["target_prot", "fam_cluster"])
        .agg(pl.len().alias("count"))
        .with_columns(pl.when(pl.col("count") > 0).then(1).otherwise(0).alias("presence"))
    )
    mat = counts.pivot(
        index="target_prot", columns="fam_cluster", values="presence", aggregate_function="max"
    ).fill_null(0)
    binmat = mat.drop("target_prot")
    dist = pdist(binmat.to_numpy(), metric="jaccard")
    linkage = hierarchy.linkage(dist, method="single")

    # Map target_prot to unique_id if all_neigh is available
    labels = mat["target_prot"].to_list()
    if (
        all_neigh is not None
        and "target_prot" in all_prot.columns
        and "unique_id" in all_prot.columns
    ):
        prot_to_uid = dict(all_prot.select(["target_prot", "unique_id"]).unique().iter_rows())
        labels = [prot_to_uid.get(tp, tp) for tp in labels]

    return _linkage_to_newick(linkage, labels)


def _make_neigh_phylo_tree(records, all_prot, all_neigh=None, all_gff=None):
    # Need to join with all_gff to get positions and calculate relative positions
    if all_gff is None or all_gff.height == 0:
        # Fallback to simple similarity tree if we don't have position info
        return _make_neigh_similarity_tree(all_prot, all_neigh)

    # Join all_prot with all_gff to get start/end positions
    # Assuming all_gff has 'id' that matches all_prot 'id', and 'start'/'end' columns
    if (
        "id" not in all_gff.columns
        or "start" not in all_gff.columns
        or "end" not in all_gff.columns
    ):
        return _make_neigh_similarity_tree(all_prot, all_neigh)

    gff_pos = all_gff.select(["id", "start", "end", "unique_id"]).unique()
    prot_with_pos = all_prot.join(gff_pos, on="id", how="left", suffix="_gff")

    # Calculate target positions for each unique_id
    target_positions = (
        prot_with_pos.filter(pl.col("id") == pl.col("target_prot"))
        .select(["unique_id", "start", "end"])
        .rename({"start": "target_start", "end": "target_end"})
        .unique()
    )

    # Join to get target positions and calculate relative positions
    prot_with_pos = prot_with_pos.join(target_positions, on="unique_id", how="left").with_columns(
        [
            (pl.col("start") - pl.col("target_start")).alias("rel_start"),
            (pl.col("end") - pl.col("target_end")).alias("rel_end"),
        ]
    )

    pa = prot_with_pos.filter(pl.col("id") != pl.col("target_prot")).with_columns(
        ((pl.col("rel_start") + pl.col("rel_end")) / 2).alias("rel_pos")
    )
    weights = pa.with_columns((1.0 / (1 + (pl.col("rel_pos")).abs())).alias("w"))
    mat = (
        weights.group_by(["target_prot", "fam_cluster"])
        .agg(pl.col("w").sum().alias("w"))
        .pivot(index="target_prot", columns="fam_cluster", values="w", aggregate_function="sum")
        .fill_null(0.0)
    )
    feature_cols = [c for c in mat.columns if c != "target_prot"]
    mat = mat.with_columns(pl.sum_horizontal(pl.col(feature_cols)).alias("row_sum"))
    mat = mat.with_columns(
        [
            (pl.col(c) / pl.when(pl.col("row_sum") > 0).then(pl.col("row_sum")).otherwise(1)).alias(
                c
            )
            for c in feature_cols
        ]
    ).drop("row_sum")

    norm_vals = mat.select(feature_cols).to_numpy()
    dist = pdist(norm_vals, metric="cosine")
    linkage = hierarchy.linkage(dist, method="single")

    # Map target_prot to unique_id if all_neigh is available
    labels = mat["target_prot"].to_list()
    if (
        all_neigh is not None
        and "target_prot" in all_prot.columns
        and "unique_id" in all_prot.columns
    ):
        prot_to_uid = dict(all_prot.select(["target_prot", "unique_id"]).unique().iter_rows())
        labels = [prot_to_uid.get(tp, tp) for tp in labels]

    return _linkage_to_newick(linkage, labels)


def _pairwise_to_matrix(
    pairwise_df: pl.DataFrame,
    ids: Iterable[str] | None = None,
    qcol: str = "qseqid",
    scol: str = "sseqid",
    valcol: str = "pident",
    value_is_identity: bool = True,
) -> pl.DataFrame:
    """Convert pairwise table to full square distance matrix DataFrame.

    - pairwise_df: rows with qcol, scol, valcol
    - ids: iterable of identifiers that must be present (will be included even if missing)
    - value_is_identity: if True, converts identity -> distance as (100 - val)

    Missing pairs are left as NaN (caller will decide fill strategy).
    """
    if pairwise_df is None or pairwise_df.height == 0:
        base_ids = []
    else:
        base_ids = [str(x) for x in pl.unique(pairwise_df[[qcol, scol]].to_numpy().ravel())]

    if ids is not None:
        ids_list = [str(x) for x in ids]
        extras = [x for x in base_ids if x not in ids_list]
        all_ids = ids_list + extras
    else:
        all_ids = sorted(set(base_ids))

    n = len(all_ids)
    mat_np = np.full((n, n), np.nan, dtype=float)
    np.fill_diagonal(mat_np, 0.0)

    if pairwise_df is not None and pairwise_df.height > 0:
        for row in pairwise_df.iter_rows(named=True):
            try:
                q = str(row[qcol])
                s = str(row[scol])
                v = float(row[valcol])
            except Exception:
                continue
            d = 100.0 - v if value_is_identity else v
            if q not in all_ids:
                all_ids.append(q)
                mat_np = np.pad(mat_np, ((0, 1), (0, 1)), constant_values=np.nan)
                mat_np[-1, -1] = 0.0
                n = len(all_ids)
            if s not in all_ids:
                all_ids.append(s)
                mat_np = np.pad(mat_np, ((0, 1), (0, 1)), constant_values=np.nan)
                mat_np[-1, -1] = 0.0
                n = len(all_ids)
            qi = all_ids.index(q)
            si = all_ids.index(s)
            mat_np[qi, si] = d
            mat_np[si, qi] = d

    rows = []
    for i, rid in enumerate(all_ids):
        row = {"id": rid}
        for j, cid in enumerate(all_ids):
            row[cid] = mat_np[i, j]
        rows.append(row)

    return pl.DataFrame(rows)


def aai_tree(
    pairwise_aai: pl.DataFrame,
    valid_uids: Iterable[str] | None = None,
    qcol: str = "qseqid",
    scol: str = "sseqid",
    pcol: str = "pident",
    algorithm: str = "nj",
    threads: int = 1,
    mode: str | None = None,
    subset_mode: str | None = None,
) -> str:
    """Build a tree from AAI pairwise table (pident) using DecentTree.

    Missing pairs (and ids not present in the table but in valid_uids) are
    filled with (max_observed_distance + 2*std_observed).
    """
    from hoodini.pipeline.helpers.decenttree_builder import run_decenttree_from_table

    if pairwise_aai is None or pairwise_aai.height == 0:
        raise ValueError("pairwise_aai is empty")

    df = (
        pairwise_aai.rename({qcol: "A", scol: "B", pcol: "AAI"})
        .with_columns(pl.col("AAI").cast(pl.Float64, strict=False))
        .drop_nulls("AAI")
    )
    if df.height == 0:
        raise ValueError("No valid AAI values available to build tree")

    if (df["AAI"].max() or 0) <= 1.0:
        df = df.with_columns((pl.col("AAI") * 100.0).alias("AAI"))

    if mode and mode.lower() == "hyper":
        raise ValueError("'hyper' mode is not supported for AAI trees")

    return run_decenttree_from_table(
        df,
        qcol="A",
        tcol="B",
        dcol="distance",
        algorithm=algorithm,
        threads=threads,
        fill_missing="max+2std",
        ids=valid_uids,
    )


def ani_tree(
    pairwise_ani: pl.DataFrame,
    valid_uids: Iterable[str] | None = None,
    qcol: str = "A",
    scol: str = "B",
    pcol: str = "ANI",
    algorithm: str = "nj",
    threads: int = 1,
    mode: str | None = None,
) -> str:
    """Build a tree from ANI pairwise table using DecentTree.

    If the ANI table stores percent identity, distances are computed as
    (100 - ANI). Missing pairs are filled with max+2std as for AAI.
    """
    from hoodini.pipeline.helpers.decenttree_builder import run_decenttree_from_table

    if pairwise_ani is None or pairwise_ani.height == 0:
        raise ValueError("pairwise_ani is empty")

    df = (
        pairwise_ani[[qcol, scol, pcol]]
        .rename({qcol: "A", scol: "B", pcol: "ANI"})
        .with_columns(pl.col("ANI").cast(pl.Float64, strict=False))
        .drop_nulls("ANI")
    )
    if df.height == 0:
        raise ValueError("No valid ANI values available to build tree")

    if (df["ANI"].max() or 0) <= 1.0:
        df = df.with_columns((pl.col("ANI") * 100.0).alias("ANI"))

    if valid_uids is not None:
        uid_list = [str(x) for x in valid_uids]
        uid_set = set(uid_list)

        def _map_val(x: str) -> str:
            s = str(x)
            if s in uid_set:
                return s
            b = Path(s).name
            if b in uid_set:
                return b
            for uid in uid_list:
                if uid in s:
                    return uid
            return s

        df = df.with_columns(
            [
                pl.col(c).map_elements(_map_val, return_dtype=pl.Utf8)
                for c in ("A", "B")
                if c in df.columns
            ]
        )

    return run_decenttree_from_table(
        df,
        qcol="A",
        tcol="B",
        dcol="ANI",
        algorithm=algorithm,
        threads=threads,
        fill_missing="max+2std",
        ids=valid_uids,
    )


def _make_foldmason_tree(records, all_prot, output_dir, threads):
    """Build tree using FoldMason structural alignment.
    
    Uses AlphaFold structures for proteins that can be mapped to UniProt IDs.
    Proteins without structures are added via MAFFT --add.
    
    Priority for UniProt ID resolution:
    1. Use existing uniprot_id from input (preserved from inputsheet)
    2. Map NCBI protein_id → UniProt using ProtMapper
    """
    from hoodini.utils.logging_utils import info, warn
    
    # Normalize inputs to Polars for consistency through the pipeline
    records = records.collect() if isinstance(records, pl.LazyFrame) else records
    records = records if isinstance(records, pl.DataFrame) else pl.from_pandas(records)
    all_prot_pl = all_prot.collect() if isinstance(all_prot, pl.LazyFrame) else all_prot
    all_prot_pl = all_prot_pl if isinstance(all_prot_pl, pl.DataFrame) else pl.from_pandas(all_prot)

    valid = records.filter(pl.col("failed").is_null())
    
    # Build protein_id -> uniprot_id mapping
    # Priority 1: Use existing uniprot_id from records (from original input)
    has_uniprot = "uniprot_id" in valid.columns
    existing_uniprot = {}
    if has_uniprot:
        uniprot_rows = valid.filter(
            pl.col("uniprot_id").is_not_null() & (pl.col("uniprot_id") != "")
        ).select(["protein_id", "uniprot_id"]).unique()
        for row in uniprot_rows.iter_rows(named=True):
            if row["protein_id"] and row["uniprot_id"]:
                existing_uniprot[row["protein_id"]] = row["uniprot_id"]
    
    info(f"Found {len(existing_uniprot)} proteins with existing UniProt IDs from input")
    
    # Get all target protein IDs
    targets = valid.select("protein_id").drop_nulls().unique().to_series().to_list()
    
    # Priority 2: Map remaining proteins without uniprot_id
    needs_mapping = [t for t in targets if t not in existing_uniprot]
    
    import pandas as pd
    
    mapped_from_api = {}
    no_map = []
    
    if needs_mapping:
        info(f"Mapping {len(needs_mapping)} proteins to UniProt via API...")
        mapper = ProtMapper()
        try:
            mapped_pd, failed_ids = mapper.get(
                ids=needs_mapping, from_db="EMBL-GenBank-DDBJ_CDS", to_db="UniProtKB"
            )
            no_map = failed_ids if failed_ids else []
            
            if mapped_pd is not None and not mapped_pd.empty:
                for _, row in mapped_pd.iterrows():
                    mapped_from_api[row["From"]] = row["Entry"]
        except Exception as e:
            warn(f"UniProt mapping failed: {e}")
            no_map = needs_mapping
    
    # Combine mappings (existing takes priority)
    protein_to_uniprot = {**mapped_from_api, **existing_uniprot}  # existing overwrites
    
    # Proteins that have UniProt IDs (either from input or mapped)
    proteins_with_uniprot = [t for t in targets if t in protein_to_uniprot]
    uniprot_ids = [protein_to_uniprot[t] for t in proteins_with_uniprot]
    
    # Proteins without UniProt mapping
    no_uniprot = [t for t in targets if t not in protein_to_uniprot]
    
    info(f"Proteins with UniProt: {len(proteins_with_uniprot)}, without: {len(no_uniprot)}")
    
    # Check if we have any UniProt IDs to work with
    if not uniprot_ids:
        console.print(
            "[bold red]Error: None of the protein IDs could be mapped to UniProt IDs.[/bold red]"
        )
        console.print(
            "[yellow]Foldmason tree requires AlphaFold structures, which need UniProt IDs.[/yellow]"
        )
        console.print(
            "[yellow]Falling back to standard alignment tree (FAMSA + VeryFastTree).[/yellow]"
        )
        return _make_tree(records, all_prot, output_dir, threads)

    # Fetch AlphaFold structures
    fetcher = AlphaFetcher(base_savedir=f"{output_dir}/struct")
    unique_uniprot = list(set(uniprot_ids))
    
    info(f"Fetching AlphaFold structures for {len(unique_uniprot)} UniProt IDs...")
    fetcher.add_proteins(unique_uniprot)
    fetcher.fetch_metadata(multithread=True, workers=threads)
    
    # Track which proteins don't have PDB structures
    failed_uniprot = set(fetcher.failed_ids) if fetcher.failed_ids else set()
    no_pdb_proteins = [t for t in targets if protein_to_uniprot.get(t) in failed_uniprot]
    no_pdb_proteins.extend(no_uniprot)
    
    # Download structures for successful ones
    successful_uniprot = [u for u in unique_uniprot if u not in failed_uniprot]
    if successful_uniprot:
        info(f"Downloading {len(successful_uniprot)} PDB structures...")
        fetcher.download_all_files(pdb=True, cif=False, multithread=True, workers=threads)

    # Check if foldmason is available
    import shutil
    if shutil.which("foldmason") is None:
        warn("foldmason not found in PATH. Falling back to FAMSA + VeryFastTree.")
        return _make_tree(records, all_prot, output_dir, threads)

    # Run foldmason
    info("Running FoldMason structural alignment...")
    try:
        subprocess.run(
            [
                "foldmason",
                "easy-msa",
                f"{output_dir}/struct/pdb_files",
                f"{output_dir}/foldmason",
                f"{output_dir}/temp",
            ],
            check=True,
        )
    except subprocess.CalledProcessError as e:
        warn(f"FoldMason failed: {e}. Falling back to FAMSA + VeryFastTree.")
        return _make_tree(records, all_prot, output_dir, threads)

    # Handle proteins without PDB structures (add via MAFFT)
    if no_pdb_proteins:
        info(f"Adding {len(no_pdb_proteins)} proteins without structures via MAFFT...")
        missing_df = (
            all_prot_pl.filter(pl.col("id").is_in(no_pdb_proteins))
            .select(["id", "sequence"])
            .unique(subset=["id"])
        )
        if missing_df.height > 0:
            missing_df.to_fasta("id", "sequence", f"{output_dir}/no_pdb.fasta")
            cmd = [
                "mafft",
                "--keeplength",
                "--add",
                f"{output_dir}/no_pdb.fasta",
                "--reorder",
                f"{output_dir}/foldmason_aa.fa",
            ]
            with open(f"{output_dir}/target_prots.aln", "w") as out:
                subprocess.run(cmd, check=True, stdout=out)
        else:
            os.rename(f"{output_dir}/foldmason_aa.fa", f"{output_dir}/target_prots.aln")
    else:
        os.rename(f"{output_dir}/foldmason_aa.fa", f"{output_dir}/target_prots.aln")

    # Build mapping: UniProt ID -> protein_id -> unique_id
    # The alignment has UniProt IDs (from foldmason) or protein IDs (from mafft --add)
    uniprot_to_protein = {v: k for k, v in protein_to_uniprot.items()}
    
    # Build protein_id -> unique_id mapping from records
    protein_to_uid = {}
    for row in valid.select(["protein_id", "unique_id"]).drop_nulls().unique().iter_rows(named=True):
        if row["protein_id"] and row["unique_id"]:
            protein_to_uid[row["protein_id"]] = str(row["unique_id"])
    
    aln_df = read_fasta(f"{output_dir}/target_prots.aln")
    
    # Map alignment IDs to unique_id:
    # 1. UniProt ID -> protein_id -> unique_id
    # 2. protein_id -> unique_id (if already a protein_id)
    def map_to_uid(x):
        # First try as UniProt -> protein_id
        prot_id = uniprot_to_protein.get(x, x)
        # Then protein_id -> unique_id
        return protein_to_uid.get(prot_id, prot_id)
    
    aln_df = aln_df.with_columns(
        pl.col("id").map_elements(map_to_uid, return_dtype=pl.Utf8).alias("unique_id")
    )
    aln_df = aln_df.unique(subset=["unique_id"])
    aln_df.to_fasta("unique_id", "sequence", f"{output_dir}/target_prots.aln")

    info("Building tree with VeryFastTree...")
    result = subprocess.run(
        ["VeryFastTree", f"{output_dir}/target_prots.aln"], capture_output=True, text=True
    )
    return result.stdout


def _linkage_to_newick(linkage, labels):
    """Convert scipy linkage to Newick format string using iterative approach.

    Avoids recursion depth limits by using explicit stack.
    Includes branch lengths from the linkage matrix.
    """
    tree = hierarchy.to_tree(linkage)

    # Iterative post-order traversal
    stack = [(tree, False)]
    result_stack = []

    while stack:
        node, visited = stack.pop()
        if node.is_leaf():
            # Leaf nodes: just label (branch length added by parent)
            result_stack.append((str(labels[node.id]), 0.0))
        elif visited:
            # Post-order: both children have been processed
            right_str, right_dist = result_stack.pop()
            left_str, left_dist = result_stack.pop()

            # Calculate branch lengths from this node to children
            left_branch = node.dist - left_dist if left_dist < node.dist else 0.0
            right_branch = node.dist - right_dist if right_dist < node.dist else 0.0

            # Format with branch lengths
            subtree = f"({left_str}:{left_branch:.6f},{right_str}:{right_branch:.6f})"
            result_stack.append((subtree, node.dist))
        else:
            # Pre-order: mark for post-processing and push children
            stack.append((node, True))
            stack.append((node.right, False))
            stack.append((node.left, False))

    return result_stack[0][0] + ";"


def calculate_taxid_distances(taxids, update_db=False):
    """Calculate pairwise taxonomic distances between taxids using ete3.

    Distance is computed as the number of steps to the lowest common ancestor (LCA)
    from both taxa combined.
    """
    import itertools

    ncbi = NCBITaxa()
    taxids_int = [int(taxid) for taxid in taxids]

    # Build lineage cache for each taxid
    lineage_cache = {}
    for taxid in taxids_int:
        try:
            lineage = ncbi.get_lineage(taxid)
            # get_lineage returns lineage from root to taxid, reverse it
            lineage_cache[taxid] = list(reversed(lineage)) if lineage else [taxid, 1]
        except Exception:
            # Fallback to unclassified if taxid not found
            lineage_cache[taxid] = [taxid, 1]  # 1 is root

    distances = {}
    for taxid1, taxid2 in itertools.combinations(taxids_int, 2):
        lineage1 = set(lineage_cache[taxid1])
        lineage2 = set(lineage_cache[taxid2])

        # Find lowest common ancestor (first shared taxid in lineages)
        common = lineage1 & lineage2
        if not common:
            # No common ancestor found, use max distance
            distance = len(lineage_cache[taxid1]) + len(lineage_cache[taxid2])
        else:
            # Distance = steps from taxid1 to LCA + steps from taxid2 to LCA
            # LCA is the one with the highest index (closest to the taxa)
            lca = None
            for t in lineage_cache[taxid1]:
                if t in common:
                    lca = t
                    break

            dist1 = (
                lineage_cache[taxid1].index(lca)
                if lca in lineage_cache[taxid1]
                else len(lineage_cache[taxid1])
            )
            dist2 = (
                lineage_cache[taxid2].index(lca)
                if lca in lineage_cache[taxid2]
                else len(lineage_cache[taxid2])
            )
            distance = dist1 + dist2

        distances[(taxid1, taxid2)] = distance

    return distances
