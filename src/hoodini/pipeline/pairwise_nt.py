import contextlib
import math
import os
import subprocess
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import mappy as mp
import polars as pl
from Bio import SeqIO
from rich.progress import (
    BarColumn,
    Progress,
    SpinnerColumn,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)

from hoodini.utils.logging_utils import info, warn
from hoodini.utils.seq_io import to_fasta  # noqa: F401 - registers pl.DataFrame.to_fasta


def _is_na(x):
    try:
        return x is None or (isinstance(x, float) and math.isnan(x))
    except Exception:
        return x is None


def _fasta_ids(fasta_path: str):
    ids = [rec.id for rec in SeqIO.parse(fasta_path, "fasta")]
    if not ids:
        raise ValueError("No sequences found.")
    dup = [k for k, v in Counter(ids).items() if v > 1]
    if dup:
        raise ValueError(f"Duplicate FASTA IDs: {dup[:10]}{'...' if len(dup) > 10 else ''}")
    return ids


def _fasta_lengths(path: str):
    idx = SeqIO.index(path, "fasta")
    try:
        return {k: len(idx[k].seq) for k in idx}
    finally:
        idx.close()


def _write_intergenic_fasta(
    all_gff, fasta_path: str, out_fasta: str, all_neigh: pl.DataFrame | None = None
):
    """
    Create a FASTA of intergenic regions from neighborhoods.fasta using a GFF (path or DataFrame).
    Returns a DataFrame with columns ['temp_seqid','seqid','start_win'].
    """

    if all_gff is None:
        raise ValueError("all_gff is required for intergenic FASTA generation")
    if isinstance(all_gff, pl.DataFrame):
        gff = all_gff.clone()
    else:
        cols = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
        gff = pl.read_csv(
            str(all_gff), separator="\t", has_header=False, new_columns=cols, comment="#"
        )

    coding_types = ["CDS", "gene"]
    gff_coding = gff.filter(pl.col("type").is_in(coding_types))

    temp_to_seqid = {}
    start_map_init = {}
    if all_neigh is not None and all_neigh.height > 0:
        for nr in all_neigh.iter_rows(named=True):
            temp = (
                str(nr.get("temp_seqid"))
                if ("temp_seqid" in all_neigh.columns and nr.get("temp_seqid") is not None)
                else None
            )
            seqid = (
                str(nr.get("seqid"))
                if ("seqid" in all_neigh.columns and nr.get("seqid") is not None)
                else temp
            )
            if temp:
                temp_to_seqid[temp] = seqid
            if seqid and "start_win" in nr and nr.get("start_win") is not None:
                start_map_init[seqid] = int(nr.get("start_win"))

    meta_rows = []
    with open(out_fasta, "w") as outfh:
        for rec in SeqIO.parse(str(fasta_path), "fasta"):
            hdr = rec.id
            seq = str(rec.seq)
            seqlen = len(seq)
            canonical = temp_to_seqid.get(hdr, hdr)

            cand = gff_coding.filter(
                pl.col("seqid").is_in([canonical, Path(canonical).name, Path(canonical).stem])
            )
            coding_intervals = []

            start_win = start_map_init.get(canonical)
            if start_win is not None:
                for r in cand.iter_rows(named=True):
                    try:
                        s_global = int(r.get("start"))
                        e_global = int(r.get("end"))
                    except Exception:
                        continue
                    s = s_global - int(start_win) + 1
                    e = e_global - int(start_win) + 1
                    s = max(1, min(s, seqlen))
                    e = max(1, min(e, seqlen))
                    if e >= s:
                        coding_intervals.append((s, e))
            else:
                warn(f"No start_win for {canonical}; treating GFF coords as neighborhood-local")
                for r in cand.iter_rows(named=True):
                    try:
                        s = int(r.get("start"))
                        e = int(r.get("end"))
                    except Exception:
                        continue
                    s = max(1, min(s, seqlen))
                    e = max(1, min(e, seqlen))
                    if e >= s:
                        coding_intervals.append((s, e))

            coding_intervals.sort()
            merged = []
            for iv in coding_intervals:
                if not merged:
                    merged.append(list(iv))
                elif iv[0] <= merged[-1][1] + 1:
                    merged[-1][1] = max(merged[-1][1], iv[1])
                else:
                    merged.append([iv[0], iv[1]])

            intergenic = []
            pos = 1
            for s, e in merged:
                if pos < s:
                    intergenic.append((pos, s - 1))
                pos = e + 1
            if pos <= seqlen:
                intergenic.append((pos, seqlen))
            if not merged:
                intergenic = [(1, seqlen)]

            for s, e in intergenic:
                if e < s:
                    continue
                sub = seq[s - 1 : e]
                temp_id = f"{hdr}__inter_{s}_{e}"
                outfh.write(f">{temp_id}\n")
                for i in range(0, len(sub), 80):
                    outfh.write(sub[i : i + 80] + "\n")
                start_win_val = start_map_init.get(canonical, 0)
                meta_rows.append(
                    {
                        "temp_seqid": temp_id,
                        "seqid": canonical,
                        "start_win": int(s) + int(start_win_val) - 1,
                    }
                )

    return pl.DataFrame(meta_rows, schema=["temp_seqid", "seqid", "start_win"], orient="row")


def _skani_like_from_blast(
    hits_df: pl.DataFrame, seq_lengths: dict, round_digits=3, overlap_on="query", overlap_tol=0
) -> pl.DataFrame:
    def _norm_iv(a, b):
        a = int(a)
        b = int(b)
        lo, hi = (a - 1, b) if a <= b else (b - 1, a)
        return lo, hi

    def _overlap_len(iv1, iv2):
        a1, a2 = iv1
        b1, b2 = iv2
        return max(0, min(a2, b2) - max(a1, b1))

    def _merge_coverage_len(intervals):
        if not intervals:
            return 0
        ints = sorted(intervals)
        total = 0
        s, e = ints[0]
        for x, y in ints[1:]:
            if x <= e:
                e = max(e, y)
            else:
                total += e - s
                s, e = x, y
        total += e - s
        return total

    def _select_non_overlapping_blast(hits, on="query", tol=0):
        key_iv = ("qstart", "qend") if on == "query" else ("sstart", "send")
        order = sorted(
            range(len(hits)),
            key=lambda i: (
                -(hits[i].get("bitscore") or 0.0),
                -(hits[i].get("length") or 0),
                -(hits[i].get("pident") or 0.0),
            ),
        )
        selected, ivs = [], []
        for i in order:
            a = hits[i].get(key_iv[0])
            b = hits[i].get(key_iv[1])
            if a is None or b is None:
                continue
            iv = _norm_iv(a, b)
            if all(_overlap_len(iv, siv) <= tol for siv in ivs):
                selected.append(i)
                ivs.append(iv)
        return selected

    if hits_df.height == 0:
        return pl.DataFrame(
            {
                "Ref_name": [],
                "Query_name": [],
                "ANI": [],
                "Align_fraction_ref": [],
                "Align_fraction_query": [],
            }
        )

    df = hits_df.clone()
    float_cols = ["pident", "bitscore"]
    int_cols = ["length", "qstart", "qend", "sstart", "send", "qlen", "slen"]
    cast_exprs = []
    for c in float_cols:
        if c in df.columns:
            cast_exprs.append(pl.col(c).cast(pl.Float64))
    for c in int_cols:
        if c in df.columns:
            cast_exprs.append(pl.col(c).cast(pl.Int64))
    if cast_exprs:
        df = df.with_columns(cast_exprs)

    pairs = df.select(["qseqid", "sseqid"]).unique()
    out = []
    for q_id, t_id in pairs.iter_rows():
        g = df.filter((pl.col("qseqid") == q_id) & (pl.col("sseqid") == t_id))
        hits = list(g.iter_rows(named=True))
        sel_idx = _select_non_overlapping_blast(hits, on=overlap_on, tol=overlap_tol)
        sel = [hits[i] for i in sel_idx]

        sum_blen = sum(int(h["length"]) for h in sel if not _is_na(h.get("length")))
        sum_mlen = sum(
            int(round((h["pident"] / 100.0) * h["length"]))
            for h in sel
            if not _is_na(h.get("pident")) and not _is_na(h.get("length"))
        )
        ani = (100.0 * sum_mlen / sum_blen) if sum_blen > 0 else float("nan")

        q_ivs = [
            _norm_iv(h["qstart"], h["qend"])
            for h in sel
            if not _is_na(h.get("qstart")) and not _is_na(h.get("qend"))
        ]
        t_ivs = [
            _norm_iv(h["sstart"], h["send"])
            for h in sel
            if not _is_na(h.get("sstart")) and not _is_na(h.get("send"))
        ]
        covered_q = _merge_coverage_len(q_ivs)
        covered_t = _merge_coverage_len(t_ivs)

        def _first(df_grp: pl.DataFrame, col: str):
            if col not in df_grp.columns or df_grp.height == 0:
                return None
            try:
                return df_grp.select(pl.first(col)).to_series().to_list()[0]
            except Exception:
                return None

        q_len = _first(g, "qlen") or seq_lengths.get(q_id)
        t_len = _first(g, "slen") or seq_lengths.get(t_id)

        frac_q = (covered_q / q_len) * 100 if q_len else float("nan")
        frac_t = (covered_t / t_len) * 100 if t_len else float("nan")

        out.append(
            {
                "Ref_name": t_id,
                "Query_name": q_id,
                "ANI": None if math.isnan(ani) else round(ani, round_digits),
                "Align_fraction_ref": None if math.isnan(frac_t) else round(frac_t, round_digits),
                "Align_fraction_query": None if math.isnan(frac_q) else round(frac_q, round_digits),
            }
        )

    return pl.DataFrame(out)


def _run_mappy_target_block(args):
    fasta_path, preset, min_mapq, mm2_threads_per_worker, tid, q_ids = args
    idx = SeqIO.index(fasta_path, "fasta")
    try:
        t_seq = str(idx[tid].seq)
        al = mp.Aligner(seq=t_seq, preset=preset, n_threads=mm2_threads_per_worker)
        if not al:
            return []

        out = []
        for qid in q_ids:
            q_seq = str(idx[qid].seq)
            for h in al.map(q_seq):
                if getattr(h, "mapq", 0) < min_mapq:
                    continue

                mlen = getattr(h, "mlen", None)
                blen = getattr(h, "blen", None)
                pid_mlen_blen = (100.0 * mlen / blen) if (mlen is not None and blen) else None

                NM = getattr(h, "NM", None)
                q_st = getattr(h, "q_st", None)
                q_en = getattr(h, "q_en", None)
                aln_len_q = (q_en - q_st) if (q_en is not None and q_st is not None) else None
                pid_nm = (
                    (100.0 * (1.0 - NM / aln_len_q))
                    if (NM is not None and aln_len_q and aln_len_q > 0)
                    else None
                )

                out.append(
                    {
                        "query_id": qid,
                        "target_id": tid,
                        "r_st": getattr(h, "r_st", None),
                        "r_en": getattr(h, "r_en", None),
                        "q_st": q_st,
                        "q_en": q_en,
                        "strand": getattr(h, "strand", 1),
                        "mapq": getattr(h, "mapq", None),
                        "cigar": getattr(h, "cigar_str", getattr(h, "cigar", None)),
                        "is_primary": getattr(h, "is_primary", None),
                        "ctg_len": getattr(h, "ctg_len", None),
                        "qry_len": getattr(h, "qry_len", None),
                        "blen": blen,
                        "mlen": mlen,
                        "NM": NM,
                        "pid": pid_mlen_blen if pid_mlen_blen is not None else pid_nm,
                    }
                )
        return out
    finally:
        idx.close()


def _skani_like_from_mappy(
    hits_df: pl.DataFrame,
    fasta_path: str,
    round_digits: int = 3,
    overlap_on: str = "query",
    overlap_tol: int = 0,
    require_primary: bool = True,
    return_percent: bool = True,
) -> pl.DataFrame:
    if hits_df.height == 0:
        return pl.DataFrame(
            columns=["Ref_name", "Query_name", "ANI", "Align_fraction_ref", "Align_fraction_query"]
        )

    df = hits_df.clone()
    need = [
        "query_id",
        "target_id",
        "q_st",
        "q_en",
        "r_st",
        "r_en",
        "mlen",
        "blen",
        "mapq",
        "is_primary",
        "NM",
    ]
    add_exprs = []
    for c in need:
        if c not in df.columns:
            add_exprs.append(pl.lit(None).alias(c))
    if add_exprs:
        df = df.with_columns(add_exprs)
    df = df.with_columns(
        [
            pl.col("q_st").cast(pl.Int64),
            pl.col("q_en").cast(pl.Int64),
            pl.col("r_st").cast(pl.Int64),
            pl.col("r_en").cast(pl.Int64),
            pl.col("mlen").cast(pl.Int64),
            pl.col("blen").cast(pl.Int64),
            pl.col("mapq").cast(pl.Int64),
            pl.col("NM").cast(pl.Int64),
        ]
    )

    if require_primary and "is_primary" in df.columns:
        df = df.filter(pl.col("is_primary").fill_null(True))

    lengths = _fasta_lengths(fasta_path)

    def _norm_iv(a, b):
        a = int(a)
        b = int(b)
        return (a, b) if a <= b else (b, a)

    def _overlap_len(iv1, iv2):
        a1, a2 = iv1
        b1, b2 = iv2
        return max(0, min(a2, b2) - max(a1, b1))

    def _merge_coverage_len(intervals):
        if not intervals:
            return 0
        ints = [_norm_iv(a, b) for a, b in intervals if a is not None and b is not None]
        if not ints:
            return 0
        ints.sort()
        total = 0
        s, e = ints[0]
        for x, y in ints[1:]:
            if x <= e:
                e = max(e, y)
            else:
                total += e - s
                s, e = x, y
        total += e - s
        return total

    def _select_non_overlapping(hits, on="query", tol=0):
        key_iv = ("q_st", "q_en") if on == "query" else ("r_st", "r_en")
        order = sorted(
            range(len(hits)),
            key=lambda i: (
                -(hits[i].get("mapq") or 0),
                -(hits[i].get("mlen") or 0),
                -(hits[i].get("blen") or 0),
            ),
        )
        selected, ivs = [], []
        for i in order:
            a = hits[i].get(key_iv[0])
            b = hits[i].get(key_iv[1])
            if a is None or b is None:
                continue
            iv = _norm_iv(a, b)
            if all(_overlap_len(iv, siv) <= tol for siv in ivs):
                selected.append(i)
                ivs.append(iv)
        return selected

    out = []
    # Use polars partition_by for grouping (returns dict of group_key -> DataFrame)
    for group_df in df.partition_by(["query_id", "target_id"], maintain_order=True):
        q_id = group_df["query_id"][0]
        t_id = group_df["target_id"][0]
        hits = group_df.to_dicts()
        sel_idx = _select_non_overlapping(hits, on=overlap_on, tol=overlap_tol)
        sel = [hits[i] for i in sel_idx]

        sum_blen, sum_mlen = 0, 0
        for h in sel:
            bl = h.get("blen")
            if _is_na(bl) or bl is None:
                qlen = (
                    (h.get("q_en") - h.get("q_st"))
                    if (h.get("q_en") is not None and h.get("q_st") is not None)
                    else None
                )
                rlen = (
                    (h.get("r_en") - h.get("r_st"))
                    if (h.get("r_en") is not None and h.get("r_st") is not None)
                    else None
                )
                bl = qlen if (qlen is not None) else rlen
            ml = h.get("mlen")
            if (_is_na(ml) or ml is None) and bl is not None and h.get("NM") is not None:
                ml = max(0, int(bl) - int(h.get("NM")))
            if bl is not None:
                sum_blen += int(bl)
            if ml is not None:
                sum_mlen += int(ml)

        ani = (100.0 * sum_mlen / sum_blen) if sum_blen > 0 else float("nan")

        q_ivs = [
            _norm_iv(h["q_st"], h["q_en"])
            for h in sel
            if not _is_na(h.get("q_st")) and not _is_na(h.get("q_en"))
        ]
        t_ivs = [
            _norm_iv(h["r_st"], h["r_en"])
            for h in sel
            if not _is_na(h.get("r_st")) and not _is_na(h.get("r_en"))
        ]
        covered_q = _merge_coverage_len(q_ivs)
        covered_t = _merge_coverage_len(t_ivs)

        q_len = lengths.get(q_id)
        t_len = lengths.get(t_id)
        frac_q = (covered_q / q_len) if q_len else float("nan")
        frac_t = (covered_t / t_len) if t_len else float("nan")
        if return_percent:
            frac_q *= 100.0
            frac_t *= 100.0

        out.append(
            {
                "Ref_name": t_id,
                "Query_name": q_id,
                "ANI": None if math.isnan(ani) else round(ani, round_digits),
                "Align_fraction_ref": None if math.isnan(frac_t) else round(frac_t, round_digits),
                "Align_fraction_query": None if math.isnan(frac_q) else round(frac_q, round_digits),
            }
        )

    return pl.DataFrame(out)


def run_pairwise_nt(
    all_neigh: pl.DataFrame | None,
    all_gff: object | None = None,
    output_dir: str = "pairwise_nt_out",
    nt_aln_mode: object | None = None,
    ani_mode: str | None = None,
    nt_links: bool = True,
    threads: int | None = None,
    blast_task: str = "blastn",
    evalue: float = 1e-5,
    perc_identity: int = 0,
    word_size: int | None = None,
    soft_masking: str = "false",
    dust: str = "no",
    overlap_on: str = "query",
    overlap_tol: int = 0,
    write_outputs: bool = True,
    mm2_preset: str = "asm20",
    mm2_min_mapq: int = 0,
    mm2_threads_per_worker: int = 1,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """
    Run pairwise nucleotide comparisons from a single FASTA and return (skani_like_df, alignment_rows_df).
    """

    def _map_and_write_pairwise_ani_uid(df: pl.DataFrame, mode_suffix: str | None = None):
        if df is None or df.height == 0:
            return None
        if all_neigh is None or all_neigh.height == 0:
            return None

        temp_map = {}
        if "temp_seqid" in all_neigh.columns and "unique_id" in all_neigh.columns:
            for r in (
                all_neigh.select([pl.col("temp_seqid").cast(pl.Utf8), "unique_id"])
                .drop_nulls()
                .iter_rows()
            ):
                temp_map[str(r[0])] = r[1]

        seq_map = {}
        if "seqid" in all_neigh.columns and "unique_id" in all_neigh.columns:
            seq_counts = (
                all_neigh.group_by("seqid")
                .agg(pl.col("unique_id").n_unique().alias("n"))
                .filter(pl.col("n") == 1)
            )
            good_seqids = set(seq_counts.select("seqid").to_series().to_list())
            if good_seqids:
                for r in (
                    all_neigh.filter(pl.col("seqid").is_in(list(good_seqids)))
                    .select(["seqid", "unique_id"])
                    .iter_rows()
                ):
                    seq_map[str(r[0])] = r[1]

        def _map_name_to_uid(name):
            if _is_na(name):
                return None
            s = str(name)
            if s in temp_map:
                return temp_map[s]
            if s in seq_map:
                return seq_map[s]
            return None

        dfc = df.clone()
        if "Ref_name" not in dfc.columns or "Query_name" not in dfc.columns:
            return None

        dfc = dfc.with_columns(
            [
                pl.col("Ref_name").map_elements(_map_name_to_uid).alias("ref_uid"),
                pl.col("Query_name").map_elements(_map_name_to_uid).alias("qry_uid"),
            ]
        ).drop_nulls(subset=["ref_uid", "qry_uid"])
        if dfc.height == 0:
            return None

        dfc = dfc.with_columns(
            [
                pl.min_horizontal(["ref_uid", "qry_uid"]).alias("A"),
                pl.max_horizontal(["ref_uid", "qry_uid"]).alias("B"),
            ]
        )

        agg = dfc.group_by(["A", "B"]).agg(
            [
                pl.col("ANI").mean().alias("ANI"),
                pl.col("Align_fraction_ref").mean().alias("Align_fraction_ref"),
                pl.col("Align_fraction_query").mean().alias("Align_fraction_query"),
            ]
        )

        if write_outputs:
            suf = mode_suffix or str(nt_aln_mode)
            out_path = out / f"pairwise_ani_uid_{suf}.tsv"
            with contextlib.suppress(Exception):
                agg.write_csv(out_path, separator="\t", include_header=False)

        return agg

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    tmp_dir = out / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    threads = threads or max(1, os.cpu_count() or 1)

    nb_dir = out / "neighborhood"
    nb_dir.mkdir(parents=True, exist_ok=True)
    fasta_path = nb_dir / "neighborhoods.fasta"
    if not fasta_path.exists():
        if all_neigh is None or all_neigh.height == 0:
            raise FileNotFoundError(
                f"neighborhood FASTA not found at {fasta_path} and no `all_neigh` provided"
            )

        if "temp_seqid" in all_neigh.columns and "sequence" in all_neigh.columns:
            df_for_fasta = (
                all_neigh[["temp_seqid", "sequence"]].drop_nulls().unique(subset=["temp_seqid"])
            )
            df_for_fasta.to_fasta("temp_seqid", "sequence", str(fasta_path))
        elif "seqid" in all_neigh.columns and "sequence" in all_neigh.columns:
            df_for_fasta = all_neigh[["seqid", "sequence"]].drop_nulls().unique(subset=["seqid"])
            df_for_fasta.to_fasta("seqid", "sequence", str(fasta_path))

    ids = _fasta_ids(str(fasta_path))
    id_pos = {cid: i for i, cid in enumerate(ids)}
    seq_lengths = _fasta_lengths(str(fasta_path))

    if ani_mode is not None and str(ani_mode).lower() == "skani":
        info(f"Running skani triangle on {fasta_path}")
        cmd = [
            "skani",
            "triangle",
            "-i",
            str(fasta_path),
            "-E",
            "--small-genomes",
            "-t",
            str(threads),
        ]
        info(f"Running: {' '.join(cmd)}")
        skani_log = out / "skani_triangle.log"
        with open(skani_log, "w") as logfh:
            subprocess.run(cmd, check=True, stdout=logfh, stderr=subprocess.DEVNULL, text=True)
        try:
            df = pl.read_csv(skani_log, separator="\t", has_header=True)
        except Exception:
            empty_align = pl.DataFrame(
                columns=["query", "query_start", "query_end", "ref", "ref_start", "ref_end", "ani"]
            )
            return (
                pl.DataFrame(
                    columns=[
                        "Ref_name",
                        "Query_name",
                        "ANI",
                        "Align_fraction_ref",
                        "Align_fraction_query",
                    ]
                ),
                empty_align,
            )

        if "Align_fraction_query" not in df.columns and "Align_fraction_ref" in df.columns:
            df = df.with_columns(pl.col("Align_fraction_ref").alias("Align_fraction_query"))
        if "Ref_name" in df.columns:
            df = df.with_columns(
                pl.col("Ref_name").cast(pl.Utf8).str.split("/").list.last().alias("Ref_name")
            )
        if "Query_name" in df.columns:
            df = df.with_columns(
                pl.col("Query_name").cast(pl.Utf8).str.split("/").list.last().alias("Query_name")
            )
        cast_exprs = []
        if "ANI" in df.columns:
            cast_exprs.append(pl.col("ANI").cast(pl.Float64))
        if "Align_fraction_ref" in df.columns:
            cast_exprs.append(pl.col("Align_fraction_ref").cast(pl.Float64))
        if "Align_fraction_query" in df.columns:
            cast_exprs.append(pl.col("Align_fraction_query").cast(pl.Float64))
        if cast_exprs:
            df = df.with_columns(cast_exprs)
        if "Ref_name" in df.columns and "Query_name" in df.columns:
            df = df.filter(pl.col("Ref_name") != pl.col("Query_name"))
        required = ["Ref_name", "Query_name", "ANI", "Align_fraction_ref", "Align_fraction_query"]
        pairwise_ani = (
            df.select(required)
            if all(c in df.columns for c in required)
            else pl.DataFrame({c: [] for c in required})
        )
        try:
            agg_uid = _map_and_write_pairwise_ani_uid(pairwise_ani, mode_suffix="skani")
            if agg_uid is not None:
                pairwise_ani = agg_uid
        except Exception:
            pass
        empty_align = pl.DataFrame(
            columns=["query", "query_start", "query_end", "ref", "ref_start", "ref_end", "ani"]
        )
        return pairwise_ani, empty_align

    if (
        nt_aln_mode
        and nt_aln_mode.lower() in ("blastn", "blast")
        or (ani_mode and str(ani_mode).lower() == "blastn")
    ):
        db_prefix = tmp_dir / "db"
        blast_tsv = tmp_dir / "allvsall.outfmt6.tsv"
        info("Building BLAST database for nucleotide comparisons...")
        subprocess.run(
            [
                "makeblastdb",
                "-in",
                fasta_path,
                "-dbtype",
                "nucl",
                "-parse_seqids",
                "-out",
                str(db_prefix),
            ],
            check=True,
            capture_output=True,
        )

        outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
        cmd = [
            "blastn",
            "-query",
            str(fasta_path),
            "-db",
            str(db_prefix),
            "-task",
            blast_task,
            "-evalue",
            str(evalue),
            "-soft_masking",
            soft_masking,
            "-dust",
            dust,
            "-outfmt",
            outfmt,
            "-num_threads",
            str(threads),
        ]
        if perc_identity and perc_identity > 0:
            cmd += ["-perc_identity", str(perc_identity)]
        if word_size:
            cmd += ["-word_size", str(word_size)]
        info("Running BLAST nucleotide all-vs-all search...")
        with open(blast_tsv, "w") as fh:
            subprocess.run(cmd, check=True, stdout=fh, stderr=subprocess.PIPE)

        blast_cols = [
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
            "qlen",
            "slen",
        ]
        raw = pl.read_csv(
            blast_tsv,
            separator="\t",
            has_header=False,
            new_columns=blast_cols,
            dtypes={"qseqid": pl.Utf8, "sseqid": pl.Utf8},
        )
        raw = raw.filter(pl.col("qseqid") != pl.col("sseqid"))
        raw = raw.with_columns(
            [
                pl.col("qseqid")
                .map_elements(lambda x: id_pos.get(str(x)), return_dtype=pl.Int64)
                .alias("q_idx"),
                pl.col("sseqid")
                .map_elements(lambda x: id_pos.get(str(x)), return_dtype=pl.Int64)
                .alias("s_idx"),
            ]
        ).drop_nulls(subset=["q_idx", "s_idx"])
        raw = raw.filter(pl.col("q_idx") < pl.col("s_idx"))

        hits_blast_df = raw.select(blast_cols)
        if (
            all_neigh is not None
            and all_neigh.height > 0
            and "temp_seqid" in all_neigh.columns
            and "seqid" in all_neigh.columns
        ):
            hits_blast_df = hits_blast_df.join(
                all_neigh[["temp_seqid", "seqid"]],
                left_on="qseqid",
                right_on="temp_seqid",
                how="left",
            )

        info(f"Raw BLAST HSPs (unique unordered, no self): {hits_blast_df.height}")
        if write_outputs:
            hits_blast_df.write_csv(
                tmp_dir / "pairwise_hits_blast.tsv", separator="\t", include_header=False
            )

        skani_df = _skani_like_from_blast(
            hits_blast_df,
            seq_lengths,
            round_digits=3,
            overlap_on=overlap_on,
            overlap_tol=overlap_tol,
        )
        if not nt_links:
            empty_align = pl.DataFrame(
                {
                    "query": [],
                    "query_start": [],
                    "query_end": [],
                    "ref": [],
                    "ref_start": [],
                    "ref_end": [],
                    "ani": [],
                }
            )
            try:
                agg_uid = _map_and_write_pairwise_ani_uid(skani_df, mode_suffix="blastn")
                if agg_uid is not None:
                    skani_df = agg_uid
            except Exception:
                pass
            return skani_df, empty_align

        align_rows = hits_blast_df.rename(
            {
                "qseqid": "query",
                "qstart": "query_start",
                "qend": "query_end",
                "sseqid": "ref",
                "sstart": "ref_start",
                "send": "ref_end",
                "pident": "ani",
            }
        ).select(["query", "query_start", "query_end", "ref", "ref_start", "ref_end", "ani"])

        if all_neigh is not None and all_neigh.height > 0:
            start_map, id_map = {}, {}
            for nr in all_neigh.iter_rows(named=True):
                temp = (
                    str(nr.get("temp_seqid"))
                    if ("temp_seqid" in all_neigh.columns and nr.get("temp_seqid") is not None)
                    else None
                )
                seqid = (
                    str(nr.get("seqid"))
                    if ("seqid" in all_neigh.columns and nr.get("seqid") is not None)
                    else temp
                )
                start_win = (
                    int(nr.get("start_win"))
                    if ("start_win" in all_neigh.columns and nr.get("start_win") is not None)
                    else 0
                )
                for key in (temp, seqid):
                    if not key:
                        continue
                    start_map[key] = start_win
                    start_map[Path(key).name] = start_win
                    start_map[Path(key).stem] = start_win
                    id_map[key] = seqid
                    id_map[Path(key).name] = seqid
                    id_map[Path(key).stem] = seqid

            def _map_to_seqid(name: str) -> str:
                if _is_na(name):
                    return name
                name = str(name)
                return (
                    id_map.get(name)
                    or id_map.get(Path(name).name)
                    or id_map.get(Path(name).stem)
                    or name
                )

            def _find_offset(name: str) -> int:
                if _is_na(name):
                    return 0
                canonical = _map_to_seqid(name)
                return (
                    start_map.get(canonical)
                    or start_map.get(Path(canonical).name)
                    or start_map.get(Path(canonical).stem)
                    or 0
                )

            align_rows = align_rows.with_columns(
                [
                    pl.col("query").map_elements(_map_to_seqid).alias("query"),
                    pl.col("ref").map_elements(_map_to_seqid).alias("ref"),
                ]
            )

            align_rows = align_rows.with_columns(
                [
                    pl.col("query").map_elements(_find_offset).alias("q_offset"),
                    pl.col("ref").map_elements(_find_offset).alias("r_offset"),
                ]
            )
            align_rows = align_rows.with_columns(
                [
                    (pl.col("query_start").cast(pl.Int64) + pl.col("q_offset")).alias(
                        "query_start"
                    ),
                    (pl.col("query_end").cast(pl.Int64) + pl.col("q_offset")).alias("query_end"),
                    (pl.col("ref_start").cast(pl.Int64) + pl.col("r_offset")).alias("ref_start"),
                    (pl.col("ref_end").cast(pl.Int64) + pl.col("r_offset")).alias("ref_end"),
                ]
            ).drop(["q_offset", "r_offset"])

            if skani_df.height > 0:
                skani_df = skani_df.with_columns(
                    [
                        pl.col("Ref_name").map_elements(_map_to_seqid).alias("Ref_name"),
                        pl.col("Query_name").map_elements(_map_to_seqid).alias("Query_name"),
                    ]
                )
            try:
                agg_uid = _map_and_write_pairwise_ani_uid(skani_df, mode_suffix="blastn")
                if agg_uid is not None:
                    skani_df = agg_uid
            except Exception:
                pass

        return skani_df, align_rows

    if nt_aln_mode and nt_aln_mode.lower() in ("intergenic_blast", "intergenic_blastn"):
        inter_fa = out / "intergenic.fasta"
        info(f"Writing intergenic FASTA to {inter_fa}")
        meta_df = _write_intergenic_fasta(
            all_gff, str(fasta_path), str(inter_fa), all_neigh=all_neigh
        )
        if meta_df.height == 0:
            info("No intergenic sequences produced; exiting")
            return (
                pl.DataFrame(
                    columns=[
                        "Ref_name",
                        "Query_name",
                        "ANI",
                        "Align_fraction_ref",
                        "Align_fraction_query",
                    ]
                ),
                pl.DataFrame(),
            )

        db_prefix = tmp_dir / "intergenic_db"
        info(f"Making BLAST DB for intergenic sequences at {db_prefix}")
        subprocess.run(
            ["makeblastdb", "-in", str(inter_fa), "-dbtype", "nucl", "-out", str(db_prefix)],
            check=True,
            capture_output=True,
        )

        blast_tsv = tmp_dir / "intergenic_allvsall.outfmt6.tsv"
        outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
        cmd = [
            "blastn",
            "-query",
            str(inter_fa),
            "-db",
            str(db_prefix),
            "-task",
            "blastn-short",
            "-evalue",
            "1e-3",
            "-word_size",
            "4",
            "-reward",
            "1",
            "-penalty",
            "-1",
            "-gapopen",
            "2",
            "-gapextend",
            "2",
            "-dust",
            "no",
            "-soft_masking",
            "false",
            "-outfmt",
            outfmt,
            "-num_threads",
            str(threads),
        ]
        info("Running intergenic BLAST")
        with open(blast_tsv, "w") as fh:
            subprocess.run(cmd, check=True, stdout=fh, stderr=subprocess.PIPE)

        blast_cols = [
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
            "qlen",
            "slen",
        ]
        raw = pl.read_csv(
            blast_tsv,
            separator="\t",
            has_header=False,
            new_columns=blast_cols,
            dtypes={"qseqid": pl.Utf8, "sseqid": pl.Utf8},
        )

        meta = meta_df.clone().rename({"temp_seqid": "temp_id", "seqid": "canonical_seqid"})
        # Join for query side
        raw = raw.join(
            meta.select(["temp_id", "canonical_seqid", "start_win"]),
            left_on="qseqid",
            right_on="temp_id",
            how="left",
        )
        raw = raw.rename({"canonical_seqid": "q_canonical", "start_win": "q_start_win"})
        # Drop temp_id only if it exists (Polars may drop it automatically in join)
        if "temp_id" in raw.columns:
            raw = raw.drop(["temp_id"])

        # Join for subject side
        raw = raw.join(
            meta.select(["temp_id", "canonical_seqid", "start_win"]),
            left_on="sseqid",
            right_on="temp_id",
            how="left",
        )
        raw = raw.rename({"canonical_seqid": "s_canonical", "start_win": "s_start_win"})
        # Drop temp_id only if it exists
        if "temp_id" in raw.columns:
            raw = raw.drop(["temp_id"])

        raw = raw.filter(pl.col("q_canonical").is_not_null() & pl.col("s_canonical").is_not_null())
        raw = raw.filter(pl.col("q_canonical") != pl.col("s_canonical"))

        raw = raw.with_columns(
            [
                (pl.col("qstart").cast(pl.Int64) + pl.col("q_start_win").cast(pl.Int64) - 1).alias(
                    "qstart_g"
                ),
                (pl.col("qend").cast(pl.Int64) + pl.col("q_start_win").cast(pl.Int64) - 1).alias(
                    "qend_g"
                ),
                (pl.col("sstart").cast(pl.Int64) + pl.col("s_start_win").cast(pl.Int64) - 1).alias(
                    "sstart_g"
                ),
                (pl.col("send").cast(pl.Int64) + pl.col("s_start_win").cast(pl.Int64) - 1).alias(
                    "send_g"
                ),
            ]
        )

        cond = pl.col("q_canonical") > pl.col("s_canonical")
        hits_blast_df = raw.with_columns(
            [
                pl.when(cond)
                .then(pl.col("s_canonical"))
                .otherwise(pl.col("q_canonical"))
                .alias("qseqid"),
                pl.when(cond)
                .then(pl.col("q_canonical"))
                .otherwise(pl.col("s_canonical"))
                .alias("sseqid"),
                pl.when(cond)
                .then(pl.col("sstart_g"))
                .otherwise(pl.col("qstart_g"))
                .alias("qstart"),
                pl.when(cond).then(pl.col("send_g")).otherwise(pl.col("qend_g")).alias("qend"),
                pl.when(cond)
                .then(pl.col("qstart_g"))
                .otherwise(pl.col("sstart_g"))
                .alias("sstart"),
                pl.when(cond).then(pl.col("qend_g")).otherwise(pl.col("send_g")).alias("send"),
            ]
        ).select(
            [
                "qseqid",
                "sseqid",
                "pident",
                "length",
                "qstart",
                "qend",
                "sstart",
                "send",
                "qlen",
                "slen",
                "mismatch",
                "gapopen",
                "evalue",
                "bitscore",
            ]
        )

        info(f"Intergenic BLAST HSPs (mapped to seqid): {hits_blast_df.height}")
        if write_outputs:
            hits_blast_df.write_csv(
                out / "pairwise_hits_intergenic_blast.tsv", separator="\t", include_header=False
            )

        seq_lengths_map = {}
        if (
            all_neigh is not None
            and all_neigh.height > 0
            and "seqid" in all_neigh.columns
            and "sequence" in all_neigh.columns
        ):
            for nr in all_neigh.iter_rows(named=True):
                sid = str(nr.get("seqid")) if nr.get("seqid") is not None else None
                if not sid:
                    continue
                seq = nr.get("sequence") if "sequence" in nr else None
                if seq is not None:
                    seq_lengths_map[sid] = len(seq)

        skani_df = _skani_like_from_blast(
            hits_blast_df,
            seq_lengths_map or seq_lengths,
            round_digits=3,
            overlap_on=overlap_on,
            overlap_tol=overlap_tol,
        )
        if write_outputs:
            skani_df.write_csv(
                out / "skani_like_intergenic_blast.tsv", separator="\t", include_header=False
            )

        align_rows = hits_blast_df.rename(
            {
                "qseqid": "query",
                "qstart": "query_start",
                "qend": "query_end",
                "sseqid": "ref",
                "sstart": "ref_start",
                "send": "ref_end",
                "pident": "ani",
            }
        ).select(["query", "query_start", "query_end", "ref", "ref_start", "ref_end", "ani"])

        return skani_df, align_rows

    if nt_aln_mode and nt_aln_mode.lower() == "fastani":
        split_dir = out / "ani_split"
        split_dir.mkdir(parents=True, exist_ok=True)

        genome_files = []
        file_to_seqid = {}
        temp_to_seqid = {}
        start_map_init = {}
        if all_neigh is not None and all_neigh.height > 0:
            for nr in all_neigh.iter_rows(named=True):
                temp = (
                    str(nr.get("temp_seqid"))
                    if ("temp_seqid" in all_neigh.columns and nr.get("temp_seqid") is not None)
                    else None
                )
                seqid = (
                    str(nr.get("seqid"))
                    if ("seqid" in all_neigh.columns and nr.get("seqid") is not None)
                    else temp
                )
                if temp:
                    temp_to_seqid[temp] = seqid
                if seqid and "start_win" in nr and nr.get("start_win") is not None:
                    start_map_init[seqid] = int(nr.get("start_win"))

        for idx, record in enumerate(SeqIO.parse(str(fasta_path), "fasta")):
            header = record.id
            filename = (
                "".join(c if c.isalnum() or c in "-._" else "_" for c in header)
                + f"_idx{idx}.fasta"
            )
            out_path = split_dir / filename
            SeqIO.write(record, str(out_path), "fasta")
            genome_files.append(str(out_path))
            canonical = temp_to_seqid.get(header, header)
            file_to_seqid[out_path.name] = canonical
            file_to_seqid[out_path.stem] = canonical
            san_hdr = "".join(c if c.isalnum() or c in "-._" else "_" for c in header)
            file_to_seqid[san_hdr] = canonical

        file_list_path = out / "fastani_genome_list.txt"
        with open(file_list_path, "w") as fh:
            for p in genome_files:
                fh.write(p + "\n")

        fastani_output = out / "fastani_output.tsv"
        cmd = [
            "fastANI",
            "--ql",
            str(file_list_path),
            "--rl",
            str(file_list_path),
            "-o",
            str(fastani_output),
            "-t",
            str(threads),
        ]
        fastani_log = out / "fastani_all.log"
        with open(fastani_log, "w") as logfh:
            subprocess.run(cmd, check=True, stdout=logfh, stderr=logfh, text=True)

        if not fastani_output.exists():
            raise FileNotFoundError(f"fastANI did not produce expected output at {fastani_output}")

        df = pl.read_csv(
            fastani_output,
            separator="\t",
            has_header=False,
            new_columns=["query", "reference", "ani", "frags_matched", "frags_total_query"],
        ).with_columns(
            [
                pl.col("ani").cast(pl.Float64),
                pl.col("frags_matched").cast(pl.Int64),
                pl.col("frags_total_query").cast(pl.Int64),
            ]
        )

        df_rev = df.rename(
            {"query": "reference", "reference": "query", "frags_total_query": "frags_total_ref"}
        )
        dfj = df.join(
            df_rev.select(["query", "reference", "frags_total_ref"]),
            on=["query", "reference"],
            how="left",
        )

        skani_like_df = (
            dfj.with_columns(
                [
                    (pl.col("frags_matched") / pl.col("frags_total_query"))
                    .fill_null(0.0)
                    .alias("Align_fraction_query"),
                    (pl.col("frags_matched") / pl.col("frags_total_ref"))
                    .fill_null(pl.col("frags_matched") / pl.col("frags_total_query"))
                    .alias("Align_fraction_ref"),
                ]
            )
            .with_columns(
                [
                    (pl.col("Align_fraction_query") * 100.0).alias("Align_fraction_query"),
                    (pl.col("Align_fraction_ref") * 100.0).alias("Align_fraction_ref"),
                    pl.col("reference").str.split("/").list.last().alias("Ref_name"),
                    pl.col("query").str.split("/").list.last().alias("Query_name"),
                ]
            )
            .select(["Ref_name", "Query_name", "ani", "Align_fraction_ref", "Align_fraction_query"])
        )

        skani_like_df = skani_like_df.filter(pl.col("Ref_name") != pl.col("Query_name"))

        agg = (
            skani_like_df.with_columns(
                [
                    pl.min_horizontal(["Ref_name", "Query_name"]).alias("A"),
                    pl.max_horizontal(["Ref_name", "Query_name"]).alias("B"),
                ]
            )
            .group_by(["A", "B"])
            .agg(
                [
                    pl.col("ani").mean().alias("ANI"),
                    pl.col("Align_fraction_ref").mean().alias("Align_fraction_ref"),
                    pl.col("Align_fraction_query").mean().alias("Align_fraction_query"),
                ]
            )
            .with_columns([pl.col("A").alias("Ref_name"), pl.col("B").alias("Query_name")])
            .select(["Ref_name", "Query_name", "ANI", "Align_fraction_ref", "Align_fraction_query"])
        )
        pairwise_ani = agg

        if write_outputs:
            pairwise_ani.write_csv(
                out / "pairwise_ani_fastani.tsv", separator="\t", include_header=False
            )

        try:
            agg_uid = _map_and_write_pairwise_ani_uid(pairwise_ani, mode_suffix="fastani")
            if agg_uid is not None:
                pairwise_ani = agg_uid
        except Exception:
            pass

        if not nt_links:
            empty_vis = pl.DataFrame(
                columns=["query", "ref", "ani", "query_start", "query_end", "ref_start", "ref_end"]
            )
            return pairwise_ani, empty_vis

        work_dir = out / "fastani_pairwise_visual"
        work_dir.mkdir(parents=True, exist_ok=True)

        def _resolve_path(name: str) -> Path:
            p = Path(name)
            if p.exists():
                return p
            candidates = [
                out / "ani_split" / name,
                out / "ani_split" / (name + ".fasta"),
                out / "ani_split" / (name + ".fa"),
                out / "neighborhood" / name,
                out / "neighborhood" / (name + ".fasta"),
            ]
            for c in candidates:
                if c.exists():
                    return c
            raise FileNotFoundError(f"Could not resolve path for genome identifier '{name}'")

        dfp = pairwise_ani.clone()
        if "Align_fraction_ref" in dfp.columns or "Align_fraction_query" in dfp.columns:
            dfp = dfp.drop_nulls(
                subset=[
                    c for c in ["Align_fraction_ref", "Align_fraction_query"] if c in dfp.columns
                ]
            )

        visual_rows = []
        seen_pairs = set()

        progress = Progress(
            TextColumn(
                f"[grey53][[/grey53][light_slate_grey]{datetime.now():%H:%M:%S}[/light_slate_grey][grey53]][/grey53]"
            ),
            SpinnerColumn(),
            TextColumn("{task.description}"),
            BarColumn(bar_width=40),
            TimeElapsedColumn(),
            TimeRemainingColumn(),
        )
        task = progress.add_task("[cyan]Running pairwise FastANI visualize...", total=len(dfp))
        progress.start()

        for row in dfp.iter_rows():
            q_name = row[dfp.columns.index("Query_name")]
            r_name = row[dfp.columns.index("Ref_name")]
            pair_key = tuple(sorted([str(q_name), str(r_name)]))
            if pair_key in seen_pairs:
                progress.advance(task)
                continue
            seen_pairs.add(pair_key)

            raw0, raw1 = str(pair_key[0]), str(pair_key[1])
            out_file = work_dir / f"{raw0}__vs__{raw1}.visual"
            if not out_file.exists():
                q_path = _resolve_path(pair_key[0])
                r_path = _resolve_path(pair_key[1])
                temp_base = work_dir / f"{raw0}__vs__{raw1}.fastani"
                cmd = [
                    "fastANI",
                    "-q",
                    str(q_path),
                    "-r",
                    str(r_path),
                    "--visualize",
                    "-o",
                    str(temp_base),
                    "-t",
                    str(threads),
                ]
                temp_log = Path(str(temp_base) + ".fastani.log")
                try:
                    with open(temp_log, "w") as logfh:
                        subprocess.run(cmd, check=True, stdout=logfh, stderr=logfh, text=True)
                except subprocess.CalledProcessError:
                    continue
                visual_file = Path(str(temp_base) + ".visual")
                if visual_file.exists():
                    visual_file.rename(out_file)

            if out_file.exists():
                parsed = pl.read_csv(
                    out_file,
                    separator="\t",
                    has_header=False,
                    new_columns=[
                        "query",
                        "ref",
                        "ani",
                        "na1",
                        "na2",
                        "na3",
                        "query_start",
                        "query_end",
                        "ref_start",
                        "ref_end",
                        "na4",
                        "na5",
                    ],
                    dtypes={"query": pl.Utf8, "ref": pl.Utf8},
                )
                parsed = parsed.select(
                    ["query", "ref", "ani", "query_start", "query_end", "ref_start", "ref_end"]
                ).with_columns(
                    [
                        pl.col("query")
                        .map_elements(lambda s: Path(str(s)).name if s is not None else s)
                        .alias("query"),
                        pl.col("ref")
                        .map_elements(lambda s: Path(str(s)).name if s is not None else s)
                        .alias("ref"),
                    ]
                )
                visual_rows.append(parsed)
            progress.advance(task)

        progress.stop()
        visual_df = (
            pl.concat(visual_rows, how="vertical")
            if visual_rows
            else pl.DataFrame(
                columns=["query", "ref", "ani", "query_start", "query_end", "ref_start", "ref_end"]
            )
        )

        if "ani" in visual_df.columns:
            visual_df = visual_df.with_columns(pl.col("ani").cast(pl.Float64))

        if all_neigh is not None and all_neigh.height > 0:
            id_map = {}
            start_map = {}
            for k, v in file_to_seqid.items():
                id_map[k] = v
                id_map[Path(k).name] = v
                id_map[Path(k).stem] = v
            for seqid, st in start_map_init.items():
                start_map[seqid] = st
                start_map[Path(seqid).name] = st
                start_map[Path(seqid).stem] = st
            for nr in all_neigh.iter_rows(named=True):
                temp = (
                    str(nr.get("temp_seqid"))
                    if ("temp_seqid" in all_neigh.columns and nr.get("temp_seqid") is not None)
                    else None
                )
                seqid = (
                    str(nr.get("seqid"))
                    if ("seqid" in all_neigh.columns and nr.get("seqid") is not None)
                    else temp
                )
                start_win = (
                    int(nr.get("start_win"))
                    if ("start_win" in all_neigh.columns and nr.get("start_win") is not None)
                    else 0
                )
                for key in (temp, seqid):
                    if not key:
                        continue
                    id_map[key] = seqid
                    id_map[Path(key).name] = seqid
                    id_map[Path(key).stem] = seqid
                    start_map[key] = start_win
                    start_map[Path(key).name] = start_win
                    start_map[Path(key).stem] = start_win

            def _map_to_seqid_fast(name: str) -> str:
                if _is_na(name):
                    return name
                name = str(name)
                return (
                    id_map.get(name)
                    or id_map.get(Path(name).name)
                    or id_map.get(Path(name).stem)
                    or name
                )

            def _find_offset_fast(name: str) -> int:
                if _is_na(name):
                    return 0
                name = str(name)
                canonical = _map_to_seqid_fast(name)
                return (
                    start_map.get(canonical)
                    or start_map.get(Path(canonical).name)
                    or start_map.get(Path(canonical).stem)
                    or 0
                )

            pairwise_ani = pairwise_ani.with_columns(
                [
                    pl.col("Ref_name").map_elements(_map_to_seqid_fast).alias("Ref_name"),
                    pl.col("Query_name").map_elements(_map_to_seqid_fast).alias("Query_name"),
                ]
            )

            if visual_df.height > 0:
                visual_df = visual_df.with_columns(
                    [
                        pl.col("query").map_elements(_find_offset_fast).alias("q_off"),
                        pl.col("ref").map_elements(_find_offset_fast).alias("r_off"),
                    ]
                )
                visual_df = visual_df.with_columns(
                    [
                        (pl.col("query_start").cast(pl.Int64) + pl.col("q_off")).alias(
                            "query_start"
                        ),
                        (pl.col("query_end").cast(pl.Int64) + pl.col("q_off")).alias("query_end"),
                        (pl.col("ref_start").cast(pl.Int64) + pl.col("r_off")).alias("ref_start"),
                        (pl.col("ref_end").cast(pl.Int64) + pl.col("r_off")).alias("ref_end"),
                        pl.col("query").map_elements(_map_to_seqid_fast).alias("query"),
                        pl.col("ref").map_elements(_map_to_seqid_fast).alias("ref"),
                    ]
                ).drop(["q_off", "r_off"])

        return pairwise_ani, visual_df

    if nt_aln_mode and nt_aln_mode.lower() in ("minimap2", "mappy"):
        plans = [
            (
                str(fasta_path),
                mm2_preset,
                mm2_min_mapq,
                mm2_threads_per_worker,
                tid,
                ids[: id_pos[tid]],
            )
            for tid in ids
        ]
        rows = []

        with Progress(
            TextColumn(
                f"[grey53][[/grey53][light_slate_grey]{datetime.now():%H:%M:%S}[/light_slate_grey][grey53]][/grey53]"
            ),
            SpinnerColumn(),
            TextColumn("{task.description}"),
            BarColumn(bar_width=40),
            TimeElapsedColumn(),
            TimeRemainingColumn(),
        ) as progress:
            task = progress.add_task("Running minimap2/mappy...", total=len(plans))
            with ProcessPoolExecutor(max_workers=max(1, min(threads, len(plans)))) as ex:
                futs = [ex.submit(_run_mappy_target_block, p) for p in plans]
                for fut in as_completed(futs):
                    rows.extend(fut.result())
                    progress.advance(task)

        hits_df = pl.DataFrame(rows)
        if write_outputs:
            hits_df.write_csv(
                tmp_dir / "pairwise_hits_mappy.tsv", separator="\t", include_header=False
            )

        skani_df = _skani_like_from_mappy(
            hits_df,
            str(fasta_path),
            round_digits=3,
            overlap_on=overlap_on,
            overlap_tol=overlap_tol,
            require_primary=True,
            return_percent=True,
        )
        if write_outputs:
            skani_df.write_csv(
                tmp_dir / "skani_like_mappy.tsv", separator="\t", include_header=False
            )

        try:
            agg_uid = _map_and_write_pairwise_ani_uid(skani_df, mode_suffix="mappy")
            if agg_uid is not None:
                skani_df = agg_uid
        except Exception:
            pass

        align_rows = hits_df.rename(
            {
                "query_id": "query",
                "q_st": "query_start",
                "q_en": "query_end",
                "target_id": "ref",
                "r_st": "ref_start",
                "r_en": "ref_end",
                "pid": "ani",
            }
        ).select(["query", "query_start", "query_end", "ref", "ref_start", "ref_end", "ani"])

        if all_neigh is not None and all_neigh.height > 0:
            start_map, id_map = {}, {}
            for nr in all_neigh.iter_rows(named=True):
                temp = (
                    str(nr.get("temp_seqid"))
                    if ("temp_seqid" in all_neigh.columns and nr.get("temp_seqid") is not None)
                    else None
                )
                seqid = (
                    str(nr.get("seqid"))
                    if ("seqid" in all_neigh.columns and nr.get("seqid") is not None)
                    else temp
                )
                start_win = (
                    int(nr.get("start_win"))
                    if ("start_win" in all_neigh.columns and nr.get("start_win") is not None)
                    else 0
                )
                for key in (temp, seqid):
                    if not key:
                        continue
                    start_map[key] = start_win
                    start_map[Path(key).name] = start_win
                    start_map[Path(key).stem] = start_win
                    id_map[key] = seqid
                    id_map[Path(key).name] = seqid
                    id_map[Path(key).stem] = seqid

            def _find_offset(name: str) -> int:
                if _is_na(name):
                    return 0
                name = str(name)
                canonical = (
                    id_map.get(name)
                    or id_map.get(Path(name).name)
                    or id_map.get(Path(name).stem)
                    or name
                )
                return (
                    start_map.get(canonical)
                    or start_map.get(Path(canonical).name)
                    or start_map.get(Path(canonical).stem)
                    or 0
                )

            def _map_to_seqid(name: str) -> str:
                if _is_na(name):
                    return name
                name = str(name)
                return (
                    id_map.get(name)
                    or id_map.get(Path(name).name)
                    or id_map.get(Path(name).stem)
                    or name
                )

            align_rows = (
                align_rows.with_columns(
                    [
                        pl.col("query").map_elements(_find_offset).alias("q_offset"),
                        pl.col("ref").map_elements(_find_offset).alias("r_offset"),
                    ]
                )
                .with_columns(
                    [
                        (pl.col("query_start").cast(pl.Int64) + pl.col("q_offset")).alias(
                            "query_start"
                        ),
                        (pl.col("query_end").cast(pl.Int64) + pl.col("q_offset")).alias(
                            "query_end"
                        ),
                        (pl.col("ref_start").cast(pl.Int64) + pl.col("r_offset")).alias(
                            "ref_start"
                        ),
                        (pl.col("ref_end").cast(pl.Int64) + pl.col("r_offset")).alias("ref_end"),
                        pl.col("query").map_elements(_map_to_seqid).alias("query"),
                        pl.col("ref").map_elements(_map_to_seqid).alias("ref"),
                    ]
                )
                .drop(["q_offset", "r_offset"])
            )

        return skani_df, align_rows

    raise ValueError("nt_aln_mode must be 'blastn', 'intergenic_blast', 'fastani', or 'mappy'")
