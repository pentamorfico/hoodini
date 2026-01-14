import contextlib
import re
import subprocess
from pathlib import Path

import polars as pl

from hoodini.utils.logging_utils import console


def _safe_name(s: str) -> str:
    s = str(s)
    s = s.strip()
    if "/" in s or "\\" in s:
        s = Path(s).name
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s)


def _resolve_path(name: str, output_dir: Path) -> Path:
    p = Path(name)
    if p.exists():
        return p
    candidates = [
        output_dir / "ani_split" / name,
        output_dir / "ani_split" / (name + ".fasta"),
        output_dir / "ani_split" / (name + ".fa"),
        output_dir / "neighborhood" / name,
        output_dir / "neighborhood" / (name + ".fasta"),
    ]
    for c in candidates:
        if c.exists():
            return c
    raise FileNotFoundError(f"Could not resolve path for genome identifier '{name}'")


def run_nt_links(
    pairwise_ani: pl.DataFrame,
    output_dir: str,
    all_neigh: pl.DataFrame | None = None,
    threads: int = 8,
    evalue: float | None = None,
    keep_temp: bool = False,
) -> pl.DataFrame:
    """Generate pairwise nucleotide visual blocks using fastANI --visualize.

    This function expects a skani-like DataFrame with columns 'Ref_name' and 'Query_name'.
    It will canonicalize reciprocal pairs, run fastANI in a single consistent direction per
    pair, cache the produced .visual files under <output_dir>/fastani_pairwise_visual, parse
    them and return a concatenated DataFrame of alignment blocks.

    Parameters
    - pairwise_ani: DataFrame with at least 'Ref_name' and 'Query_name' columns (skani-like)
    - output_dir: project output folder where run_ani may have produced 'ani_split'
    - all_neigh: optional DataFrame with sequences (not used directly but kept for API parity)
    - threads: threads to pass to fastANI (when supported)
    - keep_temp: if True, don't remove temporary files produced by fastANI

    Returns
    - pandas.DataFrame with columns ['query','ref','ani','query_start','query_end','ref_start','ref_end'] aggregated from all visual files
    """
    out = Path(output_dir)
    work_dir = out / "fastani_pairwise_visual"
    work_dir.mkdir(parents=True, exist_ok=True)

    df = pairwise_ani.copy()
    col_map = {c.lower(): c for c in df.columns}
    q_col = col_map.get("query_name") or col_map.get("query") or "Query_name"
    r_col = col_map.get("ref_name") or col_map.get("reference") or "Ref_name"

    if "Align_fraction_ref" in df.columns or "Align_fraction_query" in df.columns:
        df = df.dropna(
            subset=[c for c in ["Align_fraction_ref", "Align_fraction_query"] if c in df.columns]
        )

    visual_rows = []
    seen_pairs = set()

    def _in_jupyter():
        try:
            from IPython import get_ipython

            shell = get_ipython().__class__.__name__
            return shell == "ZMQInteractiveShell"
        except Exception:
            return False

    total = len(df)
    iterator = None
    if _in_jupyter():
        try:
            from tqdm.notebook import tqdm

            iterator = tqdm(df.itertuples(include_header=False), total=total)
        except Exception:
            iterator = df.itertuples(include_header=False)
    else:
        try:
            from rich.progress import Progress

            progress = Progress()
            task = progress.add_task("fastANI visualize", total=total)
            progress.start()

            def iterator_gen():
                for r in df.itertuples(include_header=False):
                    yield r
                    progress.advance(task)

            iterator = iterator_gen()
        except Exception:
            iterator = df.itertuples(include_header=False)

    for row in iterator:
        try:
            q_name = (
                getattr(row, q_col, None) if hasattr(row, q_col) else row[df.columns.get_loc(q_col)]
            )
        except Exception:
            q_name = row[df.columns.get_loc(q_col)] if q_col in df.columns else None
        try:
            r_name = (
                getattr(row, r_col, None) if hasattr(row, r_col) else row[df.columns.get_loc(r_col)]
            )
        except Exception:
            r_name = row[df.columns.get_loc(r_col)] if r_col in df.columns else None

        if q_name is None or r_name is None:
            continue

        pair_key = tuple(sorted([str(q_name), str(r_name)]))
        if pair_key in seen_pairs:
            continue
        seen_pairs.add(pair_key)

        raw0 = str(pair_key[0])
        raw1 = str(pair_key[1])
        out_file = work_dir / f"{raw0}__vs__{raw1}.visual"
        out_file.parent.mkdir(parents=True, exist_ok=True)

        if out_file.exists():
            visual_file = out_file
        else:
            try:
                q_path = _resolve_path(pair_key[0], out)
                r_path = _resolve_path(pair_key[1], out)
            except FileNotFoundError as exc:
                console.log(f"Skipping pair {pair_key}: {exc}")
                continue

            temp_base = work_dir / f"{raw0}__vs__{raw1}.fastani"
            temp_base.parent.mkdir(parents=True, exist_ok=True)
            temp_out = temp_base

            cmd = [
                "fastANI",
                "-q",
                str(q_path),
                "-r",
                str(r_path),
                "--visualize",
                "-o",
                str(temp_out),
                "-t",
                str(threads),
            ]
            temp_log = Path(str(temp_out) + ".fastani.log")
            try:
                with open(temp_log, "w") as logfh:
                    subprocess.run(cmd, check=True, stdout=logfh, stderr=logfh, text=True)
            except subprocess.CalledProcessError as e:
                console.log(f"fastANI failed for pair {pair_key}: {e}; see {temp_log}")
                continue

            visual_file = Path(str(temp_out) + ".visual")
            if not visual_file.exists():
                console.log(
                    f"fastANI did not create .visual for {pair_key}, expected {visual_file}"
                )
                continue

            try:
                visual_file.rename(out_file)
                visual_file = out_file
            except Exception:
                pass

            if not keep_temp:
                for ext in [".frag", ".fastani", ".log"]:
                    p = Path(str(temp_out) + ext)
                    if p.exists():
                        with contextlib.suppress(Exception):
                            p.unlink()

        try:
            parsed = pl.read_csv(
                visual_file,
                separator="\t",
                header=None,
                names=[
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
                dtype={"query": str, "ref": str},
            )
        except Exception as e:
            console.log(f"Failed to parse visual file {visual_file}: {e}")
            continue

        parsed = parsed[["query", "ref", "ani", "query_start", "query_end", "ref_start", "ref_end"]]

        parsed["query"] = parsed["query"].apply(lambda s: Path(str(s)).name if pl.notna(s) else s)
        parsed["ref"] = parsed["ref"].apply(lambda s: Path(str(s)).name if pl.notna(s) else s)

        if all_neigh is not None and all_neigh.height > 0:
            start_map = {}
            id_map = {}
            for nr in all_neigh.iter_rows(named=True):
                temp = None
                if "temp_seqid" in nr and pl.notna(nr.get("temp_seqid")):
                    temp = str(nr.get("temp_seqid"))
                if "seqid" in nr and pl.notna(nr.get("seqid")):
                    seqid = str(nr.get("seqid"))
                else:
                    seqid = temp

                start_win = (
                    int(nr.get("start_win"))
                    if "start_win" in nr and pl.notna(nr.get("start_win"))
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
                if name is None:
                    return 0
                name = str(name)
                return (
                    start_map.get(name)
                    or start_map.get(Path(name).name)
                    or start_map.get(Path(name).stem)
                    or 0
                )

            def _map_to_seqid(name: str) -> str:
                if name is None:
                    return name
                name = str(name)
                return (
                    id_map.get(name)
                    or id_map.get(Path(name).name)
                    or id_map.get(Path(name).stem)
                    or name
                )

            parsed["q_offset"] = parsed["query"].apply(_find_offset)
            parsed["r_offset"] = parsed["ref"].apply(_find_offset)

            parsed["query_start"] = parsed["query_start"].astype(float) + parsed["q_offset"]
            parsed["query_end"] = parsed["query_end"].astype(float) + parsed["q_offset"]
            parsed["ref_start"] = parsed["ref_start"].astype(float) + parsed["r_offset"]
            parsed["ref_end"] = parsed["ref_end"].astype(float) + parsed["r_offset"]

            parsed["query"] = parsed["query"].apply(_map_to_seqid)
            parsed["ref"] = parsed["ref"].apply(_map_to_seqid)

            parsed = parsed.drop(["q_offset", "r_offset"])

        visual_rows.append(parsed)

    if visual_rows:
        visual_df = pl.concat(visual_rows, how="vertical")
    else:
        visual_df = pl.DataFrame(
            columns=["query", "ref", "ani", "query_start", "query_end", "ref_start", "ref_end"]
        )

    return visual_df
