#!/usr/bin/env python3
# fragment_generation.py
# Copyright (c) 2026 Eliška Korbová ORCID 0009-0004-1247-0808
#
# Library fragment generation - with cluster and hotspot competition logic
# Output: fillered.csv

from __future__ import annotations
import argparse
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np
import pandas as pd

_GLOBAL_EDICT: dict[str, tuple[int, int]] | None = None

#loading enzyme combinations one per script
def init_worker(edict_path: str) -> None:
    global _GLOBAL_EDICT
    _GLOBAL_EDICT = load_enzyme_dict(Path(edict_path))

#parsing size-range
def parse_size_range(s: str) -> tuple[int, int]:
    parts = s.split("-")
    if len(parts) != 2:
        sys.exit(f"[ERROR] Invalid size range '{s}'. Use LOW-HIGH format e.g. 150-350")
    try:
        low, high = int(parts[0]), int(parts[1])
    except ValueError:
        sys.exit(f"[ERROR] Size range must have whole numbers, not '{s}'")
    if low >= high:
        sys.exit(f"[ERROR] LOW ({low}) must be less than HIGH ({high})")
    return low, high


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="fragment_generation.py",
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    p.add_argument("--workdir", required=True)
    p.add_argument("--size", required=True)
    p.add_argument("--parallel", type=int, default=2)
    return p.parse_args()

def load_enzyme_dict(path: Path) -> dict[str, tuple[int, int]]:
    df = pd.read_csv(path, dtype=str)
    df.columns = df.columns.str.strip()

    required = {"enzyme_name", "expanded_sequence", "cut_offset"}
    missing  = required - set(df.columns)
    if missing:
        sys.exit(f"[ERROR] enzyme_dict missing columns: {missing}")

    df["enzyme_name"]      = df["enzyme_name"].astype(str).str.strip()
    df["expanded_sequence"] = df["expanded_sequence"].astype(str).str.strip()
    df["cut_offset"]       = df["cut_offset"].astype(int)
    df["motif_length"]     = df["expanded_sequence"].str.len()

    out: dict[str, tuple[int, int]] = {}
    for enz, _seq, cutoff, mlen in df[
        ["enzyme_name", "expanded_sequence", "cut_offset", "motif_length"]
    ].itertuples(index=False, name=None):
        enz    = str(enz)
        cutoff = int(cutoff)
        mlen   = int(mlen)
        if (enz not in out) or (mlen > out[enz][1]):
            out[enz] = (cutoff, mlen)

    return out

# Fragment geometry helpers - Vektorizovaný výpočet motif_start / motif_end pro právě dva enzymy
# vraci - (motif_start, motif_end) jako int64 arrays stejné délky jako cut_positions
def compute_pair_intervals(
    cut_positions: np.ndarray,
    enzymes: np.ndarray,
    enzyme1: str,
    enzyme2: str,
    edict: dict[str, tuple[int, int]],
) -> tuple[np.ndarray, np.ndarray]:

    if enzyme1 not in edict:
        raise KeyError(f"'{enzyme1}' not in enzyme dictionary")
    if enzyme2 not in edict:
        raise KeyError(f"'{enzyme2}' not in enzyme dictionary")

    motif_start = np.empty(len(cut_positions), dtype=np.int64)
    motif_end   = np.empty(len(cut_positions), dtype=np.int64)

    for enz in (enzyme1, enzyme2):
        offset, motif_len = edict[enz]
        mask = enzymes == enz
        motif_start[mask] = cut_positions[mask] - int(offset)
        motif_end[mask]   = cut_positions[mask] + (int(motif_len) - int(offset))

    return motif_start, motif_end


def intervals_overlap(
    s1: np.ndarray | int,
    e1: np.ndarray | int,
    s2: np.ndarray | int,
    e2: np.ndarray | int,
) -> np.ndarray | bool:
    return np.maximum(s1, s2) < np.minimum(e1, e2)  # atleast 1bp overlap not edges

def size_ok(
    fragment_length: np.ndarray | int,
    min_size: int,
    max_size: int,
) -> np.ndarray | bool:
    return (fragment_length >= min_size) & (fragment_length <= max_size)  #true if in closed interval mix and max size arg

# start and end indices of blocks of same accessions
def accession_blocks(
    accessions: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:

    n = len(accessions)
    if n == 0:
        return np.array([], dtype=np.int64), np.array([], dtype=np.int64)

    changes = np.flatnonzero(accessions[1:] != accessions[:-1]) + 1
    starts  = np.concatenate(([0],  changes))
    ends    = np.concatenate((changes, [n]))
    return starts.astype(np.int64), ends.astype(np.int64)

#block processing per-accession
def process_accession_block(
    accession: str,
    pos: np.ndarray,
    enz: np.ndarray,
    mstart: np.ndarray,
    mend: np.ndarray,
    min_size: int,
    max_size: int,
) -> tuple[
    list[tuple[str, int, str, str, int]],   # filtered fragments
    list[tuple[int, int]],                  # cluster spans (local indices)
]:
    # for one accession block
    # detection of overlaps clusters/hotspots
    # emits fragments outside of clusters
    # emits boundary candidates for clusters
    # for one accession block
    # detects overlap clusters/hotspots
    # emits fragments outside clusters
    # for each cluster boundary tries up to 2 innermost cuts per side,
    # writes only the first that passes size filter (no competition)
    filtered: list[tuple[str, int, str, str, int]] = []
    clusters: list[tuple[int, int]] = []

    n = len(pos)
    if n < 2:
        return filtered, clusters

    overlap_adj = intervals_overlap(mstart[:-1], mend[:-1], mstart[1:], mend[1:])

    in_cluster = np.zeros(n, dtype=bool)
    if np.any(overlap_adj):
        ov      = overlap_adj.astype(np.int8, copy=False)
        padded  = np.pad(ov, (1, 1), constant_values=0)
        changes = np.diff(padded)

        run_starts = np.flatnonzero(changes ==  1)
        run_ends   = np.flatnonzero(changes == -1) - 1  # in coordinates of overlap_adj

        for rs, re in zip(run_starts, run_ends):
            cs = int(rs)
            ce = int(re + 1)  # from overlap run to row span
            clusters.append((cs, ce))
            in_cluster[cs : ce + 1] = True

    # normal fragments besides clusters/hotspots
    normal_mask = (
        (~in_cluster[:-1]) &
        (~in_cluster[1:])  &
        (enz[:-1] != enz[1:]) &
        (~overlap_adj)
    )
    idx = np.flatnonzero(normal_mask)
    if idx.size:
        frag_len = pos[idx + 1] - pos[idx]
        for k, i in enumerate(idx):
            if size_ok(int(frag_len[k]), min_size, max_size):
                j = i + 1
                filtered.append((accession, int(pos[i]), str(enz[i]), str(enz[j]), int(frag_len[k])))

    # borderline cluster fragments
    for cs, ce in clusters:
        # left borderline first that gets into library is written
        for i, j in [(cs, cs + 1), (cs + 1, cs + 2)]:
            if j > ce:
                break  
            if enz[i] != enz[j]:
                frag_len = int(pos[j] - pos[i])
                if size_ok(frag_len, min_size, max_size):
                    filtered.append((accession, int(pos[i]), str(enz[i]), str(enz[j]), frag_len))
                    break 

        # right borderline first that gets into library is written
        for i, j in [(ce - 1, ce), (ce - 2, ce - 1)]:
            if i < cs:
                break  
            if enz[i] != enz[j]:
                frag_len = int(pos[j] - pos[i])
                if size_ok(frag_len, min_size, max_size):
                    filtered.append((accession, int(pos[i]), str(enz[i]), str(enz[j]), frag_len))
                    break  

    return filtered, clusters

#job processing
def _fail(combo: str, error: str) -> dict:
    return {
        "combo":     combo,
        "cuts":      0,
        "uncertain": 0,
        "filtered":  0,
        "clusters":  0,
        "ok":        False,
        "error":     error,
    }

# worker entry point - one cuts.csv, returns summary dict
def process_job(
    job: dict,
    min_size: int,
    max_size: int,
) -> dict:
    
    global _GLOBAL_EDICT

    try:
        if _GLOBAL_EDICT is None:
            return _fail(job["combo"], "Enzyme dictionary not loaded in worker.")

        combo      = job["combo"]
        enzyme1    = job["enzyme1"]
        enzyme2    = job["enzyme2"]
        input_path = Path(job["input_path"])
        outdir     = Path(job["outdir"])

        outdir.mkdir(parents=True, exist_ok=True)

        if not input_path.exists():
            return _fail(combo, f"cuts.csv not found: {input_path}")

        try:
            df = pd.read_csv(
                input_path,
                usecols=["accession", "cut_position", "enzyme"],
                dtype={"accession": "string", "cut_position": "int64", "enzyme": "string"},
                low_memory=False,
            )
        except Exception as e:
            return _fail(combo, f"Failed to read CSV: {e}")

        if df.empty:
            return _fail(combo, f"No rows in {input_path}")

        df = df[df["enzyme"].isin({enzyme1, enzyme2})].copy()
        if df.empty:
            return _fail(combo, f"No cuts for {enzyme1}/{enzyme2} in {input_path}")

        observed = set(df["enzyme"].unique().tolist())
        if enzyme1 not in observed or enzyme2 not in observed:
            return _fail(combo, f"One of the requested enzymes is missing in {input_path}")

        df.sort_values(
            ["accession", "cut_position"],
            inplace=True,
            kind="mergesort",
            ignore_index=True,
        )

        accessions    = df["accession"].to_numpy(dtype=object,    copy=False)
        cut_positions = df["cut_position"].to_numpy(dtype=np.int64, copy=False)
        enzymes       = df["enzyme"].to_numpy(dtype=object,       copy=False)

        edict = _GLOBAL_EDICT
        assert edict is not None

        try:
            motif_start, motif_end = compute_pair_intervals(
                cut_positions=cut_positions,
                enzymes=enzymes,
                enzyme1=enzyme1,
                enzyme2=enzyme2,
                edict=edict,
            )
        except KeyError as e:
            return _fail(combo, str(e))

        starts, ends = accession_blocks(accessions)

        filt_fragments: list[tuple[str, int, str, str, int]] = []
        uncertain_mask = np.zeros(len(df), dtype=bool)
        n_clusters     = 0

        for s, e in zip(starts, ends):
            block_acc  = str(accessions[s])
            block_pos  = cut_positions[s:e]
            block_enz  = enzymes[s:e]
            block_ms   = motif_start[s:e]
            block_me   = motif_end[s:e]

            block_filt, clusters = process_accession_block(
                accession=block_acc,
                pos=block_pos,
                enz=block_enz,
                mstart=block_ms,
                mend=block_me,
                min_size=min_size,
                max_size=max_size,
            )

            filt_fragments.extend(block_filt)

            if clusters:
                n_clusters += len(clusters)
                for cs, ce in clusters:
                    uncertain_mask[s + cs : s + ce + 1] = True

        _FRAG_COLS = ["accession", "start_pos", "start_enzyme", "end_enzyme", "fragment_length"]

        filt_df = pd.DataFrame(filt_fragments, columns=_FRAG_COLS)
        if not filt_df.empty:
            filt_df.drop_duplicates(inplace=True)
            filt_df.sort_values(
                ["accession", "start_pos", "fragment_length", "end_enzyme"],
                inplace=True,
                kind="mergesort",
                ignore_index=True,
            )

        uncertain_df = (
            df.loc[uncertain_mask, ["accession", "cut_position", "enzyme"]]
            .copy()
            .sort_values(["accession", "cut_position"], kind="mergesort", ignore_index=True)
        )

        # output files
        filt_df.to_csv(outdir / "filtered.csv",  index=False)
        uncertain_df.to_csv(outdir / "uncertain_cuts.csv", index=False)

        total       = len(df)
        n_uncertain = len(uncertain_df)
        n_filt      = len(filt_df)
        pct         = (n_uncertain / total * 100.0) if total > 0 else 0.0

        stats = pd.DataFrame([{
            "total_cuts":         total,
            "uncertain_cuts":     n_uncertain,
            "pct_uncertain":      round(pct, 4),
            "filtered_fragments": n_filt,
            "uncertain_clusters": n_clusters,
        }])
        stats.to_csv(outdir / "statistics_uncertain.csv", index=False)

        return {
            "combo":     combo,
            "cuts":      total,
            "uncertain": n_uncertain,
            "filtered":  n_filt,
            "clusters":  n_clusters,
            "ok":        True,
            "error":     None,
        }

    except Exception as e:
        combo = job.get("combo", "unknown")
        return _fail(combo, f"Unexpected error: {e}")

# combinations z combinations.csv
def load_combinations(path: Path) -> list[tuple[str, str]]:
    
    df = pd.read_csv(path, dtype=str)
    df.columns = df.columns.str.strip()

    if {"enzyme_a", "enzyme_b"}.issubset(df.columns):
        pairs = list(zip(df["enzyme_a"].astype(str), df["enzyme_b"].astype(str)))

    elif {"enzyme1", "enzyme2"}.issubset(df.columns):
        pairs = list(zip(df["enzyme1"].astype(str), df["enzyme2"].astype(str)))

    elif "combo" in df.columns:
        pairs = []
        for val in df["combo"].astype(str):
            parts = val.split("_")
            if len(parts) != 2:
                raise ValueError(f"Cannot parse combination '{val}' in combo column.")
            pairs.append((parts[0], parts[1]))

    elif "combination" in df.columns:
        pairs = []
        for val in df["combination"].astype(str):
            parts = val.split("_")
            if len(parts) != 2:
                raise ValueError(f"Cannot parse combination '{val}' in combination column.")
            pairs.append((parts[0], parts[1]))

    else:
        if df.shape[1] < 2:
            raise ValueError(
                "Combinations file must have at least two columns "
                "(enzyme1, enzyme2) or a combo/combination column."
            )
        pairs = list(zip(df.iloc[:, 0].astype(str), df.iloc[:, 1].astype(str)))

    seen: set[tuple[str, str]] = set()
    unique_pairs: list[tuple[str, str]] = []
    for a, b in pairs:
        if (a, b) not in seen:
            seen.add((a, b))
            unique_pairs.append((a, b))

    return unique_pairs


def make_job(
    cuts_root: Path,
    outdir_root: Path,
    cuts_file: str,
    enzyme1: str,
    enzyme2: str,
) -> dict:
    combo = f"{enzyme1}_{enzyme2}"
    return {
        "combo":      combo,
        "enzyme1":    enzyme1,
        "enzyme2":    enzyme2,
        "input_path": cuts_root / combo / cuts_file,
        "outdir":     outdir_root / combo,
    }

# batch runner
def run_batch(
    combinations_path: Path,
    cuts_root: Path,
    outdir_root: Path,
    cuts_file: str,
    enzyme_dict: Path,
    parallel: int,
    min_size: int,
    max_size: int,
) -> None:
    pairs = load_combinations(combinations_path)
    if not pairs:
        sys.exit("[ERROR] No combinations found.")

    jobs = [
        make_job(
            cuts_root=cuts_root,
            outdir_root=outdir_root,
            cuts_file=cuts_file,
            enzyme1=a,
            enzyme2=b,
        )
        for a, b in pairs
    ]

    outdir_root.mkdir(parents=True, exist_ok=True)

    print(f"[batch] Combinations : {len(jobs):,}")
    print(f"[batch] Workers      : {parallel}")
    print(f"[batch] Size window  : {min_size}-{max_size} bp")
    print(f"[batch] Cuts root    : {cuts_root}")
    print(f"[batch] Outdir root  : {outdir_root}")

    results: list[dict] = []

    with ProcessPoolExecutor(
        max_workers=max(1, parallel),
        initializer=init_worker,
        initargs=(str(enzyme_dict),),
    ) as ex:
        futures = [
            ex.submit(process_job, job, min_size, max_size)
            for job in jobs
        ]

        for fut in as_completed(futures):
            res = fut.result()
            results.append(res)

            if res["ok"]:
                print(
                    f"[done] {res['combo']}: "
                    f"cuts={res['cuts']:,}, "
                    f"uncertain={res['uncertain']:,}, "
                    f"filtered={res['filtered']:,}, "
                    f"clusters={res['clusters']:,}"
                )
            else:
                print(f"[fail] {res['combo']}: {res['error']}")

    summary = pd.DataFrame(results)
    summary.to_csv(outdir_root / "batch_summary.csv", index=False)
    print(f"[batch] Summary written → {outdir_root / 'batch_summary.csv'}")


def main() -> None:
    args        = parse_args()
    min_size, max_size = parse_size_range(args.size)

    script_dir = Path(__file__).resolve().parent
    enzyme_dict = script_dir / "enzymes.csv"

    if not enzyme_dict.exists():
        sys.exit(f"[ERROR] Enzyme dict not found: {enzyme_dict}")

    results_dir = Path(args.workdir) / "results"
    combinations_path = results_dir / "combinations.csv"
    cuts_root = results_dir
    outdir_root = results_dir

    run_batch(
        combinations_path=combinations_path,
        cuts_root=cuts_root,
        outdir_root=outdir_root,
        cuts_file="cuts.csv",
        enzyme_dict=enzyme_dict,
        parallel=args.parallel,
        min_size=min_size,
        max_size=max_size,
    )


if __name__ == "__main__":
    main()