#!/usr/bin/env python3
# postprocess_metrics.py
# Copyright (c) 2026 Eliška Korbová ORCID 0009-0004-1247-0808
#
# Filltered fragments for ddRAD library GC content, size distribution of all fragments
# Output: distribution.csv, gc_metrics.csv, contig_lengths.txt

from __future__ import annotations
import argparse
import csv
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any
import numpy as np
import pyfaidx

FragmentRow = tuple[str, int, str, str, int]

_STD_LABELS = [f"{b}-{b + 99}" for b in range(0, 1000, 100)]

_WORKER_REF: pyfaidx.Fasta | None = None


def parse_size_range(size_str: str) -> tuple[int, int]:
    parts = size_str.split("-")
    if len(parts) != 2:
        raise SystemExit(f"[ERROR] Invalid size range '{size_str}'. Use LOW-HIGH, e.g. 200-400")
    try:
        low, high = int(parts[0]), int(parts[1])
    except ValueError:
        raise SystemExit(f"[ERROR] Size range must have whole numbers, not '{size_str}'")
    if low >= high:
        raise SystemExit(f"[ERROR] LOW ({low}) must be less than HIGH ({high})")
    return low, high


def load_combinations_csv(path: Path) -> list[str]:
    with open(path, "r", encoding="utf-8", newline="") as fh:
        return [f"{r['enzyme_a']}_{r['enzyme_b']}" for r in csv.DictReader(fh)]


def read_fragments_csv(path: Path) -> list[FragmentRow]:
    with open(path, "r", encoding="utf-8", newline="") as fh:
        return [
            (
                r["accession"],
                int(r["start_pos"]),
                r["start_enzyme"],
                r["end_enzyme"],
                int(r["fragment_length"]),
            )
            for r in csv.DictReader(fh)
        ]


def ensure_fai_index(ref_path: Path, results_dir: Path) -> Path:
    fai = results_dir / (ref_path.name + ".fai")  # index goes into results/, not next to FASTA
    if not fai.exists():
        print(f"  [INFO] Creating .fai index -> {fai}")
        print(f"         (one-time only, may take a few minutes for large genomes)")
        pyfaidx.Fasta(str(ref_path), indexname=str(fai))  # creates index and closes
        print(f"  [INFO] Index created")
    else:
        print(f"  [INFO] .fai index found: {fai}")
    return fai


def write_contig_lengths(ref: pyfaidx.Fasta, results_dir: Path) -> Path:
    out = results_dir / "contig_lengths.txt"
    entries = sorted(
        ((acc, len(ref[acc])) for acc in ref.keys()),
        key=lambda x: x[1],
        reverse=True,  # longest first — consistent with plots.py --chroms N
    )
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("accession\tlength_bp\n")
        for acc, length in entries:
            fh.write(f"{acc}\t{length}\n")
    print(f"  [INFO] contig_lengths.txt -> {out}")
    return out


def std_size_label(fl: int) -> str:
    if fl >= 1000:
        return ">=1000"
    b = (fl // 100) * 100
    return f"{b}-{b + 99}"


def label_overlaps_range(label: str, low: int, high: int) -> bool:
    if label == ">=1000":
        return high >= 1000
    a, b = map(int, label.split("-"))
    return a <= high and b >= low


def build_distribution_rows(
    fragments: list[FragmentRow],
    filtered_count: int,
    low: int,
    high: int,
) -> list[list[Any]]:
    counts: dict[str, int] = defaultdict(int)
    for *_, fl in fragments:
        counts[std_size_label(fl)] += 1

    total = sum(counts.values()) or 1
    custom = f"custom_{low}-{high}"

    rows = []
    for lab in _STD_LABELS + [">=1000"]:
        cnt = counts.get(lab, 0)
        rows.append([
            lab,
            cnt,
            f"{100 * cnt / total:.2f}",
            "selected" if label_overlaps_range(lab, low, high) else "",
        ])

    rows.append([custom, filtered_count, f"{100 * filtered_count / total:.2f}", "user_window"])
    return rows


def write_distribution_csv(combo_dir: Path, rows: list[list[Any]]) -> Path:
    out = combo_dir / "distribution.csv"
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["length_range", "count", "percentage", "note"])
        w.writerows(rows)
    return out


def gc_content(seq: str) -> float:
    if not seq:
        return 0.0
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / len(seq)


def compute_gc_stats(filtered: list[FragmentRow]) -> dict[str, Any]:
    vals: list[float] = []

    for acc, start, _, _, fl in filtered:
        try:
            seq = str(_WORKER_REF[acc][start:start + fl])  # type: ignore[index]
        except Exception:
            continue
        if seq:
            vals.append(gc_content(seq))

    if not vals:
        return {
            "n_fragments":   0,
            "gc_mean_pct":   "n/a",
            "gc_median_pct": "n/a",
            "gc_min_pct":    "n/a",
            "gc_max_pct":    "n/a",
            "gc_std_pct":    "n/a",
        }

    arr = np.array(vals)
    pct = lambda x: f"{x * 100:.2f}"

    return {
        "n_fragments":   len(arr),
        "gc_mean_pct":   pct(np.mean(arr)),
        "gc_median_pct": pct(np.median(arr)),
        "gc_min_pct":    pct(np.min(arr)),
        "gc_max_pct":    pct(np.max(arr)),
        "gc_std_pct":    pct(np.std(arr)),
    }


def write_gc_csv(combo_dir: Path, stats: dict[str, Any]) -> Path:
    out = combo_dir / "gc_metrics.csv"
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["metric", "value"])
        w.writerows(stats.items())
    return out


def init_worker(ref_path: str, index_path: str) -> None:
    global _WORKER_REF
    _WORKER_REF = pyfaidx.Fasta(ref_path, indexname=index_path, rebuild=False)


def process_combination(
    combo: str,
    results_dir: str,
    low: int,
    high: int,
) -> dict[str, Any]:
    combo_dir = Path(results_dir) / combo

    if not combo_dir.is_dir():
        raise FileNotFoundError(f"Missing combination directory: {combo_dir}")
    if not (combo_dir / "fragments.csv").exists():
        raise FileNotFoundError(f"Missing fragments.csv: {combo_dir / 'fragments.csv'}")
    if not (combo_dir / "filtered.csv").exists():
        raise FileNotFoundError(f"Missing filtered.csv: {combo_dir / 'filtered.csv'}")

    fragments = read_fragments_csv(combo_dir / "fragments.csv")
    filtered  = read_fragments_csv(combo_dir / "filtered.csv")

    write_distribution_csv(
        combo_dir,
        build_distribution_rows(fragments, len(filtered), low, high),
    )

    gc_stats = compute_gc_stats(filtered)
    write_gc_csv(combo_dir, gc_stats)

    return {
        "combo":       combo,
        "fragments":   len(fragments),
        "filtered":    len(filtered),
        "gc_mean_pct": gc_stats["gc_mean_pct"],
    }

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--workdir",  required=True)
    p.add_argument("--ref",      required=True)
    p.add_argument("--parallel", type=int, default=2)
    p.add_argument("--size",     required=True)
    return p.parse_args()

def main() -> None:
    args = parse_args()

    if args.parallel < 1:
        raise SystemExit("[ERROR] --parallel must be >= 1")

    try:
        low, high = parse_size_range(args.size)
    except SystemExit:
        raise

    workdir     = Path(args.workdir)
    results_dir = workdir / "results"
    ref_path    = Path(args.ref)

    if not results_dir.exists():
        raise SystemExit(f"[ERROR] results directory not found: {results_dir}")
    if not ref_path.exists():
        raise SystemExit(f"[ERROR] reference FASTA not found: {ref_path}")

    # build index into results/
    fai = ensure_fai_index(ref_path, results_dir)
    index_path = str(fai)

    # open reference in main process for contig_lengths.txt
    print(f"\n[1/4] Building contig lengths table ...")
    ref = pyfaidx.Fasta(str(ref_path), indexname=index_path, rebuild=False)
    write_contig_lengths(ref, results_dir)
    print(f"      {len(ref.keys())} contigs")

    combinations = load_combinations_csv(results_dir / "combinations.csv")
    print(f"\n[2/4] Loading combinations: {results_dir / 'combinations.csv'}")
    print(f"      {len(combinations)} combination(s) found")

    for combo in combinations:
        if not (results_dir / combo).is_dir():
            raise SystemExit(f"[ERROR] Missing combination directory: {results_dir / combo}")

    print(f"\n[3/4] Processing {len(combinations)} combination(s) | parallel={args.parallel} | size={low}-{high}\n")

    results: list[dict[str, Any]] = []
    with ProcessPoolExecutor(
        max_workers=args.parallel,
        initializer=init_worker,
        initargs=(str(ref_path), index_path),  # index path passed to every worker
    ) as ex:
        futures = {
            ex.submit(process_combination, c, str(results_dir), low, high): c
            for c in combinations
        }

        for idx, future in enumerate(as_completed(futures), start=1):
            combo = futures[future]
            try:
                r = future.result()
                results.append(r)
                print(
                    f"  [{idx:>3}/{len(futures)}] {combo} "
                    f"fragments={r['fragments']:,} "
                    f"filtered={r['filtered']:,} "
                    f"gc_mean={r['gc_mean_pct']}"
                )
            except Exception as exc:
                print(f"  [ERROR] {combo}: {exc}")
                raise

    print(f"\n[4/4] Done. Processed={len(results)}")


if __name__ == "__main__":
    main()