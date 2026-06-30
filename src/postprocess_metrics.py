#!/usr/bin/env python3
# postprocess_metrics.py
# Copyright (c) 2026 Eliška Korbová ORCID 0009-0004-1247-0808
#
# Filltered fragments for ddRAD library GC content, size distribution of all fragments
# Output: distribution.csv, gc_metrics.csv, contig_lengths.txt

from __future__ import annotations
import argparse
import csv
import gzip
import sys
import numpy as np
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any

FragmentRow = tuple[str, int, str, str, int]

_STD_LABELS = [f"{b}-{b + 99}" for b in range(0, 1000, 100)]

_WORKER_REF_SEQS: dict[str, str] | None = None

def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r", encoding="utf-8")


def gc_content(seq: str) -> float:
    if not seq:
        return 0.0
    seq = seq.upper()
    return (seq.count("G") + seq.count("C")) / len(seq)


def read_fasta_contigs(path: Path) -> list[tuple[str, str]]:
    contigs: list[tuple[str, str]] = []
    header: str | None = None
    parts: list[str] = []

    with open_text(path) as fh:
        for raw in fh:
            line = raw.rstrip("\n").rstrip("\r")
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    contigs.append((header, "".join(parts).upper()))
                raw_header = line[1:].strip()
                if not raw_header:
                    raise SystemExit(f"[ERROR] Empty FASTA header: {line!r}")
                header = raw_header.split()[0]
                parts = []
            else:
                parts.append("".join(line.split()))

    if header is not None:
        contigs.append((header, "".join(parts).upper()))
    elif not contigs:
        raise SystemExit(f"[ERROR] FASTA is empty or invalid: {path}")

    return contigs


def sort_contigs_by_length(contigs: list[tuple[str, str]]) -> list[tuple[str, str]]:
    return sorted(contigs, key=lambda x: len(x[1]), reverse=True)


def write_contig_lengths_txt(contigs_sorted: list[tuple[str, str]], results_dir: Path) -> Path:
    out_path = results_dir / "contig_lengths.txt"
    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        fh.write("accession\tlength_bp\n")
        for acc, seq in contigs_sorted:
            fh.write(f"{acc}\t{len(seq)}\n")
    print(f"  [INFO] contig_lengths.txt -> {out_path}")
    return out_path


def parse_size_range(size_str: str) -> tuple[int, int]:
    parts = size_str.split("-")
    if len(parts) != 2:
        raise ValueError(f"Invalid size range '{size_str}'. Use LOW-HIGH, e.g. 200-400")
    low, high = int(parts[0]), int(parts[1])
    if low >= high:
        raise ValueError(f"LOW ({low}) must be smaller than HIGH ({high})")
    return low, high


def load_combinations_csv(path: Path) -> list[str]:
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        combos = [f"{row['enzyme_a']}_{row['enzyme_b']}" for row in reader]
    return combos

# Je prý ok řešení, ale mám jako scaling point zmínit, že list je náročný na RAM, takže by tě možná mohla do budoucna zajímat streamovací varianta (csv.DictReader)  
def read_fragments_csv(path: Path) -> list[FragmentRow]:
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        rows: list[FragmentRow] = [] 
        for row in reader:
            rows.append(
                (
                    row["accession"],
                    int(row["start_pos"]),
                    row["start_enzyme"],
                    row["end_enzyme"],
                    int(row["fragment_length"]),
                )
            )
    return rows


def std_size_label(fl: int) -> str:
    if fl >= 1000:
        return ">=1000"
    start = (fl // 100) * 100
    return f"{start}-{start + 99}"


def label_overlaps_range(label: str, low: int, high: int) -> bool:
    if label == ">=1000":
        return high >= 1000
    b_lo, b_hi = (int(x) for x in label.split("-"))
    return b_lo <= high and b_hi >= low


def build_distribution_rows(fragments: list[FragmentRow], filtered_count: int, size_low: int, size_high: int) -> list[list[Any]]:
    counts = defaultdict(int)
    for _, _, _, _, fl in fragments:
        counts[std_size_label(fl)] += 1

    total_all = sum(counts.values()) or 1
    custom_label = f"custom_{size_low}-{size_high}"

    rows: list[list[Any]] = []
    for lab in _STD_LABELS + [">=1000"]:
        cnt = counts.get(lab, 0)
        note = "selected" if label_overlaps_range(lab, size_low, size_high) else ""
        rows.append([lab, cnt, f"{100 * cnt / total_all:.2f}", note])

    rows.append([custom_label, filtered_count, f"{100 * filtered_count / total_all:.2f}", "user_window"])
    return rows


def write_distribution_csv(combo_dir: Path, rows: list[list[Any]]) -> Path:
    out_path = combo_dir / "distribution.csv"
    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["length_range", "count", "percentage", "note"])
        w.writerows(rows)
    return out_path


def compute_gc_stats(filtered: list[FragmentRow], ref_seqs: dict[str, str]) -> dict[str, Any]:
    gc_vals: list[float] = []

    for acc, start, _, _, fl in filtered:
        frag_seq = ref_seqs.get(acc, "")[start : start + fl]
        if frag_seq:
            gc_vals.append(gc_content(frag_seq))
            # JC: Může se stát, že start + fl > len(seq)? Pokud ano, jak se pak počítá GC?

    if not gc_vals:
        return {
            "n_fragments":   0,
            "gc_mean_pct":   "n/a",
            "gc_median_pct": "n/a",
            "gc_min_pct":    "n/a",
            "gc_max_pct":    "n/a",
            "gc_std_pct":    "n/a",
        }

    gc_array = np.array(gc_vals)
    pct = lambda v: f"{v * 100:.2f}"

    return {
        "n_fragments":   len(gc_array),
        "gc_mean_pct":   pct(np.mean(gc_array)),
        "gc_median_pct": pct(np.median(gc_array)),
        "gc_min_pct":    pct(np.min(gc_array)),
        "gc_max_pct":    pct(np.max(gc_array)),
        "gc_std_pct":    pct(np.std(gc_array)),
    }


def write_gc_csv(combo_dir: Path, gc_stats: dict[str, Any]) -> Path:
    out_path = combo_dir / "gc_metrics.csv"
    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["metric", "value"])
        for k, v in gc_stats.items():
            w.writerow([k, v])
    return out_path


def init_worker(ref_path: str) -> None:
    global _WORKER_REF_SEQS
    _WORKER_REF_SEQS = {acc: seq for acc, seq in read_fasta_contigs(Path(ref_path))}


def process_combination(combo_name: str, results_dir: str, size_low: int, size_high: int) -> dict[str, Any]:
    global _WORKER_REF_SEQS

    results_path = Path(results_dir)
    combo_dir = results_path / combo_name
    fragments_path = combo_dir / "fragments.csv"
    filtered_path = combo_dir / "filtered.csv"

    if not combo_dir.is_dir():
        raise FileNotFoundError(f"Missing combination directory: {combo_dir}")
    if not fragments_path.exists():
        raise FileNotFoundError(f"Missing fragments.csv: {fragments_path}")
    if not filtered_path.exists():
        raise FileNotFoundError(f"Missing filtered.csv: {filtered_path}")

    fragments = read_fragments_csv(fragments_path)
    filtered = read_fragments_csv(filtered_path)

    dist_rows = build_distribution_rows(
        fragments=fragments,
        filtered_count=len(filtered),
        size_low=size_low,
        size_high=size_high,
    )
    write_distribution_csv(combo_dir, dist_rows)

    if _WORKER_REF_SEQS is None:
        raise RuntimeError("Reference sequences were not initialized in worker")

    gc_stats = compute_gc_stats(filtered, _WORKER_REF_SEQS)
    write_gc_csv(combo_dir, gc_stats)

    return {
        "combo": combo_name,
        "fragments": len(fragments),
        "filtered": len(filtered),
        "gc_mean_pct": gc_stats.get("gc_mean_pct", "n/a"),
    }


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="postprocess_metrics.py",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--workdir", required=True)
    p.add_argument("--ref", required=True)
    p.add_argument("--parallel", type=int, default=2)
    p.add_argument("--size", required=True)
    return p.parse_args()


def main() -> None:
    args = parse_args()

    if args.parallel < 1:
        raise SystemExit("[ERROR] --parallel must be >= 1")

    try:
        size_low, size_high = parse_size_range(args.size)
    except ValueError as e:
        raise SystemExit(f"[ERROR] {e}")

    workdir = Path(args.workdir)
    results_dir = workdir / "results"
    ref_path = Path(args.ref)
    combinations_csv = results_dir / "combinations.csv"

    if not results_dir.exists():
        raise SystemExit(f"[ERROR] results directory not found: {results_dir}")
    if not ref_path.exists():
        raise SystemExit(f"[ERROR] reference FASTA not found: {ref_path}")

    print(f"\n[1/4] Loading reference: {ref_path}")
    contigs_sorted = sort_contigs_by_length(read_fasta_contigs(ref_path))
    total_bases = sum(len(seq) for _, seq in contigs_sorted)
    print(f"      {len(contigs_sorted)} contigs | {total_bases:,} bp total")
    write_contig_lengths_txt(contigs_sorted, results_dir)

    print(f"\n[2/4] Loading combinations: {combinations_csv}")
    combo_names = load_combinations_csv(combinations_csv)
    print(f"      {len(combo_names)} combination(s) found")

    for combo_name in combo_names:
        combo_dir = results_dir / combo_name
        if not combo_dir.is_dir():
            raise SystemExit(f"[ERROR] Missing combination directory: {combo_dir}")
            # Určitě k tomu máš důvod, takže spíš zvědavosti: proč v tomto skriptu využíváš if-raise, když jinde bylo if-continue, případně nějaká poznámka? 
            # Jsou ty chybějící věci v této části tak kritické, že se musí přerušit celý run?

    print(f"\n[3/4] Processing {len(combo_names)} combination(s) | parallel={args.parallel} | size={size_low}-{size_high}\n")

    results: list[dict[str, Any]] = []
    # Copilot upozorňuje, že FASTA takto načítá 1x main process a Nx workers, takže se může reálně stát, že bude potřeba 24 GB paměti. Tohle pouze papouškuju, takže zvaž podle sebe. 
    with ProcessPoolExecutor(max_workers=args.parallel, initializer=init_worker, initargs=(str(ref_path),)) as executor:
        futures = {
            executor.submit(process_combination, combo, str(results_dir), size_low, size_high): combo
            for combo in combo_names
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
                print(f"  [ERROR] {combo}: {exc}", file=sys.stderr)
                raise

    print(f"\n[4/4] Done. Processed={len(results)}")


if __name__ == "__main__":
    main()
