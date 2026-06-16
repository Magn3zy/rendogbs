#!/usr/bin/env python3
from __future__ import annotations
import argparse
import csv
import gzip
import shutil
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any

@dataclass(frozen=True)
class AnnotationSource:
    label: str   # např. gff:gene nebo te:LTR_Gypsy
    path: str    # merged BED


FragmentRow = tuple[str, int, str, str, int]
Interval = tuple[int, int]
IntervalDict = dict[str, list[Interval]]

def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r", encoding="utf-8")


def check_bedtools() -> str:
    bedtools = shutil.which("bedtools")
    if bedtools is None:
        raise SystemExit("[ERROR] bedtools not found in PATH")
    return bedtools


def load_combinations_csv(path: Path) -> list[str]:
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        return [f"{r['enzyme_a']}_{r['enzyme_b']}" for r in reader]


def read_filtered_csv(path: Path) -> list[FragmentRow]:
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        rows = []
        for r in reader:
            rows.append((
                r["accession"],
                int(r["start_pos"]),
                r["start_enzyme"],
                r["end_enzyme"],
                int(r["fragment_length"]),
            ))
    return rows

def merge_intervals(intervals: list[Interval]) -> list[Interval]:
    """Merge overlapping / touching half-open intervals [start, end)."""
    if not intervals:
        return []

    intervals = sorted(intervals)
    merged: list[Interval] = []
    cur_s, cur_e = intervals[0]

    for s, e in intervals[1:]:
        if s <= cur_e:
            cur_e = max(cur_e, e)
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = s, e

    merged.append((cur_s, cur_e))
    return merged


def merge_interval_dict(intervals_by_chrom: IntervalDict) -> IntervalDict:
    return {
        chrom: merge_intervals(intervals)
        for chrom, intervals in intervals_by_chrom.items()
        if intervals
    }


def interval_dict_length(intervals_by_chrom: IntervalDict) -> int:
    return sum(e - s for intervals in intervals_by_chrom.values() for s, e in intervals)


def load_bed_intervals(path: Path) -> IntervalDict:
    """Load BED (already merged or not) into chrom -> merged intervals."""
    by_chrom: IntervalDict = defaultdict(list)
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            c = line.rstrip().split("\t")
            if len(c) < 3:
                continue
            chrom = c[0]
            start = int(c[1])
            end = int(c[2])
            if end > start:
                by_chrom[chrom].append((start, end))
    return merge_interval_dict(by_chrom)


def subtract_sorted_intervals(a: list[Interval], b: list[Interval]) -> list[Interval]:
    """Return a - b for two sorted non-overlapping interval lists."""
    if not a:
        return []
    if not b:
        return a[:]

    out: list[Interval] = []
    j = 0
    n_b = len(b)

    for s, e in a:
        cur = s

        while j < n_b and b[j][1] <= cur:
            j += 1

        k = j
        while k < n_b and b[k][0] < e:
            bs, be = b[k]

            if bs > cur:
                out.append((cur, min(bs, e)))

            cur = max(cur, be)
            if cur >= e:
                break
            k += 1

        if cur < e:
            out.append((cur, e))

    return out


def subtract_interval_dict(a: IntervalDict, b: IntervalDict) -> IntervalDict:
    out: IntervalDict = {}
    for chrom, a_ints in a.items():
        b_ints = b.get(chrom, [])
        if not b_ints:
            out[chrom] = a_ints[:]
        else:
            out[chrom] = subtract_sorted_intervals(a_ints, b_ints)
    return {chrom: ints for chrom, ints in out.items() if ints}


def intersect_sorted_intervals(a: list[Interval], b: list[Interval]) -> int:
    """Return total overlap length of two sorted non-overlapping interval lists."""
    if not a or not b:
        return 0

    i = j = 0
    total = 0

    while i < len(a) and j < len(b):
        s1, e1 = a[i]
        s2, e2 = b[j]

        start = max(s1, s2)
        end = min(e1, e2)
        if start < end:
            total += end - start

        if e1 < e2:
            i += 1
        else:
            j += 1

    return total


def intersect_interval_dict(a: IntervalDict, b: IntervalDict) -> int:
    total = 0
    for chrom, a_ints in a.items():
        total += intersect_sorted_intervals(a_ints, b.get(chrom, []))
    return total


def union_interval_dict(a: IntervalDict, b: IntervalDict) -> IntervalDict:
    out: IntervalDict = defaultdict(list)
    for chrom, ints in a.items():
        out[chrom].extend(ints)
    for chrom, ints in b.items():
        out[chrom].extend(ints)
    return merge_interval_dict(out)


def read_rows_as_intervals(rows: list[FragmentRow]) -> IntervalDict:
    by_chrom: IntervalDict = defaultdict(list)
    for acc, start, _, _, fl in rows:
        by_chrom[acc].append((start, start + fl))
    return merge_interval_dict(by_chrom)


def write_filtered_bed(rows: list[FragmentRow], out: Path) -> int:
    total = 0
    with open(out, "w", encoding="utf-8") as fh:
        for acc, start, _, _, fl in sorted(rows):
            end = start + fl
            fh.write(f"{acc}\t{start}\t{end}\n")
            total += fl
    return total


def sort_merge(bedtools: str, inp: Path, out: Path):
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td) / "sorted.bed"

        with open(tmp, "w", encoding="utf-8") as out_sorted:
            import subprocess
            subprocess.run([bedtools, "sort", "-i", str(inp)],
                           stdout=out_sorted, check=True, text=True)

        with open(out, "w", encoding="utf-8") as out_final:
            import subprocess
            subprocess.run([bedtools, "merge", "-i", str(tmp)],
                           stdout=out_final, check=True, text=True)


def intersect_sum(bedtools: str, a: Path, b: Path) -> int:
    import subprocess
    proc = subprocess.run(
        [bedtools, "intersect", "-a", str(a), "-b", str(b), "-wo"],
        capture_output=True,
        text=True,
        check=True
    )

    total = 0
    for line in proc.stdout.splitlines():
        if line:
            total += int(line.split("\t")[-1])
    return total


def build_gff_sources(gff: Path, cache: Path) -> list[AnnotationSource]:
    skip = {"region", "chromosome"}
    by_label: dict[str, list[tuple[str, int, int]]] = defaultdict(list)

    with open_text(gff) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue

            c = line.split("\t")
            if len(c) < 5:
                continue

            feature = c[2]
            if feature in skip:
                continue

            chrom = c[0]
            start = int(c[3]) - 1
            end = int(c[4])

            if end <= start:
                continue

            by_label[f"gff:{feature}"].append((chrom, start, end))

    cache.mkdir(exist_ok=True, parents=True)

    sources = []
    for label in sorted(by_label):
        intervals = by_label[label]
        out = cache / f"{label.replace(':','_')}.bed"
        write_merged(intervals, out)
        sources.append(AnnotationSource(label=label, path=str(out)))

    return sources


def build_te_sources(te: Path, cache: Path) -> list[AnnotationSource]:
    by_label: dict[str, list[tuple[str, int, int]]] = defaultdict(list)

    with open_text(te) as fh:
        for i, line in enumerate(fh):
            if i < 3:
                continue

            p = line.split()
            if len(p) < 15:
                continue

            chrom = p[4]
            start = max(0, int(p[5]) - 1)
            end = int(p[6])

            cls = p[10].split("/")[0]
            fam = p[10].split("/")[1] if "/" in p[10] else p[10]

            label = f"te:{cls}_{fam}"
            by_label[label].append((chrom, start, end))

    cache.mkdir(exist_ok=True, parents=True)

    sources = []
    for label in sorted(by_label):
        intervals = by_label[label]
        out = cache / f"{label.replace(':','_')}.bed"
        write_merged(intervals, out)
        sources.append(AnnotationSource(label=label, path=str(out)))

    return sources


def write_merged(intervals, out_path: Path):
    if not intervals:
        return

    intervals.sort()
    merged = []

    c, s, e = intervals[0]
    for chrom, start, end in intervals[1:]:
        if chrom != c or start > e:
            merged.append((c, s, e))
            c, s, e = chrom, start, end
        else:
            e = max(e, end)

    merged.append((c, s, e))

    with open(out_path, "w", encoding="utf-8") as f:
        for chrom, start, end in merged:
            f.write(f"{chrom}\t{start}\t{end}\n")


def count_unique_category_bases(
    library_intervals: IntervalDict,
    sources: list[AnnotationSource],
    strand_multiplier: int = 2,
) -> dict[str, int]:
    """
    Count unique overlapping bases per category.

    Overlaps between categories in the same source list are resolved
    deterministically by source order so that a base is counted only once.
    """
    counts: dict[str, int] = {}
    claimed: IntervalDict = {}

    for source in sources:
        src = load_bed_intervals(Path(source.path))
        exclusive = subtract_interval_dict(src, claimed)
        overlap = intersect_interval_dict(library_intervals, exclusive)
        counts[source.label] = overlap * strand_multiplier
        claimed = union_interval_dict(claimed, exclusive)

    return counts


def process(combo, results_dir, sources, bedtools):
    combo_dir = Path(results_dir) / combo
    filtered = combo_dir / "filtered.csv"

    rows = read_filtered_csv(filtered)
    library_intervals = read_rows_as_intervals(rows)
    library_total = interval_dict_length(library_intervals)

    strand_multiplier = 2
    effective_total = library_total * strand_multiplier

    # Bed file is kept for compatibility/debugging.
    with tempfile.TemporaryDirectory(dir=combo_dir) as td:
        td = Path(td)
        bed = td / "filtered.bed"
        write_filtered_bed(rows, bed)

        # Count bases uniquely inside GFF and TE namespaces separately.
        # This means overlaps inside each namespace are never double-counted.
        gff_sources = [s for s in sources if s.label.startswith("gff:")]
        te_sources = [s for s in sources if s.label.startswith("te:")]

        counts = {}
        counts.update(count_unique_category_bases(library_intervals, gff_sources, strand_multiplier))
        counts.update(count_unique_category_bases(library_intervals, te_sources, strand_multiplier))

    out = combo_dir / "annotation_summary.csv"
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["category", "bases_overlapping", "pct_of_library"])

        w.writerow(["total_filtered_bases", effective_total, "100.00" if effective_total else "0.00"])

        for k, v in sorted(counts.items()):
            pct = (v / effective_total * 100) if effective_total else 0
            w.writerow([k, v, f"{pct:.2f}"])

    return combo, str(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--te")
    ap.add_argument("--annotation")
    ap.add_argument("--parallel", type=int, default=4)
    args = ap.parse_args()

    workdir = Path(args.workdir)
    results = workdir / "results"
    combos = load_combinations_csv(results / "combinations.csv")

    bedtools = check_bedtools()

    cache = results / "_annotation_cache"
    cache.mkdir(exist_ok=True)

    sources = []

    if args.annotation:
        sources += build_gff_sources(Path(args.annotation), cache)

    if args.te:
        sources += build_te_sources(Path(args.te), cache)

    print(f"[INFO] {len(sources)} annotation categories")

    from concurrent.futures import ProcessPoolExecutor, as_completed

    with ProcessPoolExecutor(max_workers=args.parallel) as ex:
        futures = [
            ex.submit(process, c, results, sources, bedtools)
            for c in combos
        ]

        for f in as_completed(futures):
            combo, out = f.result()
            print(combo, "->", out)


if __name__ == "__main__":
    main()