#!/usr/bin/env python3
# annotation.py
# Copyright (c) 2026 Eliška Korbová ORCID 0009-0004-1247-0808
#
# Annotation of ddRAD fragments using bedtools
#
# This workflow uses Bedtools.
# Bedtools: https://github.com/arq5x/bedtools2
# Copyright (c) Aaron Quinlan
# Licensed under the MIT License.
#
# Output: annotation_summary.csv

from __future__ import annotations
import argparse
import csv
import gzip
import shutil
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class AnnotationSource:
    label: str
    path: Path


FragmentRow = tuple[str, int, str, str, int]
BedRecord = tuple[str, int, int]

_WORKER_SOURCES: list[AnnotationSource] | None = None
_WORKER_BEDTOOLS: str | None = None

def _init_annotation_worker(sources: list[AnnotationSource], bedtools: str) -> None:
    global _WORKER_SOURCES, _WORKER_BEDTOOLS
    _WORKER_SOURCES = sources
    _WORKER_BEDTOOLS = bedtools


GFF_PRIORITY = [
    "CDS",
    "five_prime_UTR",
    "three_prime_UTR",
    "exon",
    "mRNA",
    "gene",
]


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r", encoding="utf-8")


def check_bedtools() -> str:
    bedtools = shutil.which("bedtools")
    if bedtools is None:
        raise SystemExit("[ERROR] bedtools not found in PATH")
    return bedtools


def file_is_empty(path: Path) -> bool:
    return (not path.exists()) or path.stat().st_size == 0


def load_combinations_csv(path: Path) -> list[str]:
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        return [f"{r['enzyme_a']}_{r['enzyme_b']}" for r in reader]


def read_filtered_csv(path: Path) -> list[FragmentRow]:
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        rows: list[FragmentRow] = []
        for r in reader:
            rows.append(
                (
                    r["accession"],
                    int(r["start_pos"]),
                    r["start_enzyme"],
                    r["end_enzyme"],
                    int(r["fragment_length"]),
                )
            )
    return rows


def write_bed_records(records: Iterable[BedRecord], out: Path) -> None:
    with open(out, "w", encoding="utf-8") as fh:
        for chrom, start, end in records:
            if end > start:
                fh.write(f"{chrom}\t{start}\t{end}\n")


def bed_length(path: Path) -> int:
    total = 0
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if not line.strip():
                continue
            c = line.rstrip().split("\t")
            if len(c) < 3:
                continue
            s = int(c[1])
            e = int(c[2])
            if e > s:
                total += e - s
    return total


def run_bedtools_to_file(bedtools: str, args: list[str], out: Path) -> None:
    with open(out, "w", encoding="utf-8") as fh:
        subprocess.run([bedtools, *args], stdout=fh, text=True, check=True)

# Wait for all processes and raise if any exited non-zero
def _check_procs(procs: list[subprocess.Popen], labels: list[str]) -> None:
    for proc, label in zip(procs, labels):
        rc = proc.wait()
        if rc != 0:
            raise SystemExit(f"[ERROR] {label} failed with exit code {rc}")

# Sort and merge a BED file - without temp files
def bedtools_sort_merge(bedtools: str, inp: Path, out: Path) -> None:
    if file_is_empty(inp):
        out.write_text("", encoding="utf-8")
        return

    # bedtools sort -i <inp> | bedtools merge -i stdin
    p_sort = subprocess.Popen(
        [bedtools, "sort", "-i", str(inp)],
        stdout=subprocess.PIPE,
        text=True,
    )
    with open(out, "w", encoding="utf-8") as fh:
        p_merge = subprocess.Popen(
            [bedtools, "merge", "-i", "stdin"],
            stdin=p_sort.stdout,
            stdout=fh,
            text=True,
        )
    # Close our copy of the write-end so p_merge receives EOF when p_sort exits.
    if p_sort.stdout:
        p_sort.stdout.close()
    _check_procs([p_sort, p_merge], ["bedtools sort", "bedtools merge"])


def bedtools_subtract(bedtools: str, a: Path, b: Path, out: Path) -> None:
    if file_is_empty(a):
        out.write_text("", encoding="utf-8")
        return
    if file_is_empty(b):
        shutil.copyfile(a, out)
        return
    run_bedtools_to_file(bedtools, ["subtract", "-sorted", "-a", str(a), "-b", str(b)], out)

# Intersect a and b, then sort and merge without temp files
def bedtools_intersect_to_file(bedtools: str, a: Path, b: Path, out: Path) -> None:
    if file_is_empty(a) or file_is_empty(b):
        out.write_text("", encoding="utf-8")
        return

    # bedtools intersect -sorted -a <a> -b <b> | bedtools sort | bedtools merge
    p_isect = subprocess.Popen(
        [bedtools, "intersect", "-sorted", "-a", str(a), "-b", str(b)],
        stdout=subprocess.PIPE,
        text=True,
    )
    p_sort = subprocess.Popen(
        [bedtools, "sort", "-i", "stdin"],
        stdin=p_isect.stdout,
        stdout=subprocess.PIPE,
        text=True,
    )
    with open(out, "w", encoding="utf-8") as fh:
        p_merge = subprocess.Popen(
            [bedtools, "merge", "-i", "stdin"],
            stdin=p_sort.stdout,
            stdout=fh,
            text=True,
        )
    if p_isect.stdout:
        p_isect.stdout.close()
    if p_sort.stdout:
        p_sort.stdout.close()
    _check_procs(
        [p_isect, p_sort, p_merge],
        ["bedtools intersect", "bedtools sort", "bedtools merge"],
    )


def read_rows_as_intervals(rows: list[FragmentRow]) -> list[BedRecord]:
    return [(acc, start, start + fl) for acc, start, _, _, fl in rows if fl > 0]


def build_library_bed(rows: list[FragmentRow], workdir: Path, bedtools: str) -> tuple[Path, int]:
    raw = workdir / "filtered.bed"
    merged = workdir / "filtered_merged.bed"

    records = read_rows_as_intervals(rows)
    write_bed_records(records, raw)
    bedtools_sort_merge(bedtools, raw, merged)

    return merged, bed_length(merged)


def build_gff_sources(gff: Path, cache: Path, bedtools: str) -> list[AnnotationSource]:
    skip = {"region", "chromosome"}
    by_label: dict[str, list[BedRecord]] = defaultdict(list)

    with open_text(gff) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue

            c = line.rstrip("\n").split("\t")
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

    ordered_labels: list[str] = []
    for feat in GFF_PRIORITY:
        label = f"gff:{feat}"
        if label in by_label:
            ordered_labels.append(label)

    ordered_labels.extend(sorted(label for label in by_label if label not in ordered_labels))

    sources: list[AnnotationSource] = []
    for label in ordered_labels:
        safe = label.replace(":", "_")
        raw = cache / f"{safe}.raw.bed"
        merged = cache / f"{safe}.bed"

        write_bed_records(by_label[label], raw)
        bedtools_sort_merge(bedtools, raw, merged)
        sources.append(AnnotationSource(label=label, path=merged))

    return sources


def build_te_sources(te: Path, cache: Path, bedtools: str) -> list[AnnotationSource]:
    by_label: dict[str, list[BedRecord]] = defaultdict(list)

    with open_text(te) as fh:
        for i, line in enumerate(fh):
            if i < 3:
                continue
            if not line.strip() or line.startswith("#"):
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
            if end > start:
                by_label[label].append((chrom, start, end))

    cache.mkdir(exist_ok=True, parents=True)

    sources: list[AnnotationSource] = []
    for label in sorted(by_label):
        safe = label.replace(":", "_")
        raw = cache / f"{safe}.raw.bed"
        merged = cache / f"{safe}.bed"

        write_bed_records(by_label[label], raw)
        bedtools_sort_merge(bedtools, raw, merged)
        sources.append(AnnotationSource(label=label, path=merged))

    return sources


def combine_beds(bedtools: str, left: Path, right: Path, out: Path) -> None:
    if file_is_empty(left):
        shutil.copyfile(right, out)
        return
    if file_is_empty(right):
        shutil.copyfile(left, out)
        return

    with tempfile.TemporaryDirectory() as td_str:
        tmp = Path(td_str) / "combined.bed"
        with open(tmp, "w", encoding="utf-8") as fh:
            with open(left, encoding="utf-8") as a:
                shutil.copyfileobj(a, fh)
            with open(right, encoding="utf-8") as b:
                shutil.copyfileobj(b, fh)
        bedtools_sort_merge(bedtools, tmp, out)


def count_unique_category_bases(
    bedtools: str,
    library_bed: Path,
    sources: list[AnnotationSource],
    workdir: Path,
) -> dict[str, int]:
    counts: dict[str, int] = {}

    # Pre-intersect every source with the library work on library-sized BED files
    library_isect: list[Path] = []
    for idx, source in enumerate(sources, start=1):
        isect = workdir / f"{idx:03d}_{source.label.replace(':', '_')}_in_library.bed"
        bedtools_intersect_to_file(bedtools, library_bed, source.path, isect)
        library_isect.append(isect)

    # Walk through categories in priority order, subtracting already-claimed bases
    claimed = workdir / "claimed.bed"
    claimed.write_text("", encoding="utf-8")

    for idx, (source, isect) in enumerate(zip(sources, library_isect), start=1):
        exclusive = workdir / f"{idx:03d}_{source.label.replace(':', '_')}_exclusive.bed"
        bedtools_subtract(bedtools, isect, claimed, exclusive)

        counts[source.label] = bed_length(exclusive)

        if file_is_empty(exclusive):
            continue

        # Update claimed by merging the new exclusive region in, small files intersect runs
        # once per category
        new_claimed = workdir / f"{idx:03d}_claimed.bed"
        combine_beds(bedtools, claimed, exclusive, new_claimed)
        claimed = new_claimed

    return counts


def process(combo: str, results_dir: Path):
    assert _WORKER_SOURCES is not None
    assert _WORKER_BEDTOOLS is not None
    combo_dir = results_dir / combo
    filtered = combo_dir / "filtered.csv"

    if not filtered.exists():
        raise SystemExit(f"[ERROR] Missing file: {filtered}")

    rows = read_filtered_csv(filtered)

    with tempfile.TemporaryDirectory() as td_str:
        td = Path(td_str)
        library_bed, library_total = build_library_bed(rows, td, _WORKER_BEDTOOLS)
        counts = count_unique_category_bases(_WORKER_BEDTOOLS, library_bed, _WORKER_SOURCES, td)

    out = combo_dir / "annotation_summary.csv"
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["category", "bases_overlapping", "pct_of_library"])
        w.writerow(["total_filtered_bases", library_total, "100.00" if library_total else "0.00"])
        for label, value in counts.items():
            pct = (value / library_total * 100) if library_total else 0.0
            w.writerow([label, value, f"{pct:.2f}"])

    return combo, str(out)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    ap.add_argument("--te")
    ap.add_argument("--annotation")
    ap.add_argument("--parallel", type=int, default=4)
    args = ap.parse_args()

    if not args.annotation and not args.te:
        raise SystemExit("[ERROR] At least one of --annotation or --te must be provided")

    workdir = Path(args.workdir)
    results = workdir / "results"
    combos = load_combinations_csv(results / "combinations.csv")

    bedtools = check_bedtools()

    cache = results / "_annotation_cache"
    cache.mkdir(exist_ok=True)

    sources: list[AnnotationSource] = []

    # GFF first, then TE.
    if args.annotation:
        sources.extend(build_gff_sources(Path(args.annotation), cache, bedtools))

    if args.te:
        sources.extend(build_te_sources(Path(args.te), cache, bedtools))

    print(f"[INFO] {len(sources)} annotation categories")

    with ProcessPoolExecutor(
        max_workers=args.parallel,
        initializer=_init_annotation_worker,
        initargs=(sources, bedtools),
    ) as ex:
        futures = [
            ex.submit(process, combo, results)
            for combo in combos
        ]
        for f in as_completed(futures):
            combo, out = f.result()
            print(combo, "->", out)

if __name__ == "__main__":
    main()