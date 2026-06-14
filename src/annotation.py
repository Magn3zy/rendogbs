#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class AnnotationSource:
    kind: str          # "gff" or "te"
    path: str          # direct GFF path or converted BED path
    label_prefix: str  # "gff:" or "te:"


FragmentRow = tuple[str, int, str, str, int]


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r", encoding="utf-8")


def check_bedtools() -> str:
    bedtools = shutil.which("bedtools")
    if bedtools is None:
        raise SystemExit(
            "[ERROR] bedtools not found in PATH.\n")
    return bedtools


def load_combinations_csv(path: Path) -> list[str]:
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise SystemExit(f"[ERROR] Missing header in {path}")
        if "enzyme_a" not in reader.fieldnames or "enzyme_b" not in reader.fieldnames:
            raise SystemExit(f"[ERROR] {path} must contain columns: enzyme_a,enzyme_b")
        return [f"{row['enzyme_a']}_{row['enzyme_b']}" for row in reader]


def read_filtered_csv(path: Path) -> list[FragmentRow]:
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise SystemExit(f"[ERROR] Missing header in {path}")
        required = {"accession", "start_pos", "start_enzyme", "end_enzyme", "fragment_length"}
        if not required.issubset(set(reader.fieldnames)):
            raise SystemExit(f"[ERROR] {path} must contain columns: {', '.join(sorted(required))}")

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


def write_filtered_bed(filtered_rows: list[FragmentRow], out_path: Path) -> int:
    total_bases = 0
    rows_sorted = sorted(filtered_rows, key=lambda x: (x[0], x[1], x[4]))
    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        for acc, start, _, _, fl in rows_sorted:
            end = start + fl
            fh.write(f"{acc}\t{start}\t{end}\t.\n")
            total_bases += fl
    return total_bases


def convert_repeatmasker_out_to_bed(te_path: Path, out_path: Path) -> None:
    with open_text(te_path) as fh, open(out_path, "w", encoding="utf-8", newline="") as out:
        for i, line in enumerate(fh):
            if i < 3:
                continue
            parts = line.split()
            if len(parts) < 15:
                continue

            try:
                chrom = parts[4]
                start = max(0, int(parts[5]) - 1)
                end = int(parts[6])

                te_class = parts[10].split("/")[0]
                te_family = parts[10].split("/")[1] if "/" in parts[10] else parts[10]
                label = f"te:{te_class}_{te_family}"

                out.write(f"{chrom}\t{start}\t{end}\t{label}\n")
            except (ValueError, IndexError):
                continue


def parse_gff_label_from_intersect_line(parts: list[str]) -> str:
    # A has 4 cols, GFF has 9 cols, overlap is last column.
    # GFF feature type is column 3 (0-based index 2) -> output index 4 + 2 = 6.
    return f"gff:{parts[6]}"


def parse_bed_label_from_intersect_line(parts: list[str]) -> str:
    # A has 4 cols, BED has 4 cols, overlap is last column.
    # Label is BED column 4 (0-based index 3) -> output index 4 + 3 = 7.
    return parts[7]


def run_bedtools_intersect(
    bedtools: str,
    filtered_bed: Path,
    annotation_path: Path,
    kind: str,
) -> dict[str, int]:
    cmd = [
        bedtools,
        "intersect",
        "-a",
        str(filtered_bed),
        "-b",
        str(annotation_path),
        "-wo",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"bedtools intersect failed for {annotation_path}:\n{proc.stderr.strip()}"
        )

    counts: dict[str, int] = defaultdict(int)
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2:
            continue
        overlap = int(parts[-1])
        if kind == "gff":
            label = parse_gff_label_from_intersect_line(parts)
        elif kind == "te":
            label = parse_bed_label_from_intersect_line(parts)
        else:
            raise ValueError(f"Unknown annotation kind: {kind}")
        counts[label] += overlap

    return dict(counts)


def write_annotation_summary(combo_dir: Path, total_bases: int, counts: dict[str, int]) -> Path:
    out_path = combo_dir / "annotation_summary.csv"
    with open(out_path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["category", "bases_overlapping", "pct_of_library"])
        w.writerow(["total_filtered_bases", total_bases, "100.00" if total_bases else "0.00"])
        for label, bases in sorted(counts.items()):
            pct = f"{(100 * bases / total_bases):.2f}" if total_bases else "0.00"
            w.writerow([label, bases, pct])
    return out_path


def process_combination(
    combo_name: str,
    results_dir: str,
    sources: list[AnnotationSource],
    bedtools: str,
) -> dict[str, Any]:
    results_path = Path(results_dir)
    combo_dir = results_path / combo_name
    filtered_csv = combo_dir / "filtered.csv"
    if not filtered_csv.exists():
        raise FileNotFoundError(f"Missing filtered.csv: {filtered_csv}")

    filtered_rows = read_filtered_csv(filtered_csv)

    with tempfile.TemporaryDirectory(dir=combo_dir) as tmpdir:
        tmpdir_path = Path(tmpdir)
        filtered_bed = tmpdir_path / "filtered.bed"
        total_bases = write_filtered_bed(filtered_rows, filtered_bed)

        counts: dict[str, int] = defaultdict(int)
        for source in sources:
            source_counts = run_bedtools_intersect(
                bedtools=bedtools,
                filtered_bed=filtered_bed,
                annotation_path=Path(source.path),
                kind=source.kind,
            )
            for label, bases in source_counts.items():
                counts[label] += bases

    summary_path = write_annotation_summary(combo_dir, total_bases, counts)

    return {
        "combo": combo_name,
        "filtered": len(filtered_rows),
        "bases": total_bases,
        "summary": str(summary_path),
    }


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="annotation.py",
        description=(
            "Annotate existing filtered.csv files in <workdir>/results/<combo>/ "
            "using bedtools intersect."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--workdir", required=True, help="Root workdir containing results/")
    p.add_argument("--te", default=None, help="RepeatMasker .out file")
    p.add_argument("--annotation", default=None, help="GFF3/GFF/GTF annotation file")
    p.add_argument("--threads", type=int, default=4, help="Number of combinations processed in parallel")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    if args.threads < 1:
        raise SystemExit("[ERROR] --threads must be >= 1")

    workdir = Path(args.workdir)
    results_dir = workdir / "results"
    combinations_csv = results_dir / "combinations.csv"

    if not results_dir.exists():
        raise SystemExit(f"[ERROR] results directory not found: {results_dir}")
    if not combinations_csv.exists():
        raise SystemExit(f"[ERROR] combinations.csv not found: {combinations_csv}")
    if not args.te and not args.annotation:
        raise SystemExit("[ERROR] Provide --te, --annotation, or both")

    bedtools = check_bedtools()
    combo_names = load_combinations_csv(combinations_csv)

    for combo_name in combo_names:
        combo_dir = results_dir / combo_name
        if not combo_dir.is_dir():
            raise SystemExit(f"[ERROR] Missing combination directory: {combo_dir}")

    sources: list[AnnotationSource] = []

    if args.te:
        te_path = Path(args.te)
        if not te_path.exists():
            raise SystemExit(f"[ERROR] RepeatMasker file not found: {te_path}")
        converted_bed = results_dir / "repeatmasker_converted.bed"
        print(f"[INFO] Converting {te_path} -> {converted_bed}")
        convert_repeatmasker_out_to_bed(te_path, converted_bed)
        sources.append(AnnotationSource(kind="te", path=str(converted_bed), label_prefix="te:"))

    if args.annotation:
        ann_path = Path(args.annotation)
        if not ann_path.exists():
            raise SystemExit(f"[ERROR] Annotation file not found: {ann_path}")
        sources.append(AnnotationSource(kind="gff", path=str(ann_path), label_prefix="gff:"))

    print(f"[INFO] Processing {len(combo_names)} combination(s) | threads={args.threads}")

    results: list[dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = {
            executor.submit(process_combination, combo, str(results_dir), sources, bedtools): combo
            for combo in combo_names
        }
        for idx, future in enumerate(as_completed(futures), start=1):
            combo = futures[future]
            try:
                r = future.result()
                results.append(r)
                print(f"  [{idx:>3}/{len(futures)}] {combo} -> {r['summary']}")
            except Exception as exc:
                print(f"  [ERROR] {combo}: {exc}", file=sys.stderr)
                raise

    print(f"[INFO] Done. Processed={len(results)}")


if __name__ == "__main__":
    main()
