#!/usr/bin/env python3
# --workdir 
import csv
import argparse
from pathlib import Path


def read_csv(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def load_stats(path: Path) -> dict:
    rows = read_csv(path)
    return rows[0] if rows else {}


def load_gc(path: Path) -> dict:
    return {r["metric"]: r["value"] for r in read_csv(path) if r.get("metric")}


def load_distribution(path: Path) -> dict:
    return {r["length_range"]: r["count"] for r in read_csv(path) if r.get("length_range")}


def load_annotation(path: Path) -> dict:
    out = {}
    for r in read_csv(path):
        cat = r.get("category", "")
        if cat:
            out[cat] = r.get("pct_of_library", "")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workdir", required=True)
    args = ap.parse_args()

    results_dir = Path(args.workdir) / "results"
    combos = [
        f"{r['enzyme_a']}_{r['enzyme_b']}"
        for r in read_csv(results_dir / "combinations.csv")
        if r.get("enzyme_a") and r.get("enzyme_b")
    ]

    rows = []
    all_cols = ["combination"]

    for combo in combos:
        combo_dir = results_dir / combo
        row = {"combination": combo}
        row.update(load_stats(combo_dir / "statistics_cutting.csv"))
        row.update(load_gc(combo_dir / "gc_metrics.csv"))
        row.update(load_distribution(combo_dir / "distribution.csv"))
        row.update(load_annotation(combo_dir / "annotation.csv"))
        rows.append(row)
        for k in row:
            if k not in all_cols:
                all_cols.append(k)

    out_path = results_dir / "summary.tsv"
    with out_path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(all_cols)
        for row in rows:
            w.writerow([row.get(col, "") for col in all_cols])

    print(f"[INFO] written -> {out_path}")


if __name__ == "__main__":
    main()