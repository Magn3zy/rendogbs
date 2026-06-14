#!/usr/bin/env python3
"""
rendogbs_plots.py  -  Post-processing: plots only

Reads  : <workdir>/results/
Writes : <workdir>/results/plots/

Generates:
  - fragment length heatmap (10 bp bins, log scale)
  - chromosome distribution heatmap + bar plot
  - GC content summary plot
  - annotation coverage stacked plot
  - size distribution (100 bp bins + window highlight)
Opravit graf na anotaci a TE rozdelit!!!
"""

import os
import csv
import argparse
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import LogNorm
import matplotlib.patches as mpatches


# global plot style
PALETTE = "viridis"
FIG_EXT = "png"
_DPI = 300

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
})


def parse_size_range(s):
    parts = s.split("-")
    return int(parts[0]), int(parts[1])


def discover_combos(results_dir):
    """Find all valid enzyme-combination directories."""
    combos = []
    for entry in sorted(os.scandir(results_dir), key=lambda e: e.name):
        if entry.is_dir() and os.path.isfile(
            os.path.join(entry.path, "fragments.csv")
        ):
            combos.append(entry.name)
    return combos


def read_csv_rows(path):
    """Safe CSV loader."""
    if not os.path.isfile(path):
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def read_distribution_csv(path):
    """length_range -> count"""
    d = {}
    for row in read_csv_rows(path):
        try:
            d[row["length_range"]] = int(row["count"])
        except (KeyError, ValueError):
            pass
    return d


def read_gc_metrics_csv(path):
    """metric -> value"""
    return {
        row["metric"]: row["value"]
        for row in read_csv_rows(path)
        if "metric" in row
    }


def read_annotation_summary_csv(path):
    """category -> percentage"""
    d = {}
    for row in read_csv_rows(path):
        cat = row.get("category", "")
        if cat == "total_filtered_bases":
            continue
        try:
            d[cat] = float(row["percentage"])
        except (KeyError, ValueError):
            pass
    return d


def read_contig_lengths_txt(results_dir):
    path = os.path.join(results_dir, "contig_lengths.txt")
    rows = []
    with open(path) as fh:
        next(fh)
        for line in fh:
            parts = line.rstrip().split("\t")
            if len(parts) == 2:
                rows.append((parts[0], int(parts[1])))
    return rows


def save(fig, plots_dir, name):
    os.makedirs(plots_dir, exist_ok=True)
    path = os.path.join(plots_dir, f"{name}.{FIG_EXT}")
    fig.savefig(path, dpi=_DPI, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] {path}")
    return path


def plot_heatmap(results_dir, combos, plots_dir):
    """Fragment length heatmap (10bp bins, log scale)."""
    bin_width = 10
    max_len = 1000
    edges = list(range(0, max_len + 1, bin_width)) + [float("inf")]

    mat = np.zeros((len(combos), len(edges) - 1), dtype=np.int64)

    for i, combo in enumerate(combos):
        rows = read_csv_rows(os.path.join(results_dir, combo, "fragments.csv"))
        lengths = np.array([int(r["fragment_length"]) for r in rows] or [0])
        clipped = np.minimum(lengths, max_len)
        mat[i], _ = np.histogram(clipped, bins=edges)

    if mat.max() == 0:
        print("[WARN] heatmap empty")
        return

    fig, ax = plt.subplots(figsize=(12, max(4, len(combos) * 0.4)))

    im = ax.imshow(
        mat,
        aspect="auto",
        cmap=PALETTE,
        norm=LogNorm(vmin=1, vmax=max(mat.max(), 1)),
        origin="lower",
    )

    ax.set_title("Fragment length heatmap")
    ax.set_xlabel("Length (bp)")
    ax.set_ylabel("Combination")

    ax.set_yticks(np.arange(len(combos)))
    ax.set_yticklabels(combos, fontsize=8)

    save(fig, plots_dir, "heatmap_fragment_lengths")


def plot_chrom_distribution(results_dir, combos, plots_dir, n_chroms):
    """Chromosome distribution heatmap + barplot."""
    contigs = read_contig_lengths_txt(results_dir)
    chroms = [c for c, _ in contigs[:n_chroms]]

    idx = {c: i for i, c in enumerate(chroms)}
    mat = np.zeros((len(combos), len(chroms)), dtype=np.int64)

    for i, combo in enumerate(combos):
        rows = read_csv_rows(os.path.join(results_dir, combo, "filtered.csv"))
        for r in rows:
            if r["accession"] in idx:
                mat[i, idx[r["accession"]]] += 1

    if mat.max() == 0:
        print("[WARN] chrom distribution empty")
        return

    # heatmap
    fig, ax = plt.subplots(figsize=(10, 4))
    im = ax.imshow(mat, aspect="auto", cmap="YlOrRd", origin="lower")
    ax.set_title("Chromosome distribution")
    ax.set_xticks(range(len(chroms)))
    ax.set_xticklabels(chroms, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(combos)))
    ax.set_yticklabels(combos, fontsize=7)

    save(fig, plots_dir, "heatmap_chrom_distribution")

    # barplot
    fig, ax = plt.subplots(figsize=(10, 4))
    x = np.arange(len(chroms))
    width = 0.8 / max(len(combos), 1)

    for j, combo in enumerate(combos):
        ax.bar(x + j * width, mat[j], width=width, label=combo)

    ax.set_title("Chromosome counts (bar)")
    ax.set_xticks(x)
    ax.set_xticklabels(chroms, rotation=45, ha="right")

    save(fig, plots_dir, "bar_chrom_distribution")


def plot_gc_distribution(results_dir, combos, plots_dir):
    labels, means, stds = [], [], []

    for combo in combos:
        gc = read_gc_metrics_csv(os.path.join(results_dir, combo, "gc_metrics.csv"))
        if not gc:
            continue
        try:
            labels.append(combo)
            means.append(float(gc["gc_mean_pct"]))
            stds.append(float(gc["gc_std_pct"]))
        except:
            continue

    if not labels:
        print("[WARN] GC empty")
        return

    fig, ax = plt.subplots(figsize=(8, 4))
    x = np.arange(len(labels))

    ax.bar(x, means)
    ax.scatter(x, means)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_title("GC content")

    save(fig, plots_dir, "gc_distribution")


def plot_annotation_coverage(results_dir, combos, plots_dir):
    cats = []
    data = {}

    for c in combos:
        ann = read_annotation_summary_csv(
            os.path.join(results_dir, c, "annotation_summary.csv")
        )
        if ann:
            data[c] = ann
            cats += list(ann.keys())

    cats = sorted(set(cats))
    if not data:
        print("[WARN] annotation empty")
        return

    fig, ax = plt.subplots(figsize=(6, max(3, len(data) * 0.4)))

    for i, combo in enumerate(data):
        left = 0
        for cat in cats:
            val = data[combo].get(cat, 0)
            ax.barh(i, val, left=left)
            left += val

    ax.set_yticks(range(len(data)))
    ax.set_yticklabels(list(data.keys()))
    ax.set_title("Annotation coverage")

    save(fig, plots_dir, "annotation_coverage")


def plot_size_distributions(results_dir, combos, plots_dir, low, high):
    std = [f"{b}-{b+99}" for b in range(0, 1000, 100)] + [">=1000"]
    x = np.arange(len(std))

    fig, ax = plt.subplots(figsize=(10, 4))

    for c in combos:
        dist = read_distribution_csv(
            os.path.join(results_dir, c, "distribution.csv")
        )
        y = [dist.get(s, 0) for s in std]
        ax.plot(x, y, label=c)

    ax.axvspan(low // 100, high // 100, alpha=0.2)
    ax.set_title("Size distribution")

    save(fig, plots_dir, "size_distributions")


def parse_args():
    p = argparse.ArgumentParser()

    p.add_argument("--workdir", required=True)
    p.add_argument("--size", required=True)
    p.add_argument("--chroms", type=int, default=10)
    p.add_argument("--dpi", type=int, default=300)

    return p.parse_args()


def main():
    args = parse_args()

    global _DPI
    _DPI = args.dpi

    low, high = parse_size_range(args.size)

    results_dir = os.path.join(args.workdir, "results")
    plots_dir = os.path.join(results_dir, "plots")

    combos = discover_combos(results_dir)

    print(f"combos: {len(combos)}")

    plot_heatmap(results_dir, combos, plots_dir)
    plot_chrom_distribution(results_dir, combos, plots_dir, args.chroms)
    plot_gc_distribution(results_dir, combos, plots_dir)
    plot_annotation_coverage(results_dir, combos, plots_dir)
    plot_size_distributions(results_dir, combos, plots_dir, low, high)


if __name__ == "__main__":
    main()