#!/usr/bin/env python3
import os
import csv
import argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import LogNorm


# global plot style
PALETTE = "viridis"
FIG_EXT = "png"
_DPI = 300

# categories that represent the whole genome (100%) – must be excluded
_SKIP_CATS = {"total_filtered_bases", "gff:region"}

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
    """Safe CSV loader (comma delimiter)."""
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
    """
    category -> percentage (pct_of_library).
    Skips total_filtered_bases and gff:region (both represent 100% of genome
    and would inflate the stacked bar beyond 100%).
    """
    d = {}
    for row in read_csv_rows(path):
        cat = row.get("category", "")
        if cat in _SKIP_CATS:
            continue
        try:
            d[cat] = float(row["pct_of_library"])
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
    ax.imshow(
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
    ax.imshow(mat, aspect="auto", cmap="YlOrRd", origin="lower")
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
    labels  = []
    means   = []
    medians = []
    stds    = []
    mins_   = []
    maxs_   = []

    for combo in combos:
        gc = read_gc_metrics_csv(os.path.join(results_dir, combo, "gc_metrics.csv"))
        if not gc or gc.get("gc_mean_pct") == "n/a":
            continue
        try:
            labels.append(combo)
            means.append(float(gc["gc_mean_pct"]))
            medians.append(float(gc["gc_median_pct"]))
            stds.append(float(gc["gc_std_pct"]))
            mins_.append(float(gc["gc_min_pct"]))
            maxs_.append(float(gc["gc_max_pct"]))
        except (KeyError, ValueError):
            continue

    if not labels:
        print("  [WARN] GC plot: no gc_metrics.csv data found, skipping.")
        return

    n = len(labels)
    fig_w = max(8, n * 0.45 + 2)
    fig, ax = plt.subplots(figsize=(fig_w, 5))
    x = np.arange(n)

    ax.bar(
        x,
        [s * 2 for s in stds],
        bottom=[m - s for m, s in zip(means, stds)],
        width=0.55,
        color="#4393c3",
        alpha=0.6,
        label="mean ± 1 SD",
    )
    ax.scatter(x, means, color="#2166ac", zorder=5, s=40, label="mean")
    ax.scatter(x, medians, color="#d6604d", zorder=5, s=40, marker="D", label="median")

    for i in range(n):
        ax.plot([x[i], x[i]], [mins_[i], maxs_[i]], color="grey", linewidth=1.0, zorder=3)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=55, ha="right", fontsize=max(6, min(9, 120 // n)))
    ax.set_ylabel("GC content (%)")
    ax.set_title("GC content of filtered fragments — mean ± SD (whiskers = min/max)")
    ax.legend(loc="upper right", fontsize=8)
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{v:.0f}%"))
    save(fig, plots_dir, "gc_distribution")


def plot_annotation_coverage(results_dir, combos, plots_dir):
    """
    Two side-by-side horizontal stacked bar plots:
      left  – gff:* categories
      right – te:*  categories

    Each bar fills exactly to 100% – remaining space shown as 'unannotated'.
    gff:region and total_filtered_bases are excluded (they represent 100% genome).

    If no annotation_summary.csv files are found, the plot is skipped entirely.
    If only one prefix (gff or te) has data, a single plot is produced.
    """
    gff_cats_set, te_cats_set = set(), set()
    data = {}

    for c in combos:
        ann = read_annotation_summary_csv(
            os.path.join(results_dir, c, "annotation_summary.csv")
        )
        if not ann:
            continue
        data[c] = ann
        for cat in ann:
            if cat.startswith("gff:"):
                gff_cats_set.add(cat)
            elif cat.startswith("te:"):
                te_cats_set.add(cat)

    if not data:
        print("[SKIP] annotation_coverage – no annotation_summary.csv files found")
        return

    gff_cats = sorted(gff_cats_set)
    te_cats  = sorted(te_cats_set)
    combos_found = list(data.keys())
    n = len(combos_found)

    has_gff = bool(gff_cats)
    has_te  = bool(te_cats)
    n_panels = has_gff + has_te

    gff_colors = list(plt.cm.tab20.colors)
    te_colors  = list(plt.cm.Set2.colors)

    fig, axes = plt.subplots(
        1, n_panels,
        figsize=(max(6, n * 0.6) * n_panels, max(4, n * 0.5)),
        sharey=True,
        squeeze=False,
    )
    axes = axes[0]

    def draw_stacked(ax, cats, colors, title):
        # draw each category as a stacked segment
        for j, cat in enumerate(cats):
            lefts = [
                sum(data[c].get(cats[k], 0) for k in range(j))
                for c in combos_found
            ]
            vals = [data[c].get(cat, 0) for c in combos_found]
            ax.barh(
                range(n), vals, left=lefts,
                color=colors[j % len(colors)],
                label=cat.split(":", 1)[1],
            )

        # fill remainder to exactly 100% as "unannotated"
        for i, c in enumerate(combos_found):
            total = sum(data[c].get(cat, 0) for cat in cats)
            remainder = max(0.0, 100.0 - total)
            if remainder > 0:
                ax.barh(i, remainder, left=total, color="#cccccc",
                        label="_nolegend_")

        # single grey entry in legend
        ax.barh([], [], color="#cccccc", label="unannotated")

        ax.set_yticks(range(n))
        ax.set_yticklabels(combos_found, fontsize=8)
        ax.set_xlabel("% of library")
        ax.set_xlim(0, 100)
        ax.set_title(title, fontsize=11)
        ax.legend(loc="lower right", fontsize=7, framealpha=0.7, ncol=1)
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))

    panel_idx = 0
    if has_gff:
        draw_stacked(axes[panel_idx], gff_cats, gff_colors, "GFF annotation coverage")
        panel_idx += 1
    if has_te:
        draw_stacked(axes[panel_idx], te_cats, te_colors, "TE annotation coverage")

    fig.suptitle("Annotation coverage", fontsize=13, y=1.01)
    fig.tight_layout()
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
    p.add_argument("--workdir", required=True,
                   help="Working directory containing results/ subfolder")
    p.add_argument("--size", required=True,
                   help="Size window to highlight, e.g. 150-600")
    p.add_argument("--chroms", type=int, default=10,
                   help="Number of top chromosomes to plot (default: 10)")
    p.add_argument("--dpi", type=int, default=300,
                   help="Output DPI (default: 300)")
    return p.parse_args()


def main():
    args = parse_args()

    global _DPI
    _DPI = args.dpi

    low, high = parse_size_range(args.size)

    results_dir = os.path.join(args.workdir, "results")
    plots_dir   = os.path.join(results_dir, "plots")

    combos = discover_combos(results_dir)
    print(f"combos found: {len(combos)}")

    plot_heatmap(results_dir, combos, plots_dir)
    plot_chrom_distribution(results_dir, combos, plots_dir, args.chroms)
    plot_gc_distribution(results_dir, combos, plots_dir)
    plot_annotation_coverage(results_dir, combos, plots_dir)
    plot_size_distributions(results_dir, combos, plots_dir, low, high)


if __name__ == "__main__":
    main()