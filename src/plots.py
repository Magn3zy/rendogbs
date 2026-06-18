#!/usr/bin/env python3
# plots.py
# Copyright (c) 2026 Eliška Korbová ORCID 0009-0004-1247-0808
#
# Plot generation as pipeline summary for easy visualization
# Output: annotation_coverage.png, bar_chom_distribution.png, gc_distribution.png,
#         heatmap_fragment_lengths.png, heatmap_chrom_distribution.png, size_distribution.png

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
_DPI    = 300

# categories that represent the whole genome (100%) – excluded from annotation plots
_SKIP_CATS = {"total_filtered_bases", "gff:region"}

# no annotation info - fix color
UNANNOTATED_COLOR = "#b0b0b0"

plt.rcParams.update({
    "font.family":     "DejaVu Sans",
    "axes.titlesize":  13,
    "axes.labelsize":  11,
    "xtick.labelsize":  9,
    "ytick.labelsize":  9,
})

def parse_size_range(s):
    parts = s.split("-")
    return int(parts[0]), int(parts[1])


def discover_combos(results_dir):
    combos = []
    for entry in sorted(os.scandir(results_dir), key=lambda e: e.name):
        if entry.is_dir() and os.path.isfile(
                os.path.join(entry.path, "fragments.csv")):
            combos.append(entry.name)
    return combos

# return path to a file inside a combination results directory
def combo_file(results_dir, combo, filename):
    return os.path.join(results_dir, combo, filename)


def read_csv_rows(path):
    if not os.path.isfile(path):
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def read_distribution_csv(path):
    d = {}
    for row in read_csv_rows(path):
        try:
            d[row["length_range"]] = int(row["count"])
        except (KeyError, ValueError):
            pass
    return d


def read_gc_metrics_csv(path):
    return {row["metric"]: row["value"] for row in read_csv_rows(path)
            if "metric" in row}


def read_annotation_summary_csv(path):
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
    print(f"  [plot] saved -> {path}")
    return path

# heatmap fragment lengths fragments.csv
def plot_heatmap(results_dir, combos, plots_dir):
    bin_width = 10
    max_len   = 1000
    edges     = list(range(0, max_len + 1, bin_width)) + [float("inf")]
    n_bins    = len(edges) - 1

    mat = np.zeros((len(combos), n_bins), dtype=np.int64)

    for i, combo in enumerate(combos):
        rows = read_csv_rows(combo_file(results_dir, combo, "fragments.csv"))
        lengths = np.array([int(r["fragment_length"]) for r in rows] or [0],
                           dtype=np.int64)
        clipped = np.minimum(lengths, max_len)
        mat[i], _ = np.histogram(clipped, bins=edges)

    if mat.max() == 0:
        print("  [WARN] heatmap: no fragment data found, skipping.")
        return

    fig_h = max(4, len(combos) * 0.40 + 1.5)
    fig_w = max(12, n_bins * 0.12 + 3)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    im = ax.imshow(
        mat, aspect="auto", cmap=PALETTE,
        norm=LogNorm(vmin=1, vmax=max(int(mat.max()), 1)),
        origin="lower",
    )

    ax.set_yticks(np.arange(len(combos)))
    ax.set_yticklabels(combos,
                       fontsize=max(6, min(9, 120 // max(len(combos), 1))))

    xt_idxs   = np.arange(0, max_len + 1, 100) // bin_width
    xt_labels = [str(x) for x in range(0, max_len + 1, 100)]
    ax.set_xticks(xt_idxs)
    ax.set_xticklabels(xt_labels, rotation=45, ha="right")

    ax.set_xlabel("Fragment length (bp)")
    ax.set_ylabel("Enzyme combination")
    ax.set_title(
        f"Fragment count heatmap — {len(combos)} combinations "
        f"(10 bp bins, log scale)")

    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("Fragment count (log scale)")

    save(fig, plots_dir, "heatmap_fragment_lengths")

# graphs for chromozome distribution
def plot_chrom_distribution(results_dir, combos, plots_dir, n_chroms):
    contig_list = read_contig_lengths_txt(results_dir)
    chroms      = [acc for acc, _ in contig_list[:n_chroms]]
    chrom_idx   = {acc: i for i, acc in enumerate(chroms)}
    n_c         = len(chroms)

    mat = np.zeros((len(combos), n_c), dtype=np.int64)
    for i, combo in enumerate(combos):
        rows = read_csv_rows(combo_file(results_dir, combo, "filtered.csv"))
        for r in rows:
            if r["accession"] in chrom_idx:
                mat[i, chrom_idx[r["accession"]]] += 1

    if mat.max() == 0:
        print(f"  [WARN] chrom plots: no filtered fragments on first "
              f"{n_chroms} contigs, skipping.")
        return

    chrom_labels = [
        f"chr{i+1}\n({chroms[i][:11]}…)" if len(chroms[i]) > 11
        else f"chr{i+1}\n({chroms[i]})"
        for i in range(n_c)
    ]

    # heatmap
    fig_h = max(4, len(combos) * 0.40 + 1.5)
    fig_w = max(8, n_c * 0.55 + 3)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    im = ax.imshow(mat, aspect="auto", cmap="YlOrRd", origin="lower")
    ax.set_yticks(np.arange(len(combos)))
    ax.set_yticklabels(combos,
                       fontsize=max(6, min(9, 120 // max(len(combos), 1))))
    ax.set_xticks(np.arange(n_c))
    ax.set_xticklabels(chrom_labels, rotation=45, ha="right", fontsize=7)
    ax.set_xlabel("Chromosome (contig, longest first)")
    ax.set_ylabel("Enzyme combination")
    ax.set_title(f"Filtered fragment count per chromosome — first {n_c} contigs")
    fig.colorbar(im, ax=ax, shrink=0.8).set_label("Fragment count")
    save(fig, plots_dir, "heatmap_chrom_distribution")

    # grouped bar chart
    fig_w2 = max(10, n_c * max(len(combos), 1) * 0.05 + 3)
    fig2, ax2 = plt.subplots(figsize=(fig_w2, 5))

    x       = np.arange(n_c)
    width   = 0.8 / max(len(combos), 1)
    colours = (plt.cm.tab20 if len(combos) > 10 else plt.cm.tab10)(
        np.linspace(0, 1, len(combos)))

    for j, combo in enumerate(combos):
        offset = (j - len(combos) / 2 + 0.5) * width
        ax2.bar(x + offset, mat[j], width=width * 0.9,
                color=colours[j], label=combo, alpha=0.85)

    ax2.set_xticks(np.arange(n_c))
    ax2.set_xticklabels(chrom_labels, rotation=45, ha="right", fontsize=7)
    ax2.set_xlabel("Chromosome")
    ax2.set_ylabel("Filtered fragment count")
    ax2.set_title(f"Filtered fragments per chromosome — first {n_c} contigs")
    ax2.yaxis.set_major_formatter(
        ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax2.legend(loc="upper right", fontsize=7,
               ncol=max(1, len(combos) // 8),
               bbox_to_anchor=(1.01, 1), borderaxespad=0)
    save(fig2, plots_dir, "bar_chrom_distribution")

# gc content
def plot_gc_distribution(results_dir, combos, plots_dir):
    labels  = []
    means   = []
    medians = []
    stds    = []
    mins_   = []
    maxs_   = []

    for combo in combos:
        gc = read_gc_metrics_csv(
            combo_file(results_dir, combo, "gc_metrics.csv"))
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

    n     = len(labels)
    fig_w = max(8, n * 0.45 + 2)
    fig, ax = plt.subplots(figsize=(fig_w, 5))
    x = np.arange(n)

    ax.bar(x, [s * 2 for s in stds],
           bottom=[m - s for m, s in zip(means, stds)],
           width=0.55, color="#4393c3", alpha=0.6, label="mean ± 1 SD")
    ax.scatter(x, means,   color="#2166ac", zorder=5, s=40, label="mean")
    ax.scatter(x, medians, color="#d6604d", zorder=5, s=40,
               marker="D", label="median")
    for i in range(n):
        ax.plot([x[i], x[i]], [mins_[i], maxs_[i]],
                color="grey", linewidth=1.0, zorder=3)

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=55, ha="right",
                       fontsize=max(6, min(9, 120 // n)))
    ax.set_ylabel("GC content (%)")
    ax.set_title("GC content of filtered fragments — mean ± SD (whiskers = min/max)")
    ax.legend(loc="upper right", fontsize=8)
    ax.yaxis.set_major_formatter(
        ticker.FuncFormatter(lambda v, _: f"{v:.0f}%"))
    save(fig, plots_dir, "gc_distribution")

# annotation
def plot_annotation_coverage(results_dir, combos, plots_dir):
    gff_cats_set, te_cats_set = set(), set()
    data = {}

    for c in combos:
        ann = read_annotation_summary_csv(
            combo_file(results_dir, c, "annotation_summary.csv"))
        if not ann:
            continue
        data[c] = ann
        for cat in ann:
            if cat.startswith("gff:"):
                gff_cats_set.add(cat)
            elif cat.startswith("te:"):
                te_cats_set.add(cat)

    if not data:
        print("  [WARN] annotation plot: no annotation_summary.csv found, skipping.")
        return

    gff_cats     = sorted(gff_cats_set)
    te_cats      = sorted(te_cats_set)
    combos_found = [c for c in combos if c in data]
    n            = len(combos_found)

    has_gff  = bool(gff_cats)
    has_te   = bool(te_cats)
    n_panels = has_gff + has_te

    gff_colors = list(plt.cm.tab20.colors)
    te_colors  = list(plt.cm.Set2.colors)

    fig, axes = plt.subplots(
        1, n_panels,
        figsize=(16, max(4, n * 0.45)),
        sharey=True,
        squeeze=False,
        gridspec_kw={"wspace": 0.02},
    ) 
    axes = axes[0]

    def draw_stacked(ax, cats, colors, title):
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
            # label in % on 10% or more
            for i, (left, val) in enumerate(zip(lefts, vals)):
                if val >= 10.0:
                    ax.text(
                        left + val / 2, i, f"{val:.1f}%",
                        va="center", ha="center", fontsize=7, color="white",
                        fontweight="bold",
                    )

        # unannotated to 100%
        unannotated_in_legend = False
        for i, combo in enumerate(combos_found):
            total     = sum(data[combo].get(cat, 0) for cat in cats)
            remainder = max(0.0, 100.0 - total)
            if remainder > 0:
                label = "unannotated" if not unannotated_in_legend else "_nolegend_"
                ax.barh(i, remainder, left=total,
                        color="#cccccc", label=label)
                if remainder >= 10.0:
                    ax.text(
                        total + remainder / 2, i, f"{remainder:.1f}%",
                        va="center", ha="center", fontsize=7, color="#333333",
                        fontweight="bold",
                    )
                unannotated_in_legend = True

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

# size distribution different
def plot_size_distributions(results_dir, combos, plots_dir, size_low, size_high):
    std_labels = [f"{b}-{b+99}" for b in range(0, 1000, 100)] + [">=1000"]
    x_pos      = np.arange(len(std_labels))
    n          = len(combos)
    colours    = (plt.cm.tab20 if n > 10 else plt.cm.tab10)(
        np.linspace(0, 1, max(n, 1)))

    fig_w = max(12, len(std_labels) * 0.7 + 3)
    fig, ax = plt.subplots(figsize=(fig_w, 5))

    for i, combo in enumerate(combos):
        dist = read_distribution_csv(
            combo_file(results_dir, combo, "distribution.csv"))
        vals = [dist.get(l, 0) for l in std_labels]
        ax.plot(x_pos, vals, marker="o", markersize=4, linewidth=1.5,
                color=colours[i % len(colours)], label=combo, alpha=0.85)

    def _bin_idx(bp):
        return len(std_labels) - 1 if bp >= 1000 else bp // 100

    ax.axvspan(
        max(0, _bin_idx(size_low) - 0.5),
        min(len(std_labels) - 1, _bin_idx(size_high) + 0.5),
        color="#f4a261", alpha=0.25,
        label=f"selected window {size_low}–{size_high} bp")

    ax.set_xticks(x_pos)
    ax.set_xticklabels(std_labels, rotation=45, ha="right")
    ax.set_xlabel("Fragment length bin")
    ax.set_ylabel("Fragment count")
    ax.set_title("Fragment length distribution per combination")
    ax.yaxis.set_major_formatter(
        ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.legend(loc="upper right", fontsize=7,
              ncol=max(1, n // 12),
              bbox_to_anchor=(1.01, 1), borderaxespad=0)
    save(fig, plots_dir, "size_distributions")


def parse_args():
    p = argparse.ArgumentParser(
        prog="rendogbs_plots.py",
        formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument("--workdir", required=True)
    p.add_argument("--size", required=True)
    p.add_argument("--chroms", type=int, required=True)
    p.add_argument("--dpi", type=int, default=300)
    return p.parse_args()

def main():
    args = parse_args()

    global _DPI
    _DPI = args.dpi

    try:
        size_low, size_high = parse_size_range(args.size)
    except ValueError as e:
        print(f"[ERROR] {e}"); sys.exit(1)

    results_dir = os.path.join(args.workdir, "results")
    plots_dir   = os.path.join(results_dir, "plots")

    if not os.path.isdir(results_dir):
        print(f"[ERROR] results directory not found: {results_dir}")
        print(f"        Run rendogbs_pipeline.py first.")
        sys.exit(1)

    combos = discover_combos(results_dir)
    if not combos:
        print(f"[ERROR] No combination sub-directories found in {results_dir}")
        sys.exit(1)

    print(f"\n[rendogbs_plots]")
    print(f"  results dir : {results_dir}")
    print(f"  combinations: {len(combos)}")
    print(f"  size window : {size_low}-{size_high} bp")
    print(f"  chromosomes : {args.chroms}")
    print(f"  dpi         : {args.dpi}\n")

    print("[1/5] Fragment-length heatmap ...")
    plot_heatmap(results_dir, combos, plots_dir)

    print(f"[2/5] Per-chromosome distribution "
          f"(first {args.chroms} contigs as chromosomes) ...")
    plot_chrom_distribution(results_dir, combos, plots_dir, args.chroms)

    print("[3/5] GC content distribution ...")
    plot_gc_distribution(results_dir, combos, plots_dir)

    print("[4/5] Annotation coverage ...")
    plot_annotation_coverage(results_dir, combos, plots_dir)

    print("[5/5] Size distribution line plot ...")
    plot_size_distributions(results_dir, combos, plots_dir, size_low, size_high)

if __name__ == "__main__":
    main()