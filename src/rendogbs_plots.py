#!/usr/bin/env python3
"""
rendogbs_plots.py  -  Post-processing: TSV summary rebuild + plots

Run after rendogbs_pipeline.py has finished (or call via rendogbs.sh which
sequences the two steps automatically).

Reads  : <workdir>/results/  (per-combination sub-directories + contig_lengths.txt)
Writes : <workdir>/results/run_summary.tsv   (rebuilt from CSVs)
         <workdir>/results/plots/            (all figures as PNG)

Usage:
    python rendogbs_plots.py --workdir ./run1 --size 200-400 [--chroms 24] [--dpi 300]

Figures produced:
    heatmap_fragment_lengths.png    fragment count heatmap, all combos x 10 bp bins
    heatmap_chrom_distribution.png  filtered fragment count: combos x chromosomes
    bar_chrom_distribution.png      grouped bar chart version of the above
    gc_distribution.png             GC% summary per combination (mean ± SD)
    annotation_coverage.png         stacked bar of annotation category coverage
    size_distributions.png          line plot of 100 bp bins with user window shaded
"""

import os
import csv
import argparse
import sys
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg") # no GUI
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import LogNorm
import matplotlib.patches as mpatches

# globalni styl grafu
PALETTE = "viridis"
FIG_EXT = "png"
_DPI    = 300   # jde zmenit --dpi

plt.rcParams.update({
    "font.family":    "DejaVu Sans",
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
})

#  Pomocna funkce - oddeleni range ktere analyzujeme pro knihovnu (200, 400)
def parse_size_range(s):
    parts = s.split("-")
    return int(parts[0]), int(parts[1])

# vraci list serazenych kombinaci enzymu z nazev podslozek
def discover_combos(results_dir):
    combos = []
    for entry in sorted(os.scandir(results_dir), key=lambda e: e.name):
        if entry.is_dir() and os.path.isfile(
                os.path.join(entry.path, "fragments.csv")):
            combos.append(entry.name) # kandidatni podslozka je ta co ma soubor fragments.csv
    return combos

# nacteni souboru a vraceni seznamu radku, kontrola existence souboru, vraci vsechny radky jako seznam slovniku {"sloupec1":"hodnota1", ...},{"sloupec1":"hodnota2", ...} DictReader prevod kazdeho radku na slovnik
#
# Creates a closure counting how many times the inner function was
# called. The inner function returns any value unchanged and
# increments count. If count modulo period is zero, prints count as a
# side-effect.
def make_idx_printer(path, period):
    idx = 0
    def the_printer(v):
        nonlocal idx
        if (idx % period) == 0:
            print(f"{path}: {idx}")
        idx = idx + 1
        return v
    return the_printer

def read_csv_rows(path):
    if not os.path.isfile(path):
        return [] # neexistujici soubor
    with open(path, newline="") as fh:
        # nacteni distribution.csv a prevede na slovnik z csv po binech, preskakuje poskozene radky
        logidx = make_idx_printer(path, 10000)
        return [logidx(d) for d in csv.DictReader(fh)]

def read_distribution_csv(path):
    d = {}
    lst = read_csv_rows(path)
    print(lst)
    for row in lst:
        try:
            print(row)
            d[row["length_range"]] = int(row["count"])
        except (KeyError, ValueError):
            pass
    return d

# nacteni gc_metrics.csv a prevede na slovnik
def read_gc_metrics_csv(path):
    return {row["metric"]: row["value"] for row in read_csv_rows(path)
            if "metric" in row}

# nacteni annotation_summary.csv a prevede na slovnik, preskakuje chyby opet, gene, SINE, LINE,...
def read_annotation_summary_csv(path):
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

# nacteni contig_lengths.txt z results vraci list serazenych radku
def read_contig_lengths_txt(results_dir):
    path = os.path.join(results_dir, "contig_lengths.txt")
    rows = []
    with open(path) as fh:
        next(fh) # preskoceni hlavicky
        for line in fh:
            parts = line.rstrip().split("\t") # odstraneni konce radku
            if len(parts) == 2: # # ocekavame dva sloupce
                rows.append((parts[0], int(parts[1]))) # accession a delka vypsani
    return rows

# ukladani matplotlib figure, plots dir - plots, jmeno souboru bez pripony - stejna na zacatku specifikovana
def save(fig, plots_dir, name):
    os.makedirs(plots_dir, exist_ok=True)
    path = os.path.join(plots_dir, f"{name}.{FIG_EXT}")
    fig.savefig(path, dpi=_DPI, bbox_inches="tight") # tight = nezobrazovat mezery, ulozeni obrazku
    plt.close(fig) # zavreni figure
    print(f"  [plot] saved -> {path}")
    return path

# run_summary.csv jeden radek jedna kombinace enzymu
def build_run_summary(results_dir, combos, size_low, size_high):
    all_ann_cats = set() # mnozina vsech kategorii z anotaci (rm.out, gff)
    for combo in combos:
        ann = read_annotation_summary_csv(
            os.path.join(results_dir, combo, "annotation_summary.csv")) # nacteni anotacniho souboru pro konkretni kombinaci
        all_ann_cats.update(ann.keys()) # pridani kategorii nalezenych
    all_ann_cats = sorted(all_ann_cats) # serazeni abecedne kategorii

    std_labels   = [f"{b}-{b+99}" for b in range(0, 1000, 100)] + [">=1000"] # standardni biny
    custom_lbl   = f"custom_{size_low}-{size_high}" # custom bin s velikosti knihovny

    header = (
        ["combination", "total_cuts", "total_fragments",
         f"filtered_{custom_lbl}"]
        + [f"bin_{l}" for l in std_labels]
        + [f"bin_{custom_lbl}"]
        + ["gc_mean_pct", "gc_median_pct", "gc_std_pct",
           "gc_min_pct", "gc_max_pct", "gc_n_fragments"]
        + [f"ann_pct_{c}" for c in all_ann_cats]
    ) # hlavicka, gc statistiky, anotacni kategorie, pocty fragmentu, notacni procenta

    rows = [] #vytvorim radky zde ukladat data
    for combo in combos: #jeden pruchod, jeden radek tsv
        d      = os.path.join(results_dir, combo) #cesta k kombinaci
        cuts   = read_csv_rows(os.path.join(d, "cuts.csv"))
        frags  = read_csv_rows(os.path.join(d, "fragments.csv"))
        filt   = read_csv_rows(os.path.join(d, "filtered.csv"))
        dist   = read_distribution_csv(os.path.join(d, "distribution.csv"))
        gc     = read_gc_metrics_csv(os.path.join(d, "gc_metrics.csv"))
        ann    = read_annotation_summary_csv(os.path.join(d, "annotation_summary.csv"))

        rows.append( #pridani radku jednoho
            [combo, len(cuts), len(frags), len(filt)]
            + [dist.get(l, 0) for l in std_labels] # pocet fragmentu v binech
            + [len(filt)] # pocet fragmentu v custom binu
            + [gc.get("gc_mean_pct",   "n/a"),
               gc.get("gc_median_pct", "n/a"),
               gc.get("gc_std_pct",    "n/a"),
               gc.get("gc_min_pct",    "n/a"),
               gc.get("gc_max_pct",    "n/a"),
               gc.get("n_fragments",   0)] #GC statistiky
            + [ann.get(c, "n/a") for c in all_ann_cats] # anotacni procenta
        )

    tsv_path = os.path.join(results_dir, "run_summary.tsv")
    with open(tsv_path, "w", newline="") as fh: # otevreni souboru pro zapis
        w = csv.writer(fh, delimiter="\t") # nastaveni oddelovacu
        w.writerow(header) # zapis hlavicky
        w.writerows(rows) # zapis radku

    print(f"  [tsv]  run_summary.tsv -> {tsv_path}")


#  Fragment-length heatmap  (all fragments, 10 bp bins, log scale). Source: fragments.csv (all fragments, not just filtered)

def plot_heatmap(results_dir, combos, plots_dir): #zadani parametru - slozka s data, kombinace enzymu, slozka pro ulozeni obrazku
    bin_width = 10
    max_len   = 1000
    edges     = list(range(0, max_len + 1, bin_width)) + [float("inf")] #vytvoreni hranic binu, inf aby se nezahodili
    n_bins    = len(edges) - 1 #pocet binu

    mat = np.zeros((len(combos), n_bins), dtype=np.int64) #matice pro heatmapu s nulami pro pocet radku a sloupce

    for i, combo in enumerate(combos): # tvorba array pro heatmapu z fragments.csv po kombinacich
        rows = read_csv_rows(os.path.join(results_dir, combo, "fragments.csv"))
        lengths = np.array([int(r["fragment_length"]) for r in rows],
                           dtype=np.int64)
        if lengths.size == 0:
            continue
        clipped         = np.minimum(lengths, max_len) # oriznuti nejdelsich fragmentu nad 1000 jako 1000
        counts, _       = np.histogram(clipped, bins=edges) # pocet fragmentu v jednotlivych binech do array matice
        mat[i]          = counts # ulozeni do matice

    if mat.max() == 0: # nalezeni nejvetsiho cisla v matici pokud je nula nejsou data pro vytvoreni grafu
        print("  [WARN] heatmap: no fragment data found, skipping.")
        return 

    fig_h = max(4, len(combos) * 0.40 + 1.5) # vysku grafu zvysuji s poctem kombinaci
    fig_w = max(12, n_bins * 0.12 + 3) # sirku grafu zvysuji s poctem binu, mela by byt konstantni
    fig, ax = plt.subplots(figsize=(fig_w, fig_h)) # matplotlib fig cely obrazek a ax graf

    im = ax.imshow( #vykresleni matice na heatmapu
        mat, aspect="auto", cmap=PALETTE, #aspect automaticka velikost bunek
        norm=LogNorm(vmin=1, vmax=max(int(mat.max()), 1)), # logaritmicka skala zvyrazneni i malych poctu
        origin="lower", # zobrazeni prvni kombinace dole
    )

    ax.set_yticks(np.arange(len(combos))) # popsany y osy kombinacemi, pomoci array kazda v radku kombinace
    ax.set_yticklabels(combos,
                       fontsize=max(6, min(9, 120 // max(len(combos), 1)))) # dynamicke zmenseni velikosti pisma dle poctu kombinaci

    xt_idxs   = np.arange(0, max_len + 1, 100) // bin_width # indexy pro x osu array, po 100, ne po binech
    xt_labels = [str(x) for x in range(0, max_len + 1, 100)] # vytvori cisla po 100
    ax.set_xticks(xt_idxs) # popsana x osa
    ax.set_xticklabels(xt_labels, rotation=45, ha="right") # rotace a zarovnani

    ax.set_xlabel("Fragment length (bp)") # popis x osy
    ax.set_ylabel("Enzyme combination") # popis y osy
    ax.set_title(
        f"Fragment count heatmap — {len(combos)} combinations "
        f"(10 bp bins, log scale)") # titulek grafu

    cbar = fig.colorbar(im, ax=ax, shrink=0.8) # barevna legenda
    cbar.set_label("Fragment count (log scale)") # popis barevne legendy log skaly

    save(fig, plots_dir, "heatmap_fragment_lengths") # ulozeni obrazku

# chromozomovy graf - zobrazeni mnozstvi fragmentu na jednotlivych chromozomech; grouped bar graf pro distribuci fragmentu po chromozomech
def plot_chrom_distribution(results_dir, combos, plots_dir, n_chroms): #zadani parametru - slozka s data, kombinace enzymu, slozka pro ulozeni obrazku, pocty chromozomu
    contig_list = read_contig_lengths_txt(results_dir)

    chroms    = [acc for acc, _ in contig_list[:n_chroms]] # beru pouze x nejdelsich beru jako chromozomu
    chrom_idx = {acc: i for i, acc in enumerate(chroms)} # indexy chromozomu 0-based
    n_c       = len(chroms) # pocet chromozomu

    mat = np.zeros((len(combos), n_c), dtype=np.int64) # matice pro ulozeni dat nuly
    for i, combo in enumerate(combos): # pro kazdou kombinaci nahradim nuly poctem
        rows = read_csv_rows(
            os.path.join(results_dir, combo, "filtered.csv")) # nacitam soubor s filtrovanymi fragmenty
        for r in rows: # jedu po radcich 
            if r["accession"] in chrom_idx: # jestli je chromozom v seznamu a neni to nenamapovany contig
                mat[i, chrom_idx[r["accession"]]] += 1 # najdi radek a sloupec a pricti 1

    if mat.max() == 0: # pokud je nejvyssi cislo 0 preskocit a vytisknout chybovou hlasku
        print(f"  [WARN] chrom plots: no filtered fragments on first "
              f"{n_chroms} contigs, skipping.")
        return

    chrom_labels = [
        f"chr{i+1}\n({chroms[i][:11]}…)" if len(chroms[i]) > 11 #zkraceni dlouhych nazvu 
        else f"chr{i+1}\n({chroms[i]})"
        for i in range(n_c)
    ]

    fig_h = max(4, len(combos) * 0.40 + 1.5) # vysku grafu zvysuji s poctem kombinaci
    fig_w = max(8, n_c * 0.55 + 3) # sirku grafu zvysuji s poctem binu, mela by byt konstantni
    fig, ax = plt.subplots(figsize=(fig_w, fig_h)) # matplotlib fig cely obrazek a ax graf

    im = ax.imshow(mat, aspect="auto", cmap="YlOrRd", origin="lower") # barevna skala, prvni kombinace dole 
    ax.set_yticks(np.arange(len(combos)))
    ax.set_yticklabels(combos,
                       fontsize=max(6, min(9, 120 // max(len(combos), 1))))
    ax.set_xticks(np.arange(n_c)) # chromozomy na ose x
    ax.set_xticklabels(chrom_labels, rotation=45, ha="right", fontsize=7) # rotace nazvu, zarovnani
    ax.set_xlabel("Chromosome (contig, longest first)") # popis x osy
    ax.set_ylabel("Enzyme combination") # popis y osy
    ax.set_title(f"Filtered fragment count per chromosome — first {n_c} contigs") # popis grafu dynamicky počet chromozomu
    fig.colorbar(im, ax=ax, shrink=0.8).set_label("Fragment count") # barevna legenda
    save(fig, plots_dir, "heatmap_chrom_distribution") # ulozeni obrazku

    # grouped bar graf - stejna data jina vizualizace
    fig_w2 = max(10, n_c * max(len(combos), 1) * 0.05 + 3)
    fig2, ax2 = plt.subplots(figsize=(fig_w2, 5))

    x       = np.arange(n_c) # pozice chromozomu
    width   = 0.8 / max(len(combos), 1) # sirka sloupcu podle poctu kombinaci
    colours = (plt.cm.tab20 if len(combos) > 10 else plt.cm.tab10)(np.linspace(0, 1, len(combos))) # barvy podle poctu kombinaci automaticky

    for j, combo in enumerate(combos): # pro kazdou kombinaci
        offset = (j - len(combos) / 2 + 0.5) * width # zarovnani sloupcu vedle sebe
        ax2.bar(x + offset, mat[j], width=width * 0.9, color=colours[j], label=combo, alpha=0.85) # vykresleni sloupce pro kazdou kombinaci

    ax.set_xticks(np.arange(n_c))
    ax.set_xticklabels(chrom_labels, rotation=45, ha="right", fontsize=7) # rotace a zarovnani
    ax2.set_xlabel("Chromosome") # popisek x osy
    ax2.set_ylabel("Filtered fragment count") # popisek y osy
    ax2.set_title(f"Filtered fragments per chromosome — first {n_c} contigs") # nazev grafu, dynamicky pocet chromozomu
    ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax2.legend(loc="upper right", fontsize=7, ncol=max(1, len(combos) // 8), bbox_to_anchor=(1.01, 1), borderaxespad=0) # ukotveni legendy a vytvoreni
    save(fig2, plots_dir, "bar_chrom_distribution") # ulozeni obrazku

# GC content distribution - One bar per combination: mean GC% ± 1 SD, whiskers = min/max, diamond = median. Reads gc_metrics.csv (summary stats computed by the pipeline for filtered fragments only).

def plot_gc_distribution(results_dir, combos, plots_dir):
    labels  = []
    means   = []
    medians = []
    stds    = []
    mins_   = []
    maxs_   = []

    for combo in combos:
        gc = read_gc_metrics_csv(
            os.path.join(results_dir, combo, "gc_metrics.csv"))
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
    ax.set_title(
        "GC content of filtered fragments — mean ± SD (whiskers = min/max)")
    ax.legend(loc="upper right", fontsize=8)
    ax.yaxis.set_major_formatter(
        ticker.FuncFormatter(lambda v, _: f"{v:.0f}%"))
    save(fig, plots_dir, "gc_distribution")

#  Annotation coverage stacked horizontal bar. categories (TE classes, gene features), values = % of filtered fragment bases overlapping that category.
def plot_annotation_coverage(results_dir, combos, plots_dir):
    all_cats  = []
    combo_ann = {}

    for combo in combos:
        ann = read_annotation_summary_csv(
            os.path.join(results_dir, combo, "annotation_summary.csv"))
        if ann:
            combo_ann[combo] = ann
            for c in ann:
                if c not in all_cats:
                    all_cats.append(c)

    if not combo_ann:
        print("  [WARN] annotation plot: no annotation_summary.csv found, "
              "skipping.")
        return

    all_cats = sorted(all_cats)
    active   = [c for c in combos if c in combo_ann]
    n        = len(active)
    colours  = plt.cm.tab20(np.linspace(0, 1, max(len(all_cats), 1)))
    cat_col  = {c: colours[i] for i, c in enumerate(all_cats)}

    fig_h = max(4, n * 0.45 + 1.5)
    fig, ax = plt.subplots(figsize=(12, fig_h))

    for i, combo in enumerate(active):
        ann  = combo_ann[combo]
        left = 0.0
        for cat in all_cats:
            val = ann.get(cat, 0.0)
            if val > 0:
                ax.barh(i, val, left=left, color=cat_col[cat],
                        height=0.65,
                        label=cat if i == 0 else "")
                if val >= 1.5:
                    ax.text(left + val / 2, i, f"{val:.1f}",
                            va="center", ha="center", fontsize=6,
                            color="white")
                left += val

    ax.set_yticks(np.arange(n))
    ax.set_yticklabels(active,
                       fontsize=max(6, min(9, 120 // max(n, 1))))
    ax.set_xlabel("% of filtered fragment bases")
    ax.set_title("Annotation coverage of filtered fragments")
    ax.set_xlim(0, 105)
    ax.xaxis.set_major_formatter(
        ticker.FuncFormatter(lambda v, _: f"{v:.0f}%"))

    handles = [mpatches.Patch(color=cat_col[c], label=c) for c in all_cats]
    ax.legend(handles=handles, loc="lower right", fontsize=7,
              ncol=max(1, len(all_cats) // 12),
              bbox_to_anchor=(1.01, 0), borderaxespad=0)
    save(fig, plots_dir, "annotation_coverage")

# Size distribution line plot  (100 bp bins, user window highlighted). One line per combination across standard 100 bp bins. The user's size window is highlighted as a shaded region so it is immediately obvious how many fragments fall inside vs. outside the selected gel-extraction window.
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
            os.path.join(results_dir, combo, "distribution.csv"))
        vals = [dist.get(l, 0) for l in std_labels]
        ax.plot(x_pos, vals, marker="o", markersize=4, linewidth=1.5,
                color=colours[i % len(colours)], label=combo, alpha=0.85)

    # Shade user window
    def _bin_idx(bp):
        return len(std_labels) - 1 if bp >= 1000 else bp // 100

    ax.axvspan(max(0, _bin_idx(size_low) - 0.5),
               min(len(std_labels) - 1, _bin_idx(size_high) + 0.5),
               color="#f4a261", alpha=0.25,
               label=f"selected window {size_low}-{size_high} bp")

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


#  Argument parser
def parse_args():
    p = argparse.ArgumentParser(
        prog="rendogbs_plots.py",
        description=(
            "rendogbs_plots - Post-processing for rendogbs pipeline output.\n\n"
            "Rebuilds run_summary.tsv from per-combination CSV files and\n"
            "generates six diagnostic figures saved to "
            "<workdir>/results/plots/."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Figures produced:
  heatmap_fragment_lengths.png    fragment count heatmap (all combos x 10 bp bins,
                                  log scale) -- matches your original heatmap script
  heatmap_chrom_distribution.png  filtered fragments per chromosome (heatmap)
  bar_chrom_distribution.png      filtered fragments per chromosome (bar chart)
  gc_distribution.png             GC% per combination: mean ± SD, whiskers = min/max
  annotation_coverage.png         annotation category coverage stacked bar
  size_distributions.png          100 bp bin line plot with user window highlighted

Examples:
  python rendogbs_plots.py --workdir ./run1 --size 200-400
  python rendogbs_plots.py --workdir ./run1 --size 150-350 --chroms 24 --dpi 600
  python rendogbs_plots.py --workdir ./run1 --size 200-400 --skip-plots
""")

    p.add_argument("--workdir",
                   required=True,
                   help="rendogbs working directory (results/ sub-dir is read/written)")
    p.add_argument("--size",
                   required=True,
                   help="Size window used in the pipeline run, e.g. 200-400")
    p.add_argument("--chroms",
                   type=int, default=10,
                   help="Number of longest contigs treated as chromosomes (default: 10)")
    p.add_argument("--dpi",
                   type=int, default=300,
                   help="Output figure resolution in DPI (default: 300)")
    p.add_argument("--skip-tsv",
                   action="store_true",
                   help="Skip rebuilding run_summary.tsv")
    p.add_argument("--skip-plots",
                   action="store_true",
                   help="Rebuild TSV only, do not generate any figures")
    return p.parse_args()


#  Main
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

    # TSV 
    if not args.skip_tsv:
        print("[1/6] Rebuilding run_summary.tsv ...")
        build_run_summary(results_dir, combos, size_low, size_high)
    else:
        print("[1/6] TSV skipped (--skip-tsv)")

    if args.skip_plots:
        print("Plots skipped (--skip-plots). Done.")
        return

    # Fragment-length heatmap
    print("[2/6] Fragment-length heatmap ...")
    plot_heatmap(results_dir, combos, plots_dir)

    # Per-chromosome distribution 
    print(f"[3/6] Per-chromosome distribution "
          f"(first {args.chroms} contigs as chromosomes) ...")
    plot_chrom_distribution(results_dir, combos, plots_dir, args.chroms)

    # GC distribution
    print("[4/6] GC content distribution ...")
    plot_gc_distribution(results_dir, combos, plots_dir)

    # Annotation coverage 
    print("[5/6] Annotation coverage ...")
    plot_annotation_coverage(results_dir, combos, plots_dir)

    # Size distribution line plot 
    print("[6/6] Size distribution line plot ...")
    plot_size_distributions(results_dir, combos, plots_dir, size_low, size_high)

    print(f"""
+--------------------------------------------------+
|  rendogbs_plots - DONE                           |
+--------------------------------------------------+
|  Figures written to:                             |
|  {plots_dir:<48} |
|                                                  |
|  heatmap_fragment_lengths.png                    |
|  heatmap_chrom_distribution.png                  |
|  bar_chrom_distribution.png                      |
|  gc_distribution.png                             |
|  annotation_coverage.png                         |
|  size_distributions.png                          |
+--------------------------------------------------+
""")


if __name__ == "__main__":
    main()
