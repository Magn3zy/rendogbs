#!/usr/bin/env python3
"""
rendogbs_pipeline.py  -  In-silico restriction digest & GBS library simulator

Invoked via wrapper : rendogbs.sh
Enzymes loaded from : endonucleases.py  (must be in the same directory)
Combinations via    : --combinations fast | all | custom  (+ --combinations-file)

Output structure:
    workdir/
        results/
            contig_lengths.txt          all contigs sorted descending by length
            run_summary.tsv             Excel-ready: one row per enzyme combo
            EcoRI_MseI/
                cuts.csv                all cut sites
                fragments.csv           all adjacent pairs of different enzymes
                filtered.csv            fragments within the user size window
                distribution.csv        length distribution (std bins + custom)
                gc_metrics.csv          GC statistics for filtered fragments
                annotation_summary.csv  TE / gene coverage (optional)
            AciI_HindIII/
                ...
"""

#zkusit dependency graph - aby to delala pbs kdyztak
import os
import sys
import csv
import gzip
import re
import bisect as _bisect
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict

#LOAD ENZYME DICTIONARY

_script_dir = os.path.dirname(os.path.abspath(__file__))
_cwd        = os.getcwd()
for _search_dir in (_script_dir, _cwd):
    if _search_dir not in sys.path:
        sys.path.insert(0, _search_dir)

ENDONUCLEASES = None

# Try all known module/variable name combinations so the pipeline works
# with both the original endonucleases.py and the new restriction_enzymes.py
_candidates = [
    ("endonucleases",       "endonucleases"),
    ("restriction_enzymes", "restriction_enzymes"),
    ("endonucleases",       "restriction_enzymes"),
    ("restriction_enzymes", "endonucleases"),
]
for _mod, _var in _candidates:
    try:
        _m = __import__(_mod)
        _candidate = getattr(_m, _var, None)
        if isinstance(_candidate, dict) and len(_candidate) > 0:
            ENDONUCLEASES = _candidate
            break
    except ImportError:
        continue

if ENDONUCLEASES is None:
    print("[ERROR] Enzyme dictionary not found.")
    print(f"        Looked in: {_script_dir}")
    if _cwd != _script_dir:
        print(f"                   {_cwd}")
    print("        Expected: endonucleases.py  (variable: endonucleases)")
    print("               or restriction_enzymes.py  (variable: restriction_enzymes)")
    sys.exit(1)

# IUPAC AMBIGUITY CODES  ->  regex character classes. Only used for enzymes whose recognition sequence contains ambiguous bases. Standard A/C/G/T motifs use the faster str.find path.
_IUPAC = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'W': '[AT]',  'R': '[AG]',  'Y': '[CT]',  'S': '[CG]',
    'K': '[GT]',  'M': '[AC]',  'B': '[CGT]', 'D': '[AGT]',
    'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]',
}
_STANDARD = set("ACGT")


def _is_degenerate(motif):
    return any(c not in _STANDARD for c in motif.upper())


def _motif_to_regex(motif):
    return re.compile("".join(_IUPAC.get(c, c) for c in motif.upper()))

#Return (name, motif, cut_offset, motif_len, is_degenerate, compiled_regex). compiled_regex is None for plain ACGT motifs (str.find used instead).
def build_enzyme_entry(name):
    data    = ENDONUCLEASES[name]
    motif   = data["sequence"].upper()
    cut_off = data["cut"]
    degen   = _is_degenerate(motif)
    return (name, motif, cut_off, len(motif), degen,
            _motif_to_regex(motif) if degen else None)

#  STANDARD SIZE BINS  (always reported in distribution.csv)
_STD_STARTS = list(range(0, 1000, 100))
_STD_LABELS = [f"{b}-{b+99}" for b in _STD_STARTS]   # "0-99" .. "900-999"


def std_size_label(fl):
    """Return the standard 100 bp bin label, or '>=1000'."""
    if fl >= 1000:
        return ">=1000"
    return f"{(fl // 100) * 100}-{(fl // 100) * 100 + 99}"


# CUSTOM SIZE RANGE PARSING. Parse 'LOW-HIGH' (e.g. '200-400').  Both ends inclusive.

def parse_size_range(size_str):
    parts = size_str.split("-")
    if len(parts) != 2:
        raise ValueError(
            f"Invalid size range '{size_str}'. Use format LOW-HIGH, e.g. 200-400")
    low, high = int(parts[0]), int(parts[1])
    if low >= high:
        raise ValueError(
            f"Size range LOW ({low}) must be less than HIGH ({high})")
    return low, high


def in_custom_range(fl, low, high):
    return low <= fl <= high

#ENZYME COMBINATION PRESETS

FAST_COMBOS = [
    ("EcoRI",   "MseI"),
    ("EcoRI",   "MboI"),
    ("EcoRI",   "NlaIII"),
    ("PstI",    "MspI"),
    ("PstI",    "MseI"),
    ("PstI",    "TaqI-v2"),
    ("PstI",    "MboI"),
    ("PstI",    "HpaII"),
    ("SbfI",    "MspI"),
    ("SbfI",    "MseI"),
    ("ApeKI",   "MseI"),
    ("BamHI",   "MseI"),
    ("BamHI",   "MboI"),
    ("HindIII", "MseI"),
    ("HindIII", "NlaIII"),
    ("NheI", "MboI"),
    ("NheI", "MseI"),
    ("XbaI",    "MboI"),
    ("SpeI", "MboI"),
    ("KpnI", "MboI"),
    ("NcoI",    "MseI"),
    ("NcoI",    "MboI"),
    ("BglII",   "MboI"),
    ("BglII",   "MseI"),
    ("ClaI",    "MboI"),
]


def load_combinations(mode, custom_path=None):
    if mode == "fast":
        combos = list(FAST_COMBOS)
    elif mode == "all":
        names  = list(ENDONUCLEASES.keys())
        combos = [(a, b) for a in names for b in names if a != b]
    elif mode == "custom":
        if not custom_path:
            print("[ERROR] --combinations custom requires --combinations-file <path>")
            sys.exit(1)
        combos = _load_combinations_txt(custom_path)
    else:
        print(f"[ERROR] Unknown combinations mode: {mode!r}")
        sys.exit(1)

    valid = []
    for ea, eb in combos:
        missing = [e for e in (ea, eb) if e not in ENDONUCLEASES]
        if missing:
            print(f"  [WARN] Enzyme(s) {missing} not in dictionary, "
                  f"skipping {ea}+{eb}")
            continue
        valid.append((ea, eb))
    return valid


def _load_combinations_txt(txt_path):
    combos = []
    with open(txt_path) as fh:
        for i, line in enumerate(fh):
            line = line.strip()
            if not line or i == 0:
                continue
            parts = line.split(",")
            if len(parts) != 2:
                print(f"  [WARN] Line {i+1} bad format, skipping: {line!r}")
                continue
            combos.append((parts[0].strip(), parts[1].strip()))
    return combos

#FASTA  -  plain + gzip

def open_file(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

# Read FASTA (plain or .gz). Returns list of (accession, sequence_uppercase). Accession = first whitespace-delimited token after '>'.
def read_fasta_contigs(path):
    contigs = []
    header  = None
    parts   = []
    with open_file(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    contigs.append((header, "".join(parts).upper()))
                header = line[1:].split()[0]
                parts  = []
            else:
                parts.append(line)
    if header is not None:
        contigs.append((header, "".join(parts).upper()))
    return contigs


def sort_contigs_by_length(contigs):
    return sorted(contigs, key=lambda x: len(x[1]), reverse=True)


def write_contig_lengths_txt(contigs_sorted, results_dir):
    out_path = os.path.join(results_dir, "contig_lengths.txt")
    with open(out_path, "w") as fh:
        fh.write("accession\tlength_bp\n")
        for acc, seq in contigs_sorted:
            fh.write(f"{acc}\t{len(seq)}\n")
    print(f"  [INFO] contig_lengths.txt -> {out_path}")
    return out_path

#PROSTOR pro optimalizaci - KMP - CUT-SITE DETECTION - Forward strand only.  Cut offset from dictionary. Degenerate motifs -> regex; plain ACGT -> str.find (faster). Find all cut positions for two enzymes on the forward strand. Returns sorted list of (abs_cut_position, enzyme_name).

def find_cuts_two_enzymes(seq, entry_a, entry_b):
    hits = []
    for (name, motif, cut_off, motif_len, degen, compiled) in (entry_a, entry_b):
        if degen:
            for m in compiled.finditer(seq):
                hits.append((m.start() + cut_off, name))
        else:
            start = 0
            while True:
                pos = seq.find(motif, start)
                if pos == -1:
                    break
                hits.append((pos + cut_off, name))
                start = pos + 1
    hits.sort(key=lambda x: x[0])
    return hits

#Consecutive cut-site pairs by DIFFERENT enzymes. These fragments carry a different adapter on each end -> selectively enriched. Returns list of (start_pos, start_enzyme, end_enzyme, fragment_length).
def find_adjacent_pairs(cuts):
    pairs = []
    for (p1, e1), (p2, e2) in zip(cuts, cuts[1:]):
        if e1 != e2:
            pairs.append((p1, e1, e2, p2 - p1))
    return pairs

#GC CONTENT
def gc_content(seq):
    if not seq:
        return 0.0
    return (seq.count("G") + seq.count("C")) / len(seq)

#SINGLE-CONTIG PROCESSING

def process_contig(acc, seq, entry_a, entry_b, size_low, size_high):
    cuts  = find_cuts_two_enzymes(seq, entry_a, entry_b)
    frags = find_adjacent_pairs(cuts)

    std_dist = defaultdict(int)
    filtered = []
    for sp, se, ee, fl in frags:
        std_dist[std_size_label(fl)] += 1
        if in_custom_range(fl, size_low, size_high):
            filtered.append((acc, sp, se, ee, fl))

    return {
        "acc":      acc,
        "cuts":     [(acc, pos, enz) for pos, enz in cuts],
        "frags":    [(acc, sp, se, ee, fl) for sp, se, ee, fl in frags],
        "filtered": filtered,
        "std_dist": dict(std_dist),
    }



# SINGLE ENZYME COMBINATION  -  parallel batches of contigs. Contigs sorted longest-first, processed in batches of n_parallel. Each contig in a batch runs in its own thread. cuts.csv  fragments.csv  filtered.csv  distribution.csv  gc_metrics.csv
# zde zvazit optimalizace RAM - otestovat vetsi genomy zaznamenat peaky
def process_combination(enz_a_name, enz_b_name, contigs_sorted,
                        combo_dir, size_low, size_high, n_parallel, ref_seqs):

    os.makedirs(combo_dir, exist_ok=True)

    entry_a = build_enzyme_entry(enz_a_name)
    entry_b = build_enzyme_entry(enz_b_name)

    all_cuts     = []
    all_frags    = []
    all_filtered = []
    total_std    = defaultdict(int)

    batches = [contigs_sorted[i : i + n_parallel]
               for i in range(0, len(contigs_sorted), n_parallel)]

    for batch in batches:
        with ThreadPoolExecutor(max_workers=len(batch)) as executor:
            futures = {
                executor.submit(
                    process_contig, acc, seq,
                    entry_a, entry_b, size_low, size_high
                ): acc
                for acc, seq in batch
            }
            for future in as_completed(futures):
                r = future.result()
                all_cuts.extend(r["cuts"])
                all_frags.extend(r["frags"])
                all_filtered.extend(r["filtered"])
                for lab, cnt in r["std_dist"].items():
                    total_std[lab] += cnt

    key = lambda x: (x[0], x[1])
    all_cuts.sort(key=key)
    all_frags.sort(key=key)
    all_filtered.sort(key=key)

    # cuts.csv
    with open(os.path.join(combo_dir, "cuts.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["accession", "cut_position", "enzyme"])
        w.writerows(all_cuts)

    # fragments.csv
    with open(os.path.join(combo_dir, "fragments.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["accession", "start_pos", "start_enzyme",
                    "end_enzyme", "fragment_length"])
        w.writerows(all_frags)

    # filtered.csv
    with open(os.path.join(combo_dir, "filtered.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["accession", "start_pos", "start_enzyme",
                    "end_enzyme", "fragment_length"])
        w.writerows(all_filtered)

    # distribution.csv  (standard bins 0-99..>=1000  +  custom window row) 
    all_std_labels = _STD_LABELS + [">=1000"]
    total_all      = sum(total_std.values()) or 1
    custom_label   = f"custom_{size_low}-{size_high}"
    custom_count   = len(all_filtered)

    with open(os.path.join(combo_dir, "distribution.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["length_range", "count", "percentage", "note"])
        for lab in all_std_labels:
            cnt  = total_std.get(lab, 0)
            note = ("selected"
                    if _label_overlaps_range(lab, size_low, size_high) else "")
            w.writerow([lab, cnt, f"{100 * cnt / total_all:.2f}", note])
        w.writerow([custom_label, custom_count,
                    f"{100 * custom_count / total_all:.2f}", "user_window"])

    # gc_metrics.csv  (filtered fragments only)
    gc_values = []
    for acc, sp, se, ee, fl in all_filtered:
        frag_seq = ref_seqs.get(acc, "")[sp : sp + fl]
        if frag_seq:
            gc_values.append(gc_content(frag_seq))

    gc_stats = _compute_gc_stats(gc_values)
    with open(os.path.join(combo_dir, "gc_metrics.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["metric", "value"])
        for k, v in gc_stats.items():
            w.writerow([k, v])

    print(f"    cuts={len(all_cuts):>8,}  fragments={len(all_frags):>8,}  "
          f"filtered[{size_low}-{size_high}]={len(all_filtered):>7,}  "
          f"gc_mean={gc_stats.get('gc_mean_pct', 'n/a')}")

    return {
        "combo":        f"{enz_a_name}_{enz_b_name}",
        "cuts":         len(all_cuts),
        "frags":        len(all_frags),
        "filtered":     len(all_filtered),
        "std_dist":     dict(total_std),
        "gc_stats":     gc_stats,
        "custom_label": custom_label,
    }


def _label_overlaps_range(label, low, high):
    if label == ">=1000":
        return high >= 1000
    try:
        b_lo, b_hi = (int(x) for x in label.split("-"))
        return b_lo <= high and b_hi >= low
    except Exception:
        return False


def _compute_gc_stats(gc_values):
    if not gc_values:
        return {"n_fragments": 0, "gc_mean_pct": "n/a", "gc_median_pct": "n/a",
                "gc_min_pct": "n/a", "gc_max_pct": "n/a", "gc_std_pct": "n/a"}
    n      = len(gc_values)
    mean   = sum(gc_values) / n
    srt    = sorted(gc_values)
    mid    = n // 2
    median = srt[mid] if n % 2 else (srt[mid - 1] + srt[mid]) / 2
    var    = sum((x - mean) ** 2 for x in gc_values) / n
    pct    = lambda v: f"{v * 100:.2f}"
    return {
        "n_fragments":   n,
        "gc_mean_pct":   pct(mean),
        "gc_median_pct": pct(median),
        "gc_min_pct":    pct(min(gc_values)),
        "gc_max_pct":    pct(max(gc_values)),
        "gc_std_pct":    pct(var ** 0.5),
    }


# ANNOTATION  (optional: RepeatMasker .out + GFF/GTF). Applied to filtered fragments only.

def _build_index(records, label_fn):
    raw = defaultdict(list)
    for rec in records:
        raw[rec["chrom"]].append((rec["start"], rec["end"], label_fn(rec)))
    index = {}
    for chrom, intervals in raw.items():
        intervals.sort()
        starts = [iv[0] for iv in intervals]
        ends   = [iv[1] for iv in intervals]
        labels = [iv[2] for iv in intervals]
        index[chrom] = (starts, ends, labels)
    return index


def _query_index(index, chrom, fs, fe):
    if chrom not in index:
        return

    starts, ends, labels = index[chrom]

    hi = _bisect.bisect_left(starts, fe)

    for i in range(hi - 1, -1, -1):

        if ends[i] <= fs:
            break

        ov = min(fe, ends[i]) - max(fs, starts[i])

        if ov > 0:
            yield ov, labels[i]


def parse_te_out(te_path):
    records = []
    with open_file(te_path) as fh:
        for i, line in enumerate(fh):
            if i < 3:
                continue
            parts = line.split()
            if len(parts) < 15:
                continue
            try:
                records.append({
                    "chrom":     parts[4],
                    "start":     int(parts[5]),
                    "end":       int(parts[6]),
                    "te_class":  parts[10].split("/")[0],
                    "te_family": (parts[10].split("/")[1]
                                  if "/" in parts[10] else parts[10]),
                })
            except (ValueError, IndexError):
                continue
    return records


def parse_annotation_gff(ann_path):
    records = []
    with open_file(ann_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip().split("\t")
            if len(parts) < 9:
                continue
            try:
                records.append({
                    "chrom":        parts[0],
                    "start":        int(parts[3]),
                    "end":          int(parts[4]),
                    "feature_type": parts[2].lower(),
                })
            except (ValueError, IndexError):
                continue
    return records

# Base-level overlap between filtered fragments and TE/gene records. Interval indexes are built once per call (O(n log n)), then each fragment query is O(log n + hits) instead of O(n). Writes combo_dir/annotation_summary.csv. Returns dict {category -> (bases, pct_str)} for TSV summary.
def annotate_filtered(filtered_csv, te_records, ann_records, combo_dir):
    frags = []
    with open(filtered_csv) as fh:
        for row in csv.DictReader(fh):
            frags.append({
                "acc":   row["accession"],
                "start": int(row["start_pos"]),
                "end":   int(row["start_pos"]) + int(row["fragment_length"]),
                "len":   int(row["fragment_length"]),
            })

    if not frags:
        return {}

    total_bases = sum(f["len"] for f in frags)

    # Build indexes once - O(n log n)
    ann_index = _build_index(ann_records, lambda r: r["feature_type"])
    te_index  = _build_index(te_records,  lambda r: f"TE_{r['te_class']}_{r['te_family']}")

    ann_counts = defaultdict(int)
    te_counts  = defaultdict(int)

    # Query - O(m * log n) instead of O(m * n)
    for frag in frags:
        fs, fe, acc = frag["start"], frag["end"], frag["acc"]
        for ov, label in _query_index(ann_index, acc, fs, fe):
            ann_counts[label] += ov
        for ov, label in _query_index(te_index, acc, fs, fe):
            te_counts[label] += ov

    ann_result = {}
    out = os.path.join(combo_dir, "annotation_summary.csv")
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["category", "bases", "percentage"])
        w.writerow(["total_filtered_bases", total_bases, "100.00"])
        for feat, cnt in sorted(ann_counts.items()):
            pct = f"{100 * cnt / total_bases:.2f}"
            w.writerow([feat, cnt, pct])
            ann_result[feat] = (cnt, pct)
        for feat, cnt in sorted(te_counts.items()):
            pct = f"{100 * cnt / total_bases:.2f}"
            w.writerow([feat, cnt, pct])
            ann_result[feat] = (cnt, pct)

    print(f"    -> annotation_summary.csv  ({len(ann_result)} categories)")
    return ann_result


# RUN-LEVEL TSV SUMMARY

def write_run_summary_tsv(results, results_dir, size_low, size_high):
    if not results:
        return

    all_ann_cats   = sorted({cat for r in results for cat in r.get("ann_result", {}).keys()})
    all_std_labels = _STD_LABELS + [">=1000"]
    custom_label   = f"custom_{size_low}-{size_high}"

    header = (
        ["combination", "total_cuts", "total_fragments", f"filtered_{custom_label}"]
        + [f"bin_{l}" for l in all_std_labels]
        + [f"bin_{custom_label}"]
        + ["gc_mean_pct", "gc_median_pct", "gc_std_pct",
           "gc_min_pct", "gc_max_pct", "gc_n_fragments"]
        + [f"ann_pct_{c}" for c in all_ann_cats]
    )

    out_path = os.path.join(results_dir, "run_summary.tsv")
    with open(out_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        for r in results:
            gc   = r.get("gc_stats", {})
            dist = r.get("std_dist", {})
            ann  = r.get("ann_result", {})
            w.writerow(
                [r["combo"], r["cuts"], r["frags"], r["filtered"]]
                + [dist.get(l, 0) for l in all_std_labels]
                + [r["filtered"]]
                + [gc.get("gc_mean_pct",   "n/a"),
                   gc.get("gc_median_pct", "n/a"),
                   gc.get("gc_std_pct",    "n/a"),
                   gc.get("gc_min_pct",    "n/a"),
                   gc.get("gc_max_pct",    "n/a"),
                   gc.get("n_fragments",   0)]
                + [ann.get(c, ("n/a", "n/a"))[1] for c in all_ann_cats]
            )

    print(f"  [INFO] run_summary.tsv -> {out_path}")


# ARGUMENT PARSER

def parse_args():
    p = argparse.ArgumentParser(
        prog="rendogbs_pipeline.py",
        description=(
            "rendogbs - In-silico ddRAD/GBS pipeline.\n"
            "Digests a reference genome in silico, selects fragments within a\n"
            "user-defined size window, computes GC metrics, and optionally\n"
            "overlaps fragments with TE and gene annotation.\n\n"
            "All output is written to <workdir>/results/."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Combination modes:
  --combinations fast     25 built-in common ddRAD/GBS pairs (default)
  --combinations all      all ordered pairs from the enzyme dictionary (very slow)
  --combinations custom   read from --combinations-file <path>

Size range:
  --size 200-400          any integer range (both ends inclusive)
                          Standard 100 bp bins (0-99 to 900-999 + >=1000) are
                          always reported; filtered.csv keeps only fragments
                          within [LOW, HIGH].

Examples:
  rendogbs_pipeline.py --ref genome.fa --workdir ./run1 --parallel 11 --size 200-400

  rendogbs_pipeline.py --ref genome.fa.gz --workdir ./run1 \\
      --combinations custom --combinations-file my_pairs.txt \\
      --parallel 8 --size 150-350 \\
      --annotation genes.gff3.gz --te repeats.out
""")

    p.add_argument("--ref",
                   required=True,
                   help="Reference FASTA (.fa / .fasta / .fa.gz)")
    p.add_argument("--workdir",
                   required=True,
                   help="Root working directory; results written to <workdir>/results/")
    p.add_argument("--combinations",
                   default="fast", choices=["fast", "all", "custom"],
                   help="Enzyme combination set: fast (default) | all | custom")
    p.add_argument("--combinations-file",
                   default=None, dest="combinations_file",
                   help="CSV file for custom mode (header + enzyme_a,enzyme_b rows)")
    p.add_argument("--parallel",
                   type=int, default=4,
                   help="Contigs processed in parallel per batch (default: 4)")
    p.add_argument("--size",
                   required=True,
                   help="Fragment size window to retain, e.g. 200-400 (inclusive)")
    p.add_argument("--annotation",
                   default=None,
                   help="[OPTIONAL] GFF3/GFF/GTF gene annotation (plain or .gz)")
    p.add_argument("--te",
                   default=None,
                   help="[OPTIONAL] RepeatMasker .out TE annotation (plain or .gz)")
    return p.parse_args()


# MAIN - zlepsit
def main():
    args = parse_args()

    try:
        size_low, size_high = parse_size_range(args.size)
    except ValueError as e:
        print(f"[ERROR] {e}"); sys.exit(1)

    # All output goes to workdir/results/
    results_dir = os.path.join(args.workdir, "results")
    os.makedirs(results_dir, exist_ok=True)

    # load and sort reference contigs
    print(f"\n[1/4] Loading reference: {args.ref}")
    contigs_sorted = sort_contigs_by_length(read_fasta_contigs(args.ref))
    total_bases    = sum(len(seq) for _, seq in contigs_sorted)
    print(f"      {len(contigs_sorted)} contigs | {total_bases:,} bp total")
    write_contig_lengths_txt(contigs_sorted, results_dir)
    ref_seqs = {acc: seq for acc, seq in contigs_sorted}

    # load enzyme combinations
    print(f"\n[2/4] Loading combinations: {args.combinations}")
    combos = load_combinations(args.combinations, args.combinations_file)
    print(f"      {len(combos)} valid combinations")

    # optional annotation files
    te_records  = []
    ann_records = []
    if args.te:
        print(f"\n[3/4] TE annotation  : {args.te}")
        te_records = parse_te_out(args.te)
        print(f"      {len(te_records):,} TE records")
    if args.annotation:
        print(f"      GFF annotation: {args.annotation}")
        ann_records = parse_annotation_gff(args.annotation)
        print(f"      {len(ann_records):,} annotation records")
    if not args.te and not args.annotation:
        print(f"\n[3/4] Annotation skipped (--te / --annotation not provided)")

    # run pipeline for each combination
    print(f"\n[4/4] Processing {len(combos)} combination(s) | "
          f"parallel={args.parallel} | size={size_low}-{size_high}\n")

    all_results = []

    for idx, (ea, eb) in enumerate(combos, 1):
        combo_name = f"{ea}_{eb}"
        combo_dir  = os.path.join(results_dir, combo_name)
        print(f"  [{idx:>3}/{len(combos)}] {combo_name}")

        result = process_combination(
            ea, eb, contigs_sorted, combo_dir,
            size_low, size_high, args.parallel, ref_seqs
        )

        if te_records or ann_records:
            result["ann_result"] = annotate_filtered(
                os.path.join(combo_dir, "filtered.csv"),
                te_records, ann_records, combo_dir)

        all_results.append(result)

    if all_results:
        write_run_summary_tsv(all_results, results_dir, size_low, size_high)
    else:
        print("  [WARN] No combinations produced results; run_summary.tsv not written.")

    print(f"""



+--------------------------------------------------+
|  rendogbs pipeline - DONE                        |
+--------------------------------------------------+
|  Combinations : {len(combos):<5}                          |
|  Size window  : {size_low}-{size_high:<10}                |
|  Results dir  : {results_dir:<32} |
|  Summary file : results/run_summary.tsv          |
+--------------------------------------------------+
""")


if __name__ == "__main__":
    main()

#upravit tabulku at je hezka, funkci co to vykresli
