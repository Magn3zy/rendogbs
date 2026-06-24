#!/usr/bin/env bash
# validation_pe.sh
#
# Usage:
#   ./ddrad_recovery.sh <filtered.csv> <bam> <genome.fa>
#
# Optional env vars: (trimming of recognition sites)
#   START_TOL  bp tolerance at fragment start (default 10)
#   END_TOL    bp tolerance at fragment end   (default 10)
#
# filtered.csv columns:
#   accession, start_pos, start_enzyme, end_enzyme, fragment_length
#
# Requirements: samtools, bedtools, python3

set -euo pipefail

REF="${1:?Usage: $0 filtered.csv bam genome.fa}"
BAM="${2:?Usage: $0 filtered.csv bam genome.fa}"
GENOME="${3:?Usage: $0 filtered.csv bam genome.fa}"

START_TOL="${START_TOL:-10}"
END_TOL="${END_TOL:-10}"

OUT_MATCHED="matched_fragments.csv"
OUT_UNMATCHED="unmatched_fragments.csv"
OUT_GC_BINS="gc_bins.tsv"
OUT_COVERAGE="fragment_coverage.tsv"

TMP_REF="$(mktemp --suffix=.bed)"
TMP_READS="$(mktemp --suffix=.bed)"
TMP_OVERLAP="$(mktemp --suffix=.tsv)"
TMP_FASTA="$(mktemp --suffix=.tab)"
TMP_COVERAGE="$(mktemp --suffix=.tsv)"

trap 'rm -f "$TMP_REF" "$TMP_READS" "$TMP_OVERLAP" "$TMP_FASTA" "$TMP_COVERAGE"' EXIT

echo "============================================================"
echo " ddRAD fragment recovery + GC analysis"
echo "  CSV          : $REF"
echo "  BAM          : $BAM"
echo "  Genome       : $GENOME"
echo "  Start tol    : ${START_TOL} bp"
echo "  End tol      : ${END_TOL} bp"
echo "============================================================"

# genome + BAM index
if [[ ! -f "${GENOME}.fai" ]]; then
    echo "[0] Genome index missing — creating .fai ..."
    samtools faidx "$GENOME"
fi
if [[ ! -f "${BAM}.bai" ]]; then
    echo "[0] BAM index missing — creating .bai ..."
    samtools index "$BAM"
fi

# STEP 1 – CSV → BED
#
# BED cols (9):
#   chr  frag_s(0b)  frag_e(0b)  row_id  flen  strand  start_pos  start_enzyme  end_enzyme

echo ""
echo "[1/6] Building reference BED ..."

awk -F',' 'BEGIN{OFS="\t"}
NR==1 { next }
{
    chr  = $1
    spos = $2 + 0
    se   = $3
    ee   = $4
    flen = $5 + 0
    row  = NR - 1

    if (chr=="" || spos<=0 || flen<=0) next

    frag_s = spos - 1
    frag_e = frag_s + flen
    strand = (se=="EcoRI") ? "+" : "-"

    print chr, frag_s, frag_e, row, flen, strand, spos, se, ee
}' "$REF" | LC_ALL=C sort -k1,1 -k2,2n > "$TMP_REF"

TOTAL_REF=$(wc -l < "$TMP_REF")
echo "   Total predicted fragments: $TOTAL_REF"

# STEP 2 – BAM → BED  (-F 2308: drop unmapped/secondary/supplementary)
echo "[2/6] BAM to BED ..."

samtools view -F 2308 "$BAM" | awk 'BEGIN{OFS="\t"}
{
    chr==$3; pos=$4+0; cigar=$6
    chr=$3

    if (chr=="*" || pos<1 || cigar=="*") next

    ref_len=0; t=cigar
    while (match(t, /[0-9]+[MIDNSHP=X]/)) {
        tok=substr(t,RSTART,RLENGTH)
        op=substr(tok,length(tok),1)
        n=substr(tok,1,length(tok)-1)+0
        if (op~/[MDN=X]/) ref_len+=n
        t=substr(t,RSTART+RLENGTH)
    }
    if (ref_len<=0) next

    s=pos-1; e=s+ref_len
    if (s>=0 && e>s) printf "%s\t%d\t%d\n", chr, s, e
}' | LC_ALL=C sort -k1,1 -k2,2n > "$TMP_READS"

TOTAL_READS=$(wc -l < "$TMP_READS")
echo "   Total reads: $TOTAL_READS"

# STEP 3 – bedtools coverage per fragment
#
# counts reads interval is contained within the fragment (with buffer defined)
# by using -f and -F fraction, exact containment so we use our own BED reads and intersect with -f 1.0 on reads
# i.e. 100% of the read must fall within the fragment+tolerance window
#
# expand the fragment BED by tolerance on both sides for coverage counting
# then bedtools coverage -a fragments -b reads counts reads overlapping each
# fragment, containment check is checked in step 4 for matched/unmatched,
# for total read counts per locus we use bedtools coverage directly
echo "[3/6] Computing coverage per fragment (bedtools coverage) ..."

# expand fragment intervals by tolerance for coverage window
awk -v st="$START_TOL" -v et="$END_TOL" 'BEGIN{OFS="\t"}{
    s = $2 - st; if (s<0) s=0
    e = $3 + et
    print $1, s, e, $4, $5, $6
}' "$TMP_REF" | \
bedtools coverage \
    -a stdin \
    -b "$TMP_READS" \
    -counts \
    -sorted > "$TMP_COVERAGE"

# reads_on_target - reads overlapping at least one predicted fragment
# bedtools coverage col7 - number of reads covering interval
READS_ON_TARGET=$(awk '$7 > 0 {sum += $7} END {print sum+0}' "$TMP_COVERAGE")

awk -v r="$READS_ON_TARGET" -v t="$TOTAL_READS" \
    'BEGIN { printf "   Reads on predicted loci : %d / %d (%.2f %%)\n", r, t, (t>0 ? r/t*100 : 0) }'
awk -v r="$READS_ON_TARGET" -v t="$TOTAL_READS" \
    'BEGIN { printf "   Reads off-target        : %d / %d (%.2f %%)\n", t-r, t, (t>0 ? (t-r)/t*100 : 0) }'

# STEP 4 – Candidate overlaps for containment matching
echo "[4/6] Intersecting for fragment matching ..."

bedtools intersect -sorted -wa -wb \
    -a "$TMP_READS" \
    -b "$TMP_REF" > "$TMP_OVERLAP"

echo "   Candidate overlaps: $(wc -l < "$TMP_OVERLAP")"

# STEP 5 – Fragment sequences for GC
echo "[5/6] Extracting fragment sequences ..."

bedtools getfasta \
    -fi "$GENOME" \
    -bed "$TMP_REF" \
    -name \
    -tab > "$TMP_FASTA"

# STEP 6 – Python: containment filter + coverage join + GC bins + outputs
echo "[6/6] GC analysis + writing outputs ..."

python3 - "$REF" "$TMP_OVERLAP" "$TMP_FASTA" "$TMP_COVERAGE" \
          "$OUT_MATCHED" "$OUT_UNMATCHED" \
          "$OUT_GC_BINS" "$OUT_COVERAGE" \
          "$START_TOL" "$END_TOL" \
          "$TOTAL_READS" "$READS_ON_TARGET" << 'PYEOF'
import sys, csv, collections

(ref_path, overlap_path, fasta_path, coverage_path,
 out_m, out_u, out_bins, out_cov,
 start_tol, end_tol,
 total_reads_str, reads_on_target_str) = sys.argv[1:]

start_tol        = int(start_tol)
end_tol          = int(end_tol)
total_reads      = int(total_reads_str)
reads_on_target  = int(reads_on_target_str)

# ── sequences keyed by row_id ────────────────────────────────────────────────
seqs = {}
with open(fasta_path) as fh:
    for line in fh:
        line = line.rstrip("\n")
        if not line:
            continue
        name, seq = line.split("\t", 1)
        row_id = int(name.split("::")[0])
        seqs[row_id] = seq.upper()

# ── coverage per row_id from bedtools coverage ───────────────────────────────
# TMP_COVERAGE cols: chr  s  e  row_id  flen  strand  read_count
coverage = {}
with open(coverage_path) as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 7:
            continue
        row_id    = int(parts[3])
        read_count = int(parts[6])
        coverage[row_id] = read_count

def gc_content(seq):
    valid = [b for b in seq if b in "ACGT"]
    if not valid:
        return None
    return (valid.count("G") + valid.count("C")) / len(valid)

BIN_LABELS = [f"{i*5}-{i*5+5}" for i in range(20)]

def bin_idx(gc_pct):
    return min(int(gc_pct // 5), 19)

# ── containment filter ───────────────────────────────────────────────────────
matched_ids = set()
with open(overlap_path) as fh:
    for line in fh:
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 12:
            continue
        read_s = int(parts[1])
        read_e = int(parts[2])
        frag_s = int(parts[4])
        frag_e = int(parts[5])
        row_id = int(parts[6])
        if read_s >= (frag_s - start_tol) and read_e <= (frag_e + end_tol):
            matched_ids.add(row_id)

# ── read CSV + assign GC / status / coverage ────────────────────────────────
rows_matched     = []
rows_unmatched   = []
matched_counts   = collections.Counter()
unmatched_counts = collections.Counter()

with open(ref_path, newline="") as fh:
    reader = csv.reader(fh)
    next(reader)
    for row_id, row in enumerate(reader, start=1):
        seq      = seqs.get(row_id, "")
        gc       = gc_content(seq)
        gc_str   = f"{gc*100:.2f}" if gc is not None else "NA"
        status   = "matched" if row_id in matched_ids else "unmatched"
        n_reads  = coverage.get(row_id, 0)
        out_row  = row + [gc_str, str(n_reads), status]

        if row_id in matched_ids:
            rows_matched.append(out_row)
            if gc is not None:
                matched_counts[bin_idx(gc * 100)] += 1
        else:
            rows_unmatched.append(out_row)
            if gc is not None:
                unmatched_counts[bin_idx(gc * 100)] += 1

total_frags = len(rows_matched) + len(rows_unmatched)

# ── summary ──────────────────────────────────────────────────────────────────
print()
print("  ── Fragment recovery ──────────────────────────────────")
print(f"   Predicted fragments : {total_frags}")
print(f"   Matched             : {len(rows_matched)}  ({len(rows_matched)/total_frags*100:.2f} %)")
print(f"   Unmatched           : {len(rows_unmatched)}  ({len(rows_unmatched)/total_frags*100:.2f} %)")

print()
print("  ── Read assignment (bedtools coverage) ────────────────")
print(f"   Total reads (BAM)   : {total_reads}")
print(f"   Reads on predicted  : {reads_on_target}  ({reads_on_target/total_reads*100:.2f} % of all reads)")
print(f"   Reads off-target    : {total_reads - reads_on_target}  ({(total_reads-reads_on_target)/total_reads*100:.2f} % of all reads)")
print()

# ── write CSVs ───────────────────────────────────────────────────────────────
OUT_HEADER = ["accession","start_pos","start_enzyme","end_enzyme",
              "fragment_length","gc_content","reads_in_fragment","status"]

for path, rows in [(out_m, rows_matched), (out_u, rows_unmatched)]:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(OUT_HEADER)
        w.writerows(rows)

# ── coverage TSV (all fragments, sorted by reads desc) ───────────────────────
with open(out_cov, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(OUT_HEADER)
    all_rows = sorted(rows_matched + rows_unmatched,
                      key=lambda r: int(r[6]), reverse=True)
    w.writerows(all_rows)

# ── GC bins TSV ──────────────────────────────────────────────────────────────
total_m = sum(matched_counts.values())
total_u = sum(unmatched_counts.values())

with open(out_bins, "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["gc_bin_pct",
                "matched_n",   "matched_pct_of_all_matched",
                "unmatched_n", "unmatched_pct_of_all_unmatched",
                "total_n",     "recovery_rate_in_bin"])
    for i, label in enumerate(BIN_LABELS):
        mc  = matched_counts.get(i, 0)
        uc  = unmatched_counts.get(i, 0)
        tot = mc + uc
        w.writerow([
            label,
            mc,   f"{mc/total_m*100:.2f}" if total_m else "0.00",
            uc,   f"{uc/total_u*100:.2f}" if total_u else "0.00",
            tot,  f"{mc/tot*100:.2f}" if tot else "0.00",
        ])

print(f"   {out_m}: {len(rows_matched)} rows")
print(f"   {out_u}: {len(rows_unmatched)} rows")
print(f"   {out_bins}: written")
print(f"   {out_cov}: written")
PYEOF

echo ""
echo "============================================================"
echo " Done!"
echo "  matched_fragments.csv   – fragments supported by reads (+ GC + read count)"
echo "  unmatched_fragments.csv – fragments not supported       (+ GC + read count)"
echo "  fragment_coverage.tsv   – all fragments sorted by read count"
echo "  gc_bins.tsv             – 5 % GC bin recovery table"
echo "============================================================"