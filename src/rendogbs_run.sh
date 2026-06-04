#!/bin/sh
#
# Requires only POSIX.1 shell (busybox sh, dash ...)
#

# Configuration
PIPELINE=rendogbs_pipeline.py
PLOTS=rendogbs_plots.py

# Get to the script directory
cd ${0%/*}

#
# Logging functions
info() {
    printf "\e[0;36m[INFO]\e[0m  %s\n" "$*"
}
err() {
    printf "\e[0;31m[ERROR]\e[0m %s\n" "$*"
}
warn() {
    printf "\e[1;33m[WARN]\e[0m  %s\n" "$*"
}
ok() {
    printf "\e[0;32m[OK]\e[0m    %s\n" "$*"
}

#
# Prints usage instructions
usage() {
    cat <<'EOF'
rendogbs  -  In-silico ddRAD/GBS library pipeline

Runs rendogbs_pipeline.py (digest + analysis), then rendogbs_plots.py
(TSV summary + figures) sequentially.  Plots start only after the pipeline
exits successfully.  All output is written to <workdir>/results/.

REQUIRED
  --ref            <file>   Reference FASTA (.fa / .fasta / .fa.gz)
  --workdir        <dir>    Working directory (created if absent)
  --parallel       <int>    Contigs processed in parallel per batch
  --size           <range>  Fragment size window, e.g. 200-400 (both inclusive).
                            Standard 100 bp bins (0-99 .. 900-999 + >=1000) are
                            always reported; filtered.csv retains only fragments
                            within [LOW, HIGH].

COMBINATIONS  (optional, default: fast)
  --combinations   fast     25 built-in common ddRAD/GBS pairs (default)
                   all      all ordered enzyme pairs from the dictionary
                   custom   read from --combinations-file
  --combinations-file <f>   CSV for custom mode (header line + enzyme_a,enzyme_b rows)

ANNOTATION  (optional)
  --annotation     <file>   GFF3/GFF/GTF gene annotation (plain or .gz)
  --te             <file>   RepeatMasker .out TE annotation (plain or .gz)

PLOT OPTIONS  (optional)
  --chroms         <int>    Longest N contigs treated as chromosomes in
                            per-chromosome plots (default: 10)
  --dpi            <int>    Figure resolution in DPI (default: 300)
  --skip-plots              Run pipeline only, skip plot generation

OTHER
  -h, --help                Show this help and exit

OUTPUT STRUCTURE
  workdir/
    results/
      contig_lengths.txt          contigs sorted descending by length
      run_summary.tsv             Excel-ready: one row per combination,
                                  all bin counts + GC metrics + annotation %
      EcoRI_MseI/
        cuts.csv                  all cut sites (accession, position, enzyme)
        fragments.csv             all adjacent pairs of different enzymes
        filtered.csv              fragments within --size window
        distribution.csv          standard 100 bp bins + custom window row
        gc_metrics.csv            GC statistics for filtered fragments
        annotation_summary.csv    TE / gene coverage (only with --te / --annotation)
      AciI_HindIII/
        ...
      plots/
        heatmap_fragment_lengths.png    fragment count heatmap (10 bp bins, log scale)
        heatmap_chrom_distribution.png  filtered fragments per chromosome (heatmap)
        bar_chrom_distribution.png      filtered fragments per chromosome (bar chart)
        gc_distribution.png             GC% per combination: mean +/- SD
        annotation_coverage.png         annotation category coverage stacked bar
        size_distributions.png          100 bp bin line plot, user window highlighted

EXAMPLES
  # Basic run, fast preset, 200-400 bp window
  rendogbs.sh --ref genome.fa --workdir ./run1 --parallel 11 --size 200-400

  # All combinations, with annotation, 24 chromosomes in plots
  rendogbs.sh --ref genome.fa.gz --workdir ./run1 --parallel 8 --size 150-350 \
              --combinations all \
              --annotation genes.gff3.gz --te repeats.out \
              --chroms 24

  # Custom enzyme pairs, pipeline only
  rendogbs.sh --ref genome.fa --workdir ./run1 --parallel 4 --size 200-500 \
              --combinations custom --combinations-file my_pairs.txt \
              --skip-plots

EOF
}

# Print help if no arguments given

if [ -z "$1" ] ; then
    usage
    exit 1
fi

# Parse command-line arguments, mandatory args are checked in
# parse_args of rendogbs_pipeline.py
SKIP_PLOTS=0
PIPELINE_ARGS=
PLOTS_EXTRA=
WORKDIR=
SIZE=
refok=0
workdirok=0
parallelok=0
sizeok=0
while [ -n "$1" ] ; do
    case "$1" in
	--skip-plots)
	    SKIP_PLOTS=1
	    shift
	    ;;
	--chroms|--dpi)
	    PLOTS_EXTRA="$PLOTS_EXTRA $2 $3"
	    shift
	    shift
	    shift
	    ;;
	--workdir)
	    WORKDIR="$2"
	    PIPELINE_ARGS="$PIPELINE_ARGS $1 $2"
	    workdirok=1
	    shift 2
	    ;;
	--size)
	    SIZE="$2"
	    PIPELINE_ARGS="$PIPELINE_ARGS $1 $2"
	    sizeok=1
	    shift
	    shift
	    ;;
	--parallel)
	    PIPELINE_ARGS="$PIPELINE_ARGS $1"
	    parallelok=1
	    shift
	    ;;
	--ref)
	    PIPELINE_ARGS="$PIPELINE_ARGS $1"
	    refok=1
	    shift
	    ;;
	--help|-h)
	    usage
	    exit 0
	    ;;
	*)
	    PIPELINE_ARGS="$PIPELINE_ARGS $1"
	    shift
	    ;;
    esac
done

# Check mandatory arguments beforehand
for req in ref workdir parallel size ; do
    eval "ok=\$${req}ok"
    if [ $ok -eq 0 ] ; then
        err "Missing required argument: --${req%ok}"
        echo "See --help for more information."
        exit 1
    fi
done

# Run the pipeline
info "Step 1/2 — Running pipeline ..."
pipeline_start=$(date +%s)
python3 $PIPELINE ${PIPELINE_ARGS}
pipeline_exit=$?
pipeline_elapsed=$(( $(date +%s) - pipeline_start ))
echo ""
if [ $pipeline_exit -ne 0 ]; then
    err "Pipeline failed (exit code $pipeline_exit). Plots will not run."
    exit $pipeline_exit
fi
ok "Pipeline finished in ${pipeline_elapsed}s."
echo ""

plots_elapsed=0

if [ $SKIP_PLOTS -eq 1 ]; then
    warn "Plots skipped (--skip-plots)."
else
    info "Step 2/2 — Generating plots and rebuilding run_summary.tsv ..."
    echo ""

    plots_start=$(date +%s)
    echo python3 "$PLOTS" \
            --workdir "$WORKDIR" \
            --size    "$SIZE"    \
            "${PLOTS_EXTRA}"
    python3 "$PLOTS" \
            --workdir "$WORKDIR" \
            --size    "$SIZE"    \
            ${PLOTS_EXTRA}
    plots_exit=$?
    plots_elapsed=$(( $(date +%s) - plots_start ))

    echo ""
    if [ $plots_exit -ne 0 ]; then
        err "Plot generation failed (exit code $plots_exit)."
        exit $plots_exit
    fi
    ok "Plots finished in ${plots_elapsed}s."
fi

total_elapsed=$(( $(date +%s) - pipeline_start ))
echo ""
ok "All done."
echo ""
echo "  Pipeline : ${pipeline_elapsed}s"
[ $SKIP_PLOTS -eq 0 ] && echo "  Plots    : ${plots_elapsed}s"
echo "  Total    : ${total_elapsed}s"
echo "  Results  : ${WORKDIR}/results/"
echo ""
