#!/bin/sh
#
# Requires only POSIX.1 shell (busybox sh, dash ...)
#
# Starts as root, runs scripts as unprivileged user given by LUID/LGID
# environment variables.
#

# Configuration
PIPELINE=rendogbs_pipeline.py
PLOTS=rendogbs_plots.py

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
# Prints usage instructions - these are shorter than those of the
# outer wrapper script but apart from the container type they are the
# same.
usage() {
    cat <<'EOF'
rendogbs  -  In-silico ddRAD/GBS library pipeline

Runs diges and analysis, then generates plots
(TSV summary + figures) sequentially.  Plots start only after the pipeline
exits successfully.  All output is written to <workdir>/results/.

REQUIRED
  --ref            <file>   Reference FASTA (.fa / .fasta / .fa.gz)
  --workdir        <dir>    Working directory (created if absent)
  --parallel       <int>    Combinations processed in parallel per batch, don't more threads than combinations
  --size           <range>  Fragment size window, e.g. 200-400 (both inclusive).
                            Standard 100 bp bins (0-99 .. 900-999 + >=1000) are
                            always reported; filtered.csv retains only fragments
                            within --size-range.

COMBINATIONS  (optional, default: fast)
  --combinations   fast     23 built-in common ddRAD/GBS pairs (default)
                   custom   read from --combinations-file
  --combinations-file <f>   CSV for custom mode (header line: enzyme_a,enzyme_b),
                            please reffer to --enzymes to choose enzymes from the list 
                            of 171 available enzymes.

ANNOTATION  (optional)
  --annotation     <file>   GFF3/GFF/GTF gene annotation
  --te             <file>   RepeatMasker .out TE annotation

PLOT OPTIONS  (optional)
  --chroms         <int>    Longest N contigs treated as chromosomes in
                            per-chromosome plots (default: 10)
  --dpi            <int>    Figure resolution in DPI (default: 300)
  --skip-plots              Run pipeline only, skip plot generation

OTHER
  -h, --help                Show this help and exit
  --enzymes                 Show available enzymes and exit
  --fast-combinations       Show fast combinations and exit

OUTPUT STRUCTURE
  workdir/
    results/
      contig_lengths.txt          contigs sorted descending by length
      summary.tsv                 Excel-ready: one row per combination,
                                  all bin counts + GC metrics + annotation %
      EcoRI_MseI/
        cuts.csv                  all cut sites (accession, position, enzyme)
        fragments.csv             all adjacent pairs of different enzymes
        filtered.csv              fragments within --size window and with resolved multi-cutter regions
        distribution.csv          standard 100 bp bins + custom window row
        gc_metrics.csv            GC statistics for filtered fragments
        annotation_summary.csv    TE / gene coverage (only with --te / --annotation)
        uncertain_cuts.csv        fragments with ambiguous cut sites
        statistics_uncertain.csv  statistics for ambiguous cut sites
        statistics_cutting.csv    statistics for all cut sites
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
  # Don't use more threads than combinations, pipeline won't finish faster
  # Basic run, fast preset, 200-400 bp window
  rendogbs.sh --ref genome.fasta --workdir ./run1 --parallel 11 --size 200-400

  # All combinations, with annotation, 11 chromosomes in plots
  rendogbs.sh --ref genome.fasta --workdir ./run1 --parallel 11 --size 150-350 \
              --combinations fast \
              --annotation genes.gff3.gz --te repeats.out \
              --chroms 11

  # Custom enzyme pairs, pipeline only
  rendogbs.sh --ref genome.fa --workdir ./run1 --parallel 4 --size 200-500 \
              --combinations custom --combinations-file my_pairs.csv \
              --skip-plots
EOF
}

# May be disabled by --skip-plots
RUN_PLOTS=1

# Local arguments
WORKDIR=

# Passed arguments to both subtasks
S1ARGS=
S2ARGS=

# Iterate over all command-line options and their arguments accordingly
while [ -n "$1" ] ; do
    case "$1" in
	--skip-plots)
	    RUN_PLOTS=0
	    shift
	    ;;
	--chroms|--dpi)
	    PLOTS_ARGS="$PLOTS_ARGS $1 $2"
	    shift
	    shift
	    ;;
	--workdir)
	    WORKDIR="$2"
	    S1ARGS="$S1ARGS $1 $2"
	    PLOTS_ARGS="$PLOTS_ARGS $1 $2"
	    shift
	    shift
	    ;;
	--combinations)
	    S1ARGS="$S1ARGS $1 $2"
	    shift
	    shift
	    ;;
	--combinations-file)
	    S1ARGS="$S1ARGS $1 $2"
	    shift
	    shift
	    ;;
	--size)
	    PLOTS_ARGS="$PLOTS_ARGS $1 $2"
	    shift
	    shift
	    ;;
	--parallel)
	    S2ARGS="$S2ARGS $1 $2"
	    shift
	    shift
	    ;;
	--ref)
	    S2ARGS="$S2ARGS $1 $2"
	    shift
	    shift
	    ;;
	--help|-h)
	    usage
	    exit 0
	    ;;
	*)
	    shift
	    ;;
    esac
done

######## 1) Combinations Processing

info "Step 1: combination_processing.py"
tss1=$(date +%s)
python3 combination_processing.py $S1ARGS

######## 2) Finder (Rust)
tss2=$(date +%s)
info Took $((tss2 - tss1)) seconds.

info "Step 2: rendogbs_finder"
./rendogbs_finder \
  --out-dir "$WORKDIR/results/cuts/" \
  --enzymes-run "$WORKDIR/results/enzymes_run.csv" \
   $S2ARGS

######## 3
tss3=$(date +%s)
info Took $((tss3 - tss2)) seconds.

find "$WORKDIR"

################################################################

tsse=$(date +%s)

info Elapsed time $((tsse - tss1)) seconds.

echo "  Results  : ${WORKDIR}/results/"
