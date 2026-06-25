#!/bin/sh
#
# rendogbs_user.sh
#
# Copyright (c) 2026
# Eliška Korbová ORCID 0009-0004-1247-0808
# Dominik Pantůček ORCID 0009-0000-3509-0905
#
# Requires only POSIX.1 shell (busybox sh, dash ...)
#
# Runs the script under current user. See rendogbs_user.sh for
# handling root/non-root setup in docker and singularity/apptainer
# containers.
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
DO_ANNOTATION=0

# Arguments passed to different subtasks
S1ARGS=
S2ARGS=
S3ARGS=
S4ARGS=
S5ARGS=
S6ARGS=
S7ARGS=
S8ARGS=
S9ARGS=

# Iterate over all command-line options and their arguments accordingly
while [ -n "$1" ] ; do
    case "$1" in
	--skip-plots)
	    RUN_PLOTS=0
	    shift
	    ;;
	--chroms|--dpi)
	    S9ARGS="$S9ARGS $1 $2"
	    shift
	    shift
	    ;;
	--workdir)
	    WORKDIR="$2"
	    S1ARGS="$S1ARGS $1 $2"
	    S4ARGS="$S4ARGS $1 $2"
	    S6ARGS="$S6ARGS $1 $2"
	    S7ARGS="$S7ARGS $1 $2"
	    S8ARGS="$S8ARGS $1 $2"
	    S9ARGS="$S9ARGS $1 $2"
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
	    S4ARGS="$S4ARGS $1 $2"
	    S6ARGS="$S6ARGS $1 $2"
	    S9ARGS="$S9ARGS $1 $2"
	    shift
	    shift
	    ;;
	--parallel)
	    S2ARGS="$S2ARGS $1 $2"
	    S3ARGS="$S3ARGS $1 $2"
	    S4ARGS="$S4ARGS $1 $2"
	    S5ARGS="$S5ARGS $1 $2"
	    S6ARGS="$S6ARGS $1 $2"
	    S7ARGS="$S7ARGS $1 $2"
	    shift
	    shift
	    ;;
	--ref)
	    S2ARGS="$S2ARGS $1 $2"
	    S6ARGS="$S6ARGS $1 $2"
	    shift
	    shift
	    ;;
	--help|-h)
	    usage
	    exit 0
	    ;;
	--annotation|--te)
	    DO_ANNOTATION=1
	    S7ARGS="$S7ARGS $1 $2"
	    shift
	    shift
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
res=$?
if [ $res -ne 0 ] ; then
    err Step 1 failed: $res
    exit 1
fi

######## 2) Finder (Rust)
tss2=$(date +%s)
info Took $((tss2 - tss1)) seconds.

info "Step 2: rendogbs_finder"
./rendogbs_finder \
  --out-dir "$WORKDIR/results/cuts/" \
  --enzymes-run "$WORKDIR/results/enzymes_run.csv" \
   $S2ARGS
res=$?
if [ $res -ne 0 ] ; then
    err Step 2 failed: $res
    exit 1
fi

######## 3) Cut Merge (Rust)
tss3=$(date +%s)
info Took $((tss3 - tss2)) seconds.

info "Step 3: cut_merge"
./cut_merge \
  --out-dir "$WORKDIR/results/" \
  --combinations "$WORKDIR/results/combinations.csv" \
  --cuts-dir "$WORKDIR/results/cuts" \
  $S3ARGS
res=$?
if [ $res -ne 0 ] ; then
    err Step 3 failed: $res
    exit 1
fi

######## 4) Fragment Generation
tss4=$(date +%s)
info Took $((tss4 - tss3)) seconds.

info "Step 4: fragment_generation.py"
python3 fragment_generation.py $S4ARGS
res=$?
if [ $res -ne 0 ] ; then
    err Step 4 failed: $res
    exit 1
fi

######## 5) Fragment Generation (Rust)
tss5=$(date +%s)
info Took $((tss5 - tss4)) seconds.

info "Step 5: fragment_generation"
./fragment_generation \
  --out-dir "$WORKDIR/results/" \
  --combinations "$WORKDIR/results/combinations.csv" \
  $S5ARGS
res=$?
if [ $res -ne 0 ] ; then
    err Step 5 failed: $res
    exit 1
fi

######## 6) Postprocess Metrics
tss6=$(date +%s)
info Took $((tss6 - tss5)) seconds.

info "Step 6: postprocess_metrics.py"
python3 postprocess_metrics.py $S6ARGS
res=$?
if [ $res -ne 0 ] ; then
    err Step 6 failed: $res
    exit 1
fi

######## 7) Annotation
tss7=$(date +%s)
info Took $((tss7 - tss6)) seconds.

if [ $DO_ANNOTATION -eq 1 ] ; then
    info "Step 7: annotation.py"
    python3 annotation.py $S7ARGS
    res=$?
    if [ $res -ne 0 ] ; then
	err Step 7 failed: $res
	exit 1
    fi
else
    info "Step 7: annotation.py skipped (no --te / --annotation provided)"
fi

######## 8) Summary
tss8=$(date +%s)
info Took $((tss8 - tss7)) seconds.

info "Step 8: summary.py"
python3 summary.py $S8ARGS
res=$?
if [ $res -ne 0 ] ; then
    err Step 8 failed: $res
    exit 1
fi


######## 9) Plots
tss9=$(date +%s)
info Took $((tss9 - tss8)) seconds.

if [ $RUN_PLOTS -eq 1 ] ; then
    info "Step 9: plots.py"
    python3 plots.py $S9ARGS
    res=$?
    if [ $res -ne 0 ] ; then
	err Step 9 failed: $res
	exit 1
    fi
else
    info "Step 9: plots.py skipped (--skip-plots)"
fi

################################################################

tsse=$(date +%s)

summary_table() {
    echo " Step |  Begin  |   End   | Duration"
    echo "------+---------+---------+----------"
    idx=1
    zero=$1
    while [ -n "$2" ] ; do
	start=$(($1 - zero))
	end=$(($2 - zero))
	dur=$((end - start))
	printf "   %2d | %7d | %7d |  %7d\n" $idx $start $end $dur
	shift
	idx=$((idx + 1))
    done
}

echo
summary_table $tss1 $tss2 $tss3 $tss4 $tss5 $tss6 $tss7 $tss8 $tss9 $tsse
echo

echo "  Results  : ${WORKDIR}/results/"
