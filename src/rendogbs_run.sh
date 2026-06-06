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

# Parse command-line arguments, mandatory arguments to the python
# scripts were checked in the outer wrapper.

# Inner script arguments are consumed other arguments are passed on to
# pipeline or plots (or both) subtasks as needed.

# May be disabled by --skip-plots
RUN_PLOTS=1

# Passed arguments to both subtasks
PIPELINE_ARGS=
PLOTS_ARGS=

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
	    PIPELINE_ARGS="$PIPELINE_ARGS $1 $2"
	    PLOTS_ARGS="$PLOTS_ARGS $1 $2"
	    shift
	    shift
	    ;;
	--size)
	    PIPELINE_ARGS="$PIPELINE_ARGS $1 $2"
	    PLOTS_ARGS="$PLOTS_ARGS $1 $2"
	    shift
	    shift
	    ;;
	--parallel)
	    PIPELINE_ARGS="$PIPELINE_ARGS $1"
	    shift
	    ;;
	--ref)
	    PIPELINE_ARGS="$PIPELINE_ARGS $1"
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

# Run the pipeline
info "Step 1/2 - Running pipeline ..."
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
    info "Step 2/2 - Generating plots and rebuilding run_summary.tsv ..."
    echo ""

    plots_start=$(date +%s)
    echo python3 "$PLOTS" \
            "${PLOTS_ARGS}"
    python3 "$PLOTS" \
            ${PLOTS_ARGS}
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
