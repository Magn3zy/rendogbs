#!/bin/sh

# Inner environment
IHOME=/home/rendogbs
IMAGE=rendogbs-v1

# Outer environment: current working directory and script directory
wd=`pwd`
sd=`dirname $0`
sd=`readlink -f $sd`

# Nice logging
err() {
    printf "\e[0;31m[ERROR]\e[0m %s\n" "$*"
}

#
# Prints usage instructions
usage() {
    cat <<'EOF'
rendogbs  -  In-silico ddRAD/GBS library pipeline

Runs rendogbs_pipeline.py (digest + analysis), then rendogbs_plots.py
(TSV summary + figures) sequentially.  Plots start only after the pipeline
exits successfully.  All output is written to <workdir>/results/.

CONTAINER TYPE (default: --docker)
  --docker                  Run the docker container (default)
  --singularity             Run the singularity container rendogbs-v1.sif in this
                            directory

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

################################
# Handle arguments

# Generic arguments collected for passing to the inner wrapper
INNER_ARGS=

# Any file or directory arguments generate -v mappings for docker run
VMAPPINGS=

# Mask of (semi-)mandatory arguments
#  1 - workdir
#  2 - ref
#  4 - parallel
#  8 - size
# 16 - combinations-file (mandatory if --combinations custom)
ARGSMSK=0

# Required arguments mask (starts without combinations-file)
REQAMSK=15

# Default container is docker
CCMD=`which docker`
VMOPT=-v
CIDMAP="-e LUID=$(id -u) -e LGID=$(id -g) --rm"
PBSWRAP=
PBSEND=
PBSOPTS=
PBSNAME="-N rendogbs"

# Iterate through all command-line options and their arguments
while [ -n "$1" ] ; do
    case "$1" in
	--singularity)
	    CCMD=`which singularity`
	    IMAGE="$sd/rendogbs-v1.sif"
	    VMOPT=-B
	    CIDMAP=
	    shift
	    ;;
	--docker)
	    CCMD=`which docker`
	    IMAGE=rendogbs-v1
	    VMOPT=v
	    CIDMAP="-e LUID=$(id -u) -e LGID=$(id -g) --rm"
	    shift
	    ;;
	--qsub)
	    PBSWRAP="qsub"
	    PBSEND="--"
	    shift
	    ;;
	--limits|-l)
	    shift
	    PBSOPTS="$PBSOPTS -l $1"
	    shift
	    ;;
	--interactive|-I)
	    PBSOPTS="$PBSOBTS -I"
	    shift
	    ;;
	--name|-N)
	    shift
	    PBSNAME="-N $1"
	    shift
	    ;;
	--workdir)
	    OWORKDIR=`readlink -f "$wd/$2"`
	    INNER_ARGS="$INNER_ARGS $1 ./workdir"
	    shift
	    if [ -n "$1" ] ; then
		VMAPPINGS="$VMAPPINGS $VMOPT $OWORKDIR:$IHOME/workdir"
		ARGSMSK=$((ARGSMSK | 1))
		shift
	    fi
	    ;;
	--ref)
	    INNER_ARGS="$INNER_ARGS $1 ${2##*/}"
	    shift
	    if [ -n "$1" ] ; then
		OREF=`readlink -f "$wd/$1"`
		VMAPPINGS="$VMAPPINGS $VMOPT $OREF:$IHOME/${1##*/}"
		ARGSMSK=$((ARGSMSK | 2))
		shift
	    fi
	    ;;
	--combinations-file)
	    INNER_ARGS="$INNER_ARGS $1 ${2##*/}"
	    shift
	    if [ -n "$1" ] ; then
		OCFILE=`readlink -f "$wd/$1"`
		VMAPPINGS="$VMAPPINGS $VMOPT $OCFILE:$IHOME/${1##*/}"
		ARGSMSK=$((ARGSMSK | 16))
		shift
	    fi
	    ;;
	--annotation|--te)
	    INNER_ARGS="$INNER_ARGS $1 ${2##*/}"
	    shift
	    if [ -n "$1" ] ; then
		OFILE=`readlink -f "$wd/$1"`
		VMAPPINGS="$VMAPPINGS $VMOPT $OFILE:$IHOME/${1##*/}"
		shift
	    fi
	    ;;
	--parallel)
	    if [ -n "$2" ] ; then
		INNER_ARGS="$INNER_ARGS $1 $2"
		ARGSMSK=$((ARGSMSK | 4))
		shift
		shift
	    else
		shift
	    fi
	    ;;
	--size)
	    if [ -n "$2" ] ; then
		INNER_ARGS="$INNER_ARGS $1 $2"
		ARGSMSK=$((ARGSMSK | 8))
		shift
		shift
	    else
		shift
	    fi
	    ;;
	--combinations)
	    if [ -n "$2" ] ; then
		INNER_ARGS="$INNER_ARGS $1 $2"
		if [ "$2" = "custom" ] ; then
		    # Make --combinations-file mandatory for custom
		    # combinations
		    REQAMSK=$((REQAMSK | 16))
		fi
		shift
		shift
	    else
		shift
	    fi
	    ;;
	*)
	    INNER_ARGS="$INNER_ARGS $1"
	    shift
	    ;;
    esac
done

# Report all errors and exit if at least one popped up
if [ $((ARGSMSK & REQAMSK)) -ne $REQAMSK ] ; then
    # Display usage _above_ errors (it is quite long)
    usage
fi
if [ $((ARGSMSK & 1)) -eq 0 ] ; then
    err "Missing required argument: --workdir"
fi
if [ $((ARGSMSK & 2)) -eq 0 ] ; then
    err "Missing required argument: --ref"
fi
if [ $((ARGSMSK & 4)) -eq 0 ] ; then
    err "Missing required argument: --parallel"
fi
if [ $((ARGSMSK & 8)) -eq 0 ] ; then
    err "Missing required argument: --size"
fi
if [ $((ARGSMSK & 16)) -lt $((REQAMSK & 16)) ] ; then
    err "--combinations custom requires --combinations-file <path>"
fi
if [ $((ARGSMSK & REQAMSK)) -ne $REQAMSK ] ; then
    # Although usage is displayed above, this is kinda standard
    echo "See --help for more information."
    exit 1
fi

# Ensure that if someone does not specify --qsub, the PBS options are
# not applied
if [ -z "$PBSWRAP" ] ; then
    PBSEND=
    PBSOPTS=
    PBSNAME=
fi

# Run the container and pass mappings and arguments to the inner
# wrapper script.
echo $PBSWRAP $PBSOPTS $PBSNAME $PBSEND $CCMD run \
      $VMAPPINGS \
      $CIDMAP \
      $IMAGE \
      $INNER_ARGS
echo ================================================================
$PBSWRAP $PBSOPTS $PBSNAME $PBSEND $CCMD run \
      $VMAPPINGS \
      $CIDMAP \
      $IMAGE \
      $INNER_ARGS
