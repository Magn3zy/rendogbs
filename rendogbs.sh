#!/bin/sh

# Inner environment
IHOME=/home/rendogbs
IMAGE=rendogbs-v1

# Nice logging
err() {
    printf "\e[0;31m[ERROR]\e[0m %s\n" "$*"
}

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

# Iterate through all command-line options and their arguments
while [ -n "$1" ] ; do
    case "$1" in
	--workdir)
	    INNER_ARGS="$INNER_ARGS $1 ./workdir"
	    shift
	    if [ -n "$1" ] ; then
		VMAPPINGS="$VMAPPINGS -v $1:$IHOME/workdir"
		ARGSMSK=$((ARGSMSK | 1))
		shift
	    fi
	    ;;
	--ref)
	    INNER_ARGS="$1 ${2##*/}"
	    shift
	    if [ -n "$1" ] ; then
		VMAPPINGS="$VMAPPINGS -v $1:$IHOME/${1##*/}"
		ARGSMSK=$((ARGSMSK | 2))
		shift
	    fi
	    ;;
	--combinations-file)
	    INNER_ARGS="$1 ${2##*/}"
	    shift
	    if [ -n "$1" ] ; then
		VMAPPINGS="$VMAPPINGS -v $1:$IHOME/${1##*/}"
		ARGSMSK=$((ARGSMSK | 16))
		shift
	    fi
	    ;;
	--annotation|--te)
	    INNER_ARGS="$1 ${2##*/}"
	    shift
	    if [ -n "$1" ] ; then
		VMAPPINGS="$VMAPPINGS -v $1:$IHOME/${1##*/}"
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
		INNER_ARGS="$1 $2"
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
    echo "See --help for more information."
    exit 1
fi

# Run the container and pass mappings and arguments to the inner
# wrapper script.
docker run \
       $VMAPPINGS \
       --rm \
       $IMAGE \
       $INNER_ARGS
