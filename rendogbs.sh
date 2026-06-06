#!/bin/sh

# Inner environment
IHOME=/home/rendogbs
IMAGE=rendogbs-v1

# Handle arguments
INNER_ARGS=
VMAPPINGS=
ARGSMSK=0
while [ -n "$1" ] ; do
    case "$1" in
	--workdir)
	    INNER_ARGS="$INNER_ARGS $1 ./workdir"
	    VMAPPINGS="$VMAPPINGS -v $2:$IHOME/workdir"
	    ARGSMSK=$((ARGSMSK | 1))
	    shift
	    shift
	    ;;
	--ref)
	    VMAPPINGS="$VMAPPINGS -v $2:$IHOME/${2##*/}"
	    INNER_ARGS="$1 ${2##*/}"
	    ARGSMSK=$((ARGSMSK | 2))
	    shift
	    shift
	    ;;
	--combinations-file|--annotation|--te)
	    VMAPPINGS="$VMAPPINGS -v $2:$IHOME/${2##*/}"
	    INNER_ARGS="$1 ${##*/}"
	    shift
	    shift
	    ;;
	--parallel)
	    INNER_ARGS="$INNER_ARGS $1"
	    ARGSMSK=$((ARGSMSK | 4))
	    shift
	    ;;
	--size)
	    INNER_ARGS="$INNER_ARGS $1"
	    ARGSMSK=$((ARGSMSK | 8))
	    shift
	    ;;
	*)
	    INNER_ARGS="$INNER_ARGS $1"
	    shift
	    ;;
    esac
done

# Report all errors and exit if at least one popped up
if [ $((ARGSMSK & 1)) -eq 0 ] ; then
    echo you must specify workdir
fi
if [ $((ARGSMSK & 2)) -eq 0 ] ; then
    echo you must specify ref
fi
if [ $((ARGSMSK & 4)) -eq 0 ] ; then
    echo you must specify parallel
fi
if [ $((ARGSMSK & 8)) -eq 0 ] ; then
    echo you must specify size
fi
if [ $ARGSMSK -ne 15 ] ; then
    exit 1
fi

# Run the container and pass mappings and arguments to the inner
# wrapper script.
echo docker run \
       $VMAPPINGS \
       --rm \
       $IMAGE \
       $INNER_ARGS

docker run \
       $VMAPPINGS \
       --rm \
       $IMAGE \
       $INNER_ARGS
