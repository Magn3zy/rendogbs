#!/bin/sh

# Inner environment
IHOME=/home/rendogbs

# Handle arguments
INNER_ARGS=
WORKDIR=
VMAPPINGS=
while [ -n "$1" ] ; do
    case "$1" in
	--workdir)
	    INNER_ARGS="$INNER_ARGS $1 ./workdir"
	    WORKDIR=${2}
	    VMAPPINGS="$VMAPPINGS -v $WORKDIR:$IHOME/workdir"
	    shift
	    shift
	    ;;
	--ref)
	    REF=${2}
	    VMAPPINGS="$VMAPPINGS -v $REF:$IHOME/${2##*/}"
	    INNER_ARGS="$1 ${2##*/}"
	    shift
	    shift
	    ;;
	--combinations-file)
	    COMBINATIONS_FILE=${2}
	    ICOMBINATIONS_FILE=${2##*/}
	    INNER_ARGS="$1 ${##*/}"
	    ;;
	*)
	    INNER_ARGS="$INNER_ARGS $1"
	    shift
	    ;;
    esac
done

if [ -z "$WORKDIR" ] ; then
    echo you must specify workdir
fi


echo docker run \
       $VMAPPINGS \
       --rm \
       rendogbs01 \
       $INNER_ARGS


docker run \
       $VMAPPINGS \
       --rm \
       rendogbs01 \
       $INNER_ARGS
