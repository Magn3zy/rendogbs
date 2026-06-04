#!/bin/sh

# Handle arguments
INNER_ARGS=
WORKDIR=
REF=
IREF=
while [ -n "$1" ] ; do
    case "$1" in
	--workdir)
	    INNER_ARGS="$INNER_ARGS $1 ./workdir"
	    WORKDIR=${2}
	    shift
	    shift
	    ;;
	--ref)
	    REF=${2}
	    IREF=${2##*/}
	    echo REF=$REF IREF=$IREF
	    INNER_ARGS="$INNER_ARGS $1 $IREF"
	    shift
	    shift
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
       -v $WORKDIR:/home/rendogbs/workdir \
       -v $REF:/home/rendogbs/$IREF \
       --rm \
       rendogbs01 \
       $INNER_ARGS


docker run \
       -v $WORKDIR:/home/rendogbs/workdir \
       -v $REF:/home/rendogbs/$IREF \
       --rm \
       rendogbs01 \
       $INNER_ARGS
