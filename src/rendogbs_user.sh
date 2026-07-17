#!/bin/sh
#
# rendogbs_user.sh
# Copyright (c) 2026 Dominik Pantůček ORCID 0009-0000-3509-0905
#
# Docker image entrypoint for rendogbs. Handles permissions and
# running as appropriate user.
#
# Requirements:
#  - POSIX.1 shell
#

# Get to the script directory
cd ${0%/*}

#
# Logging functions
info() { printf "\e[0;36m[INFO]\e[0m  %s\n" "$*" ; }
warn() { printf "\e[1;33m[WARN]\e[0m  %s\n" "$*" ; }

# Setup permissions
if [ -n "$LUID" ] ; then
    chown $LUID -R .
    if [ -n "$LGID" ] ; then
	chgrp $LGID -R .
	info "Running as $LUID:$LGID"
    else
	warn "Running as $LUID:0 (local group not specified!)"
    fi
    # Prepare environment
    DUID=$LUID
    DGID=${LGID:-0}
    if [ "$DGID" -eq 0 ] ; then
        if ! [ "$DUID" -eq 0 ] ; then
	    adduser -G root -D -H -u $DUID rendogbs
	fi
    else
	GRP=$(getent group $DGID | cut -d: -f1)
	if [ -z "$GRP" ] ; then
	    GRP=rendogbs
	    addgroup -g $DGID rendogbs
	fi
	adduser -G $GRP -D -H -u $DUID rendogbs
    fi
    if [ $DUID -eq 0 ] ; then
        sh rendogbs_run.sh $*
    else
        su rendogbs -c "sh rendogbs_run.sh $*"
    fi
else
    # Run as root (or current user - under singularity)
    if [ -z "$SINGULARITY_CONTAINER" ] ; then
	# If $SINGULARITY_CONTAINER is set, this is not true
	warn "Running as 0:0 (local user not specified!)"
    fi
    /bin/sh rendogbs_run.sh "$@"
fi
