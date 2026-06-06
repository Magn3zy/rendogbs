#!/bin/sh
#
# Docker image entrypoint for rendogbs. Handles permissions and
# running as appropriate user.
#
# Copyright (c) 2026 Dominik Pantůček ORCID 0009-0000-3509-0905
#
# Requirements:
#  - POSIX.1 shell
#

# Get to the script directory
cd ${0%/*}

#
# Logging functions
info() { printf "\e[0;36m[INFO]\e[0m  %s\n" "$*" }
warn() { printf "\e[1;33m[WARN]\e[0m  %s\n" "$*" }

# Setup permissions
if [ -n "$LUID" ] ; then
    chown $LUID -R .
fi
if [ -n "$LGID" ] ; then
    chgrp $LGID -R .
fi
if [ -n "$LUID" -a -n "$LGID" ] ; then
    info "Running as $LUID:$LGID"
else if [ -n "$LUID" ] ; then
	 warn "Running as $LUID:0 (local group not specified!)"
     else if [ -n "$LGID" ] ; then
	      warn "Running as 0:$LGID (local user not specified!)"
	  else
	      warn "Running as 0:0 (neither local user nor group specified!)"
	  fi
     fi
fi

