#!/bin/sh

INITWAIT=10
UNIQNAME=`date +%Y%m%dT%H%M%S`-`mktemp -u XXXXXXXXXXXXXXXX`
WORKDIR=.
CTYPE=docker

# Peek into arguments and extract useful ones
peek_args() {
    while [ -n "$1" ] ; do
	case $1 in
	    --workdir)
		WORKDIR=$2
		shift
		shift
		;;
	    *)
		shift
		;;
	esac
    done
}
peek_args "$@"
echo Output directory: $WORKDIR

# Measure "real" time
rtstart=$(date +%s)

# Start the pipeline as background job which runs the docker container
sh rendogbs.sh "$@" --cname $UNIQNAME >$WORKDIR/$UNIQNAME.output 2>&1 &

# Keep the pipeline PID for cleanup upon exit
PIPID=$!
echo pipeline PID is $PIPID

# Now wait for docker to start and get the container long id
echo Waiting for docker container $UNIQNAME to start ...
LONGID=
attempt=0
SCOPE=
while [ -z "$SCOPE" -a $attempt -le $INITWAIT ] ; do
    if [ $attempt -gt 0 ] ; then
	echo Sleeping 1s ...
	sleep 1
    fi
    attempt=$((attempt + 1))
    LONGID=$(docker ps --no-trunc -f name=$UNIQNAME --quiet)
    if [ -n "$LONGID" ] ; then
	echo Docker container running as $LONGID
	SCOPE=/sys/fs/cgroup/system.slice/docker-$LONGID.scope
    fi
done

# Check
if [ -z "$SCOPE" ] ; then
    echo Cannot get cgroup scope.
    echo Cannot get running Docker container long id, bailing out.
    kill $PIPID
    wait $PIPID
    exit 1
fi

echo Using scope: $SCOPE

# Parse key-value as list of arguments
parse_cpustat() {
    while ! [ -z "$2" ] ; do
	case $1 in
	    usage_usec)
		usage_usec=$2
		;;
	    user_usec)
		user_usec=$2
		;;
	    system_usec)
		system_usec=$2
		;;
	esac
	shift
	shift
    done
}

# Periodically print the stats
(
    echo "#TS THREADS MAXTHREADS MEMORY MAXMEMORY USER SYSTEM TOTAL"
    maxthreads=0
    while ps $PIPID >/dev/null 2>&1 ; do
	ts=$(date +%s)
	nthreadss=$(wc -l $SCOPE/cgroup.threads 2>/dev/null)
	nthreads=${nthreadss%% *}
	if [ -n "$nthreads" ] ; then
	    if [ $maxthreads -lt $nthreads ] ; then
		maxthreads=$nthreads
	    fi
	fi
	mempeak=$(cat $SCOPE/memory.peak 2>/dev/null)
	memcur=$(cat $SCOPE/memory.current 2>/dev/null)
	parse_cpustat $(cat $SCOPE/cpu.stat 2>/dev/null)
	if [ -n "$nthreads" -a -n "$memcur" -a -n "$mempeak" -a -n "$user_usec" -a -n "$system_usec" -a -n "$usage_usec" ] ; then
	    echo $ts $nthreads $maxthreads $memcur $mempeak $user_usec $system_usec $usage_usec
	fi
	sleep 1
    done
) | tee $WORKDIR/$UNIQNAME.data

# Generate GNUPlot source and run it
cat <<EOF >$WORKDIR/$UNIQNAME.gnuplot
name='$UNIQNAME'
data=name . '.data'
outtmname=name . '-tm.svg'
outcname=name . '-c.svg'

set terminal svg size 960,600 font "Sans,16"
set output outtmname

set grid
set object 1 rectangle from screen 0,0 to screen 1,1 fillcolor rgb "white" behind
set title 'Threads and Memory Usage'
set xdata time
set timefmt '%s'
set y2tics
set ylabel "Threads"
set y2label "Memory [MB]"
set xlabel "Running Time"

plot data u 1:2 w l t 'Threads', \
     data u 1:3 w l t 'Max Threads', \
     data u 1:(\$4/1024/1024) w l t 'Memory' axes x1y2, \
     data u 1:(\$5/1024/1024) w l t 'Max Memory' axes x1y2

set output outcname

set title 'CPU Time Usage'

unset y2label
unset y2tics
set ylabel "CPU Time [s]"
plot data u 1:(\$6/1000000) w l t "User", \
     data u 1:(\$7/1000000) w l t "System", \
     data u 1:(\$8/1000000) w l t "Total"
EOF
owd=$(pwd)
cd $WORKDIR
gnuplot $UNIQNAME.gnuplot
cd "$owd"

# Done
rtend=$(date +%s)
wait $PIPID
echo DONE

# Summary
echo
(
    echo Ran with options: "$@"
    echo
    echo Maximum threads: $(grep . $WORKDIR/$UNIQNAME.data|tail -n 1|awk '{print $3}')
    mempeak=$(grep . $WORKDIR/$UNIQNAME.data|tail -n 1|awk '{print $5}')
    echo Memory peak usage: $((mempeak / 1024 / 1024)) MB
    user_usec=$(grep . $WORKDIR/$UNIQNAME.data|tail -n 1|awk '{print $6}')
    system_usec=$(grep . $WORKDIR/$UNIQNAME.data|tail -n 1|awk '{print $7}')
    total_usec=$(grep . $WORKDIR/$UNIQNAME.data|tail -n 1|awk '{print $8}')
    echo User time: $((user_usec / 1000000)) s
    echo System time: $((system_usec / 1000000)) s
    echo Total time: $((total_usec / 1000000)) s
    echo
    echo Real time: $((rtend - rtstart)) s
) | tee $WORKDIR/$UNIQNAME.result
echo
