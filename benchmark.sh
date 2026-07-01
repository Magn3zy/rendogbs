#!/bin/sh

INITWAIT=10
UNIQNAME=`mktemp -u XXXXXXXXXXXXXXXX`

# Start the pipeline as background job which runs the docker container
sh rendogbs.sh \
   --ref ./data/GCF_000001735.3_TAIR10_genomic.fna.gz \
   --workdir ./runbench \
   --parallel 2 \
   --size 200-400 \
   --annotation data/GCF_000001735.3_TAIR10_genomic.gff.gz \
   --chroms 5 \
   --te ./data/GCF_000001735.3_TAIR10_rm.out.gz \
   --cname $UNIQNAME >/dev/null 2>/dev/null &

# Keep the pipeline PID for cleanup upon exit
PIPID=$!
echo pipeline PID is $PIPID

# Now wait for docker to start and get the container long id
echo Waiting for docker container $UNIQNAME to start ...
LONGID=
attempt=0
while [ -z "$LONGID" -a $attempt -le $INITWAIT ] ; do
    if [ $attempt -gt 0 ] ; then
	echo Waiting 1s ...
	sleep 1
    fi
    attempt=$((attempt + 1))
    LONGID=$(docker ps --no-trunc -f name=$UNIQNAME --quiet)
done

# Check
if [ -z "$LONGID" ] ; then
    echo Cannot get running Docker container long id, bailing out.
    kill $PIPID
    wait $PIPID
    exit 1
fi
echo Docker container running as $LONGID

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
	nthreadss=$(wc -l /sys/fs/cgroup/system.slice/docker-$LONGID.scope/cgroup.threads 2>/dev/null)
	nthreads=${nthreadss%% *}
	if [ $maxthreads -lt $nthreads ] ; then
	    maxthreads=$nthreads
	fi
	mempeak=$(cat /sys/fs/cgroup/system.slice/docker-$LONGID.scope/memory.peak 2>/dev/null)
	memcur=$(cat /sys/fs/cgroup/system.slice/docker-$LONGID.scope/memory.current 2>/dev/null)
	parse_cpustat $(cat /sys/fs/cgroup/system.slice/docker-$LONGID.scope/cpu.stat 2>/dev/null)
	echo $ts $nthreads $maxthreads $memcur $mempeak $user_usec $system_usec $usage_usec
	sleep 1
    done
) | tee $UNIQNAME.data

# Generate GNUPlot source and run it
cat <<EOF >$UNIQNAME.gnuplot
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
set y2label "Memory"
set xlabel "Running Time"

plot data u 1:2 w l t 'Threads', \
     data u 1:3 w l t 'Max Threads', \
     data u 1:4 w l t 'Memory' axes x1y2, \
     data u 1:5 w l t 'Max Memory' axes x1y2

set output outcname

set title 'CPU Time Usage'

unset y2label
unset y2tics
set ylabel "CPU Time"
plot data u 1:6 w l t "User", \
     data u 1:7 w l t "System", \
     data u 1:8 w l t "Total"
EOF
gnuplot $UNIQNAME.gnuplot

# Done
wait $PIPID
echo DONE
