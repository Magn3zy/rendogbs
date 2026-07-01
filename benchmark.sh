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
echo TS THREADS MEMORY USER SYSTEM TOTAL
while ps $PIPID >/dev/null 2>&1 ; do
    ts=$(date +%s)
    nthreadss=$(wc -l /sys/fs/cgroup/system.slice/docker-$LONGID.scope/cgroup.threads 2>/dev/null)
    nthreads=${nthreadss%% *}
    mempeak=$(cat /sys/fs/cgroup/system.slice/docker-$LONGID.scope/memory.peak 2>/dev/null)
    parse_cpustat $(cat /sys/fs/cgroup/system.slice/docker-$LONGID.scope/cpu.stat 2>/dev/null)
    echo $ts $nthreads $mempeak $user_usec $system_usec $usage_usec
    sleep 1
done
      
# Done
wait $PIPID
echo DONE
