#!/bin/bash

host=`uname -n`

if [[ ${hostname:8:9} == 'sulis.hpc' ]]; then
    rank=$PMI_RANK
else
    rank=$OMPI_COMM_WORLD_RANK
fi

suffix=`printf "%3.3d" $rank`

export OMP_NUM_THREADS=1

ulimit -s unlimited

#export KMP_BLOCKTIME=10000
#export KMP_LIBRARY=turnaround
#export KMP_AFFINITY="verbose,granularity=thread,scatter"
#export KMP_AFFINITY="verbose,compact"
#export KMP_MONITOR_STACKSIZE=100m
# これがないと落ちる
#export KMP_STACKSIZE=8000m
export KMP_STACKSIZE=250m

 
#./md.g.mpi inp/HP2_repT.inp 1> log.$suffix 2> err.$suffix
./md.g.mpi inp/HP2_repT.inp data/restart 1> log.$suffix 2> err.$suffix
