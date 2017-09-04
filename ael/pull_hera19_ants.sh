#!/bin/bash

#SBATCH --mem=90G
#SBATCH -t 4:00:00
#SBATCH -n 15

pids=()
for obs in "$@"
do
     echo $obs
     pull_ants_pyuvdata.py $obs -a 80,104,96,64,53,31,65,88,9,20,89,43,105,22,81,10,72,112,97 &
     pids+=($!)
     len=${#pids[@]}
     if [[ $len -eq 30 ]]
     then
        echo ${pids[@]}
        wait
        pids=()
     fi
done
echo ${pids[@]}
wait
