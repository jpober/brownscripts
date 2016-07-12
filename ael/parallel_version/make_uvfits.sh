#!/bin/bash

### Submit the sbatch array command to generate a set of uvfits files

fname='mwa128'
paramfile='golden_1061311664_sim_1.py'
calfile='hera19_cal.py'
N=1   #Number of files
mem='10G'
time='1:00:00'


#sbatch -o /dev/null  --array=1-$N --mem=$mem -t $time -n 3 --export=N=$N,fname=$fname,paramfile=$paramfile,calfile=$calfile zeros_job.sh
sbatch  --array=1-$N --mem=$mem -t $time -n 3 --export=N=$N,fname=$fname,paramfile=$paramfile,calfile=$calfile zeros_job.sh
