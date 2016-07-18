#!/bin/bash

### Submit the sbatch array command to generate a set of uvfits files

fname='testhex1061311664'
paramfile='psa2hex_sim_cal.py'
calfile='hex_sim_1.py'
N=4   #Number of files
mem='20G'
time='1:00:00'


#sbatch -o /dev/null  --array=1-$Narr --mem=$mem -t $time -n 3 --export=N=$N,fname=$fname,paramfile=$paramfile,calfile=$calfile zeros_job.sh
sbatch  --array=1-$N --mem=$mem -t $time -n 3 --export=N=$N,fname=$fname,paramfile=$paramfile,calfile=$calfile zeros_job.sh
