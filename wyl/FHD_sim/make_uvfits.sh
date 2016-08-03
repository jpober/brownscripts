#!/bin/bash

### Submit the sbatch array command to generate a set of uvfits files

fname='hexpro1061311664'
paramfile='mwa_128_cal.py'
calfile='hex_sim_1.py'
N=94   #Number of files
mem='50G'
time='2:00:00'


#sbatch -o /dev/null  --array=1-$Narr --mem=$mem -t $time -n 3 --export=N=$N,fname=$fname,paramfile=$paramfile,calfile=$calfile zeros_job.sh
sbatch  --array=1-$N --mem=$mem -t $time -n 3 --export=N=$N,fname=$fname,paramfile=$paramfile,calfile=$calfile zeros_job.sh
