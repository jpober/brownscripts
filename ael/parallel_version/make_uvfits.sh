#!/bin/bash

### Submit the sbatch array command to generate a set of uvfits files

fname='hera37'
paramfile='HERA_sim_2.py'
calfile='hera37_cal.py'
N=94   #Number of files
mem='10G'
time='1:00:00'


sbatch -p jpober-test -o /dev/null --array=0-$N --mem=$mem -t $time -n 3 --export=N=$N,fname=$fname,paramfile=$paramfile,calfile=$calfile zeros_job.sh
