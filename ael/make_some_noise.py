#!/bin/env python

import numpy as np
import os, optparse
import sys
import pickle as pkl
import subprocess

### Submit the sbatch array command to generate a set of uvfits files

o = optparse.OptionParser()
o.set_usage('make_some_noise [options] <input files>')  #Miriad files only

o.add_option('-o', dest='opath', help='Destination directory',default='./noise_files')
o.add_option('-i', dest='ipath', help='Input directory', default=None)
o.add_option('--new', action='store_true', dest='newflag', help='Create a new file of purely noise.', default=False)
o.add_option('-N', dest='nfiles', help='Number of output files per input.', default=1)   #To create two or more noise realizations per file.
o.add_option('-A', dest='Area', help='Effective Area (cm^2)', default=55000.)  #cm^2
o.add_option('-I', '--instrument', dest='instr', help='Instrument (HERA/PAPER/MWA). Overrides A and Trec')
o.add_option('--scale', type="float", help='Arbitrary scaling on the noise level.', default=1)
o.add_option('--Trec', dest='Trec', help='Receiver temperature (K)', default=100. ) #K
o.add_option('--Tref', dest='Tref', help='Sky temperature at 150 MHz (K)', default=400 ) #K
o.add_option('--oname', help='Output file name key', default=None)
#o.add_option('-df', dest='df', help='Channel width (Hz)', default=492610.837 )  #Hz
#o.add_option('-dt', dest='dt', help='Integration length (seconds)', default=30.0 ) #sec

opts,args = o.parse_args(sys.argv[1:])


#args = MIRIAD files
if len(args) < 1:
	print 'No files given.'
	sys.exit()

basepath, file0 = os.path.split(args[0])
if basepath=='': basepath='.'
basepath = os.path.realpath(basepath)
if opts.ipath == None: opts.ipath = basepath

optfile = open('/gpfs_home/alanman/extra_scripts/opts.pkl', 'w')
pkl.dump(opts,optfile)

args = map(os.path.basename, args)
filelist = ":".join(args)

mem='15G'
time='01:30:00'


#Nf=$(ls -d  $basepath/$key* | wc -l)    #The d flag is necessary to prevent it from finding files within the MIRIAD 'folders'

Nf = len(args)

logfile_base = "slurm"
if not opts.oname is None:
     logfile_base += "_"+opts.oname
copynum=1
if os.path.exists(logfile_base+'-0.out'):
    logfile_name = logfile_base + str(copynum)
    while os.path.exists(logfile_name+"-0.out"):
        copynum+=1
        logfile_name = logfile_base + str(copynum)
else:
    logfile_name=logfile_base


batstr = 'sbatch --array=0-'+str(Nf-1) + ' -o \''+logfile_name+'-%a.out\' --mem='+mem+' -t '+time+' /gpfs_home/alanman/extra_scripts/add_noise.py '+filelist

#print batstr

subprocess.call(batstr,shell=True)


#TODO -- Remove the opts.pkl file after the full set of processes is completed.

#os.remove('/gpfs/data/jpober/alanman/NOISE/opts.pkl')

#sbatch  --array=0-$(($Nf - 1)) -o "slurm-"$obs_base"_%a.out" --mem=$mem -t $time --export=N=$N,suffix=$suffix,obs=$obs_base,path=$basepath,outpath=$outpath /gpfs/data/jpober/alanman/NOISE/noise_job.sh

#sbatch  --array=0-$(($Nf - 1)) -n 5 -o "slurm-"$key"_%a.out" --mem=$mem -t $time --export=N=$N,path=$basepath,key=$key,outpath=$outpath /gpfs/data/jpober/alanman/NOISE/noise_job.sh
