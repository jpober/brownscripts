#!/bin/env python

### Add together UV data_arrays with matching shapes.

import sys, os
import numpy as np
from copy import deepcopy
import optparse
import re
import pickle as pkl
import subprocess


o = optparse.OptionParser()
o.set_usage('e.g. miriad_arithmetic.py "dir1+0.5*dir2=outdir" -k <search key> -d <base_directory>')

o.add_option('-k', dest='search', help='Search string for locating files',default='*')
o.add_option('-d', dest='basedir', help='Folder containing input directories.',default='.')

opts,args = o.parse_args(sys.argv[1:])

command=args[0]    #

terms=re.split('([\+\-\*\/\=])',command)
if not '=' in terms:
	print 'Warning: No output directory specified. Using default.'
	terms.append('=')
	terms.append('outdir')

#print command

dir_inds=[]

#Build a list of input and output directories
for i,t in enumerate(terms):
    try:
	float(t)
    except ValueError:
	if not t in ['+','-','*','=','/']:
	     dir_inds.append(i)
             if terms[i-1]=='=': break   # The directory after the equals sign must be the output.
n_in=len(dir_inds)-1

#Find input files.

os.chdir(opts.basedir)

mirlist=[]   #list of miriad files. (n_in,n_files)

j=0
indirs=[]
for i in dir_inds[:-1]:
    if not os.path.isdir(terms[i]):
	print "Error: " + terms[i] + " could not be found."
	sys.exit()
    else:
	mirlist.append(next(os.walk(terms[i]))[1])
	indirs.append(terms[i])
	terms[i]='uvds['+str(j)+'].data_array'    			#Building the command string for sub-tasks
	j += 1
	#TODO -- Remove miriad files from the list that don't match the search criterion (-k)

#Make output directory, if necessary
if not os.path.isdir(terms[-1]):
	os.mkdir(terms[-1])

#Verify matching files
for i in range(n_in-1):
	if not set(mirlist[i]) == set(mirlist[i+1]):
		print "Warning: MIRIAD file lists do not match! Selecting those that do."
		break

setlist = map(set, mirlist)
mirlist[0] = list(set.intersection(*setlist))

command_str=''.join(terms)


opts='/gpfs_home/alanman/extra_scripts/opts.pkl'


pkl.dump({'mirlist':mirlist[0], 'indirs':indirs, 'command':command_str}, open(opts, 'w+'))

Nf=len(mirlist[0])

mem='15G'
time='01:30:00'

batstr = 'sbatch --array=0-'+str(Nf-1)+' -o \'slurm-%a.out\' --mem='+mem+' -t '+time+' /gpfs_home/alanman/extra_scripts/miriad_arithmetic_job.py'

subprocess.call(batstr,shell=True)


#Start an array job of smaller scripts that will do the actual arithmetic.
	##  In a pickle object, store -- file list, input directories (in order), output directory, command_string
	## The job script will unpickle the information and, based on the array_job_id, select a filename to work on.
	## Command_str has the form [arithmetic]=[output]. In the jobscript, split by '=' and use the arithmetic part in an eval statement. This should work provided you've loaded the various miriad files into UVData objects in an array called "uvds"
	## The job script should then create a new file with the correct file name in [output].


sys.exit()

uvd0 = UVData()
uvd1 = UVData()
uvdout = UVData()

uvd0.read_miriad(args[0])
data0 = uvd0.data_array
del(uvd0)

uvd1.read_miriad(args[1])
uvdout = deepcopy(uvd1)
uvdout.data_array = uvd1.data_array + data0
del(uvd1)
uvdout.write_miriad(opts.opath)

#uv0 = fits.open(sys.argv[1])
#uv1 = fits.open(sys.argv[2])


