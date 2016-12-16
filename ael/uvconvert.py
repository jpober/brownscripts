#!/bin/env python
import os,sys,glob,optparse,subprocess


### General script for using pyuvdata to convert among uv file formats. Submits array tasks to slurm, one for each file to convert.
# Usage:   uvconvert.py [-i <in_format>] [-o <out_format(s)>] [<obs_IDs>]
#  All parameters are optional. If run by itself, uvconvert will find FHD save files and convert them to both uvfits and MIRIAD formats.
#  Available formats in are FHD, uvfits, and MIRIAD (case insensitive). Available formats out are uvfits or miriad.
#  Example:    uvconvert.py -i "fhd,miriad"
###


o = optparse.OptionParser()

o.set_usage('uvconvert.py [-i <in_format>] [-o <out_format(s)>] [<obs_IDs>]')

o.add_option('-i', dest='in_format', help="Input file formats (miriad, fhd, uvfits)", default='fhd')
o.add_option('-o', dest='out_format', help="Output file formats (miriad, uvfits)", default='miriad,uvfits')
o.add_option('-c',  dest='clobber', action='store_true', help="Overwrite existing files", default=False)
o.add_option('-s', '--seq', dest='batch', action='store_false', help='Run as a batch job (default) or set false to run in sequence', default=True)


opts,args = o.parse_args(sys.argv[1:])


def check_file(fname, fmt):
	if fmt == 'uvfits':
		return'.'.join(fname.split('.')[0:-1])   #remove uvfits extension
	if fmt == 'miriad':
		contents=os.listdir(fname)
		if 'visdata' in contents and 'vartable' in contents: return fname
		else:
			print "Invalid MIRIAD file:", fname
			sys.exit(1)

in_formats=opts.in_format.lower().split(",")
for fmt in in_formats:
    if fmt == 'fhd':
	## Check for the vis_data directory
	if not os.path.isdir('vis_data'):
	    print 'Missing vis_data directory'
	    sys.exit(1)
	else:
	    if len(args) == 0:
		print 'Seeking all FHD save files for conversion'
		files=map(os.path.basename,glob.glob('vis_data/*flags.sav'))
		args=["_".join(obs.split("_")[0:-1]) for obs in files]
	    obslist=args
    if fmt == 'miriad':
	if len(args) == 0:
	    print 'Seeking all local MIRIAD files for conversion'
	    args = [d for d in os.listdir('.') if os.path.isdir(d)]
        obslist=[check_file(l,'miriad') for l in args]

    if fmt == 'uvfits':
	if len(args) == 0:
	    print 'Seeking all local uvfits files for conversion'
	    args = [d for d in glob.glob('./*.uvfits')]
        obslist=[check_file(l,'uvfits') for l in args]

## We now have a list "obslist" of all obs IDs to be converted.
## If they refer to MIRIAD or uvfits files, they must be in the current directory. If they refer to FHD save files, they must be in a "vis_data" subdir of the current directory.


out_formats=opts.out_format.lower().split(',')

ofiles=[]
for fmt in out_formats:
## Check for the proper output directories. If they don't exist, create them. If they do exist, check for existing files in them of the correct format.
## IF the option 'clobber' is not in place, remove found files from the obslist
    if fmt == 'miriad':
	if not os.path.isdir('MIRIAD'):
	    os.mkdir('MIRIAD')
	fils = [d for d in os.listdir('MIRIAD') if os.path.isdir("MIRIAD/"+d)]
	fils = map(os.path.normpath, fils)
	ofiles.append( map(os.path.basename, fils) )
    if fmt == 'uvfits':
	if not os.path.isdir("UVFITS"):
	   os.mkdir('UVFITS')
	fils=map(os.path.basename,[d for d in glob.glob('./UVFITS/*.uvfits')]) 
	ofiles.append([check_file(l,'uvfits') for l in fils])

if len(out_formats) == 2:
	of_mir, of_uvf = ofiles[0], ofiles[1]     #Don't flatten yet
	ofiles = [of for of in of_mir if of in of_uvf]
else:
        ofiles = [f for sub in ofiles for f in sub]  # flatten

## If outputting to both formats, we need to be sure that each file exists in both directories. If it's in one but not the other, still run the task but set the writer not to clobber. This is achieved by seeking unique entries in the ofiles list.
#if len(out_formats) == 2:
#	uniq = list(set(ofiles))  #Remove duplicates

if not opts.clobber:  obslist = [ob for ob in obslist if ob not in ofiles]


### We now have a list of obs_ids to convert, and formats to convert it to. This information will be passed along to the respective array jobs.
### One array job per input file.

mem='15G'
time='01:30:00'

Nf= len(obslist)
if Nf == 0:
	print "No files to convert"
	sys.exit(1)

obslist=":".join(obslist)

if opts.batch:

	batstr = 'sbatch --array=0-'+str(Nf-1) + ' -o \'slurm-%a.out\' --mem='+mem+' -t '+time+\
		' /gpfs_home/alanman/extra_scripts/uvconvert_job.py --formats='+",".join(in_formats)+":"+",".join(out_formats)+\
		' --clobber='+str(opts.clobber)+\
		' '+obslist

	subprocess.call(batstr,shell=True)
else:
	pids=[]
	for i in range(Nf):
	   batstr =' /gpfs_home/alanman/extra_scripts/uvconvert_job.py --formats='+",".join(in_formats)+":"+",".join(out_formats)+\
		' --clobber='+str(opts.clobber)+\
		' '+obslist

	   p=subprocess.Popen(batstr, stderr=open('subprocess-'+str(i)+'.err','w+'), stdout=open('subprocess-'+str(i)+'.out','w+'), shell=True, env=dict(os.environ,**{'SLURM_ARRAY_TASK_ID':str(i)}))
	   pids.append(p.pid)
	print pids



