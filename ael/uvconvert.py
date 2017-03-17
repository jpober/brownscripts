#!/bin/env python
import os,sys,glob,optparse,subprocess
import numpy as np

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
o.add_option('-m', '--model', dest='model', action='store_true', help='If converting from fhd files, check if _model_vis_** files exist and if so, convert them separately.', default=False)


opts,args = o.parse_args(sys.argv[1:])


def check_file(fname, fmt):
	if fmt == 'uvfits':
		if fname.endswith('.uvfits'):
		      return '.'.join(fname.split('.')[0:-1])   #remove uvfits extension
		else:
		      return fname
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
    mod = '/model' if opts.model else ''
    if fmt == 'miriad':
	if not os.path.isdir('MIRIAD'):
	    os.mkdir('MIRIAD')
	if opts.model and not os.path.isdir('MIRIAD/model'):
	    os.mkdir('MIRIAD/model')
	fils = [d for d in os.listdir('MIRIAD'+mod) if os.path.isdir("MIRIAD/"+d)]
	fils = map(os.path.normpath, fils)
	ofiles.append( map(os.path.basename, fils) )
    if fmt == 'uvfits':
	if not os.path.isdir("UVFITS"):
	   os.mkdir('UVFITS')
	if opts.model and not os.path.isdir('UVFITS/model'):
	    os.mkdir('UVFITS/model')
	fils=map(os.path.basename,[d for d in glob.glob('./UVFITS'+mod+'/*.uvfits')]) 
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

obslist = [ob for ob in obslist if not ob == ""]

### We now have a list of obs_ids to convert, and formats to convert it to. This information will be passed along to the respective array jobs.
### One array job per input file.

mem='30G'
time='06:30:00'

Nf= len(obslist)
if Nf == 0:
	print "No files to convert"
	sys.exit(1)
else:
	if 'miriad' in in_formats: testfile=obslist[0]
	elif 'uvfits' in in_formats: testfile=obslist[0]+'.uvfits'
	else: testfile='vis_data/'+obslist[0]+"_vis_XX.sav"
	size = subprocess.check_output('du -hs '+testfile+' | cut -f 1',shell=True)
	#size=p.communicate()

	scale=size.strip()[-1]
	value=np.log10(float(size.strip()[0:-1]))
#	exp= np.ceil(value)
	if scale=='K': value += 3
	if scale=='M': value += 6
	if scale=='G': value += 9

	if 'fhd' in in_formats: value += 1     #(Account for the fact that data are split among fhd save files)
	#mem is specified in megabytes
	if value < 4: mem='5G'
	else:
		mem = 7*np.power(10,value-6)   # MB
		if mem < 1e4: mem = 1.5e4    #Less than 10G... ask for more.
		mem = str(int(np.floor(mem)))+"M"
obslist=":".join(obslist)
print mem

if opts.batch:
	n=0
	if opts.model: model_exists=['True', 'False']
	else: model_exists = ['False']
	for model in model_exists:
	   mod = '_'+str(n) if opts.model else ''
	   n += 1
	   batstr = 'sbatch --array=0-'+str(Nf-1) + ' -o \'slurm'+mod+'-%a.out\' --mem='+mem+' -t '+time+\
		' /gpfs_home/alanman/extra_scripts/uvconvert_job.py --formats='+",".join(in_formats)+":"+",".join(out_formats)+\
		' --clobber='+str(opts.clobber)+\
		' --model='+model+\
		' '+obslist

	   subprocess.call(batstr,shell=True)
else:
	pids=[]
	if opts.model: model_exists=['True', 'False']
	else: model_exists = ['False']
	for model in model_exists:
	  for i in range(Nf):
	    batstr =' /gpfs_home/alanman/extra_scripts/uvconvert_job.py --formats='+",".join(in_formats)+":"+",".join(out_formats)+\
		' --clobber='+str(opts.clobber)+\
		' --model='+model+\
		' '+obslist

	    p=subprocess.Popen(batstr, stderr=open('subprocess-'+str(i)+'.err','w+'), stdout=open('subprocess-'+str(i)+'.out','w+'), shell=True, env=dict(os.environ,**{'SLURM_ARRAY_TASK_ID':str(i)}))
	    pids.append(p.pid)
	print pids



