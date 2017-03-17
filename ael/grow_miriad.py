#!/bin/env python

#SBATCH -J grow_miriad
#SBATCH --mem=30G
#SBATCH --time=1:00:00

# Take a MIRIAD uvdata file containing a single baseline and "grow" it to mimic an array of a given number of baselines.

import aipy as a, numpy as np
import sys, optparse, shutil
import pickle, os, subprocess
from astropy import constants as const

generic_cal_filepath = '/users/alanman/brownscripts/ael/generic_cal.py'

o = optparse.OptionParser()
o.set_usage('grow_miriad.py [options] <input file>')  #Miriad files only!

o.add_option('-o', dest='ofile', help='Output file')
o.add_option('-N', '--nbls', dest='nbls', help='Number of baselines for new file')
o.add_option('--nants', dest='nants', help='Instead of defining the number of baselines, define a number of antennas')
o.add_option('--no-clobber', dest='clobber', action='store_false', help='Overwrite existing files (default = yes)', default=True)
#o.add_option('-s', dest='suffix', help='Optional suffix for output files', default=None)
#o.add_option('--new', action='store_true', dest='newflag', help='Create a new file of purely noise.', default=True)

opts,args = o.parse_args(sys.argv[1:])

if opts.ofile is None:
	opts.ofile=args[0]+"_new"

#if not opts.nants is None:
#	nants = int(opts.nants)
#	opts.nbls = int(nants*(nants-1)/2.)   #Exclude autocorrs
#else:
#	opts.nbls = 1

if not opts.nbls is None:
	opts.nbls = int(opts.nbls)
	nants = opts.nbls+1
elif not opts.nants is None:
	nants = int(opts.nants)
	print nants
#	if not nants%2 == 0:
#		print "Error: Number of antennas must be even"
#		sys.exit()
#	else:
	opts.nbls = int(opts.nants)-1  #Nants must be even
else:
	print "Please specify a number of antennas or  baselines"
	sys.exit()


if opts.clobber:
  try:
	shutil.rmtree(opts.ofile)
  except OSError:
	pass

ifile = args[0]

#if '.' in opts.ofile:
#	ofile_base = ".".join(os.path.basename(opts.ofile).split('.')[:-1])
#else:
ofile_base = opts.ofile   # Assume no periods in the final result.

uvi = a.miriad.UV(ifile)
uvo = a.miriad.UV(opts.ofile, status='new')
uvo.init_from_uv(uvi)
uvi.select('auto',-1,-1,include=False)
ntimes=uvi['ntimes']
npol=uvi['npol']
nchan=uvi['nchan']


data_accumulator = {}
pol_list = []
for (uvw, t, (i, j)), d, f in uvi.all(raw=True):
    try:
        cnt = uvi['cnt']
    except(KeyError):
        cnt = np.ones(d.shape, dtype=np.float)
    ra = uvi['ra']
    dec = uvi['dec']
    lst = uvi['lst']
    source = uvi['source']
    try:
        data_accumulator[uvi['pol']].append([uvw, t, i, j, d, f, cnt,
                                            ra, dec])
    except(KeyError):
        data_accumulator[uvi['pol']] = [[uvw, t, i, j, d, f, cnt,
                                        ra, dec]]
        pol_list.append(uvi['pol'])


## Loop over data_accumulator, writing to file.
# For each baseline, duplicate the correct number of times.

uvo['antnums'] = np.arange(nants).astype(np.float64)
antenna_names = np.arange(nants).astype(np.str)
ant_name_flt = np.array([int(elem.encode("hex"), 16) for elem in antenna_names]).astype(np.float64)
#uvo.add_var('antnames', 'd')
uvo['antnames'] = ant_name_flt
uvo['nants'] = nants
uvo['nblts'] = opts.nbls*ntimes
uvo['nbls'] = opts.nbls


for pol in data_accumulator:
	for pd in data_accumulator[pol]:
	    #bl_arr = data_accumulator[pol][t]
	    d = pd[4]
	    flags = pd[5]
	    i, j = 0, 1
	    for viscnt in range(opts.nbls):
		preamble = (pd[0], pd[1], (i,j))
		uvo['ra'] = pd[7]
		uvo['dec'] = pd[8]
		uvo['pol'] = pol
		uvo.write(preamble, d, flags)
		j += 1

del uvi


### Make a cal file.

#cmd = 'extract_antpos.py '+ opts.ofile
#out = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#(stdout, stderr) = out.communicate()
#print stderr.decode().split()

#out.wait()

##Convert the uvws dictionary to a numpy array:
uvw_orig = np.array(data_accumulator[pol_list[0]][0][0])*const.c.to('m/ns').value
uvws = np.array([(-1)*uvw_orig for n in range(nants)])
uvws = np.insert(uvws,0,[0.,0.,0.],axis=0)

bls = {}
top_names = ['top_x', 'top_y', 'top_z']
for i,a in enumerate(uvo['antnums']):
	bls[a] = {}
	for j in range(3):
		bls[int(a)][top_names[j]] = uvws[i,:][j]


from astropy.coordinates import Angle
from astropy import units as u

lat = uvo['latitud']  # in units of radians
lon = uvo['longitu']
tloc = (lat, lon)

del(uvo)

lat, lon = Angle(tloc[0],u.radian).signed_dms, Angle(tloc[1], u.radian).signed_dms
lat = str(lat[0]*lat[1])+":"+str(lat[2])+str(lat[3])   # signed D:M:S format
lon = str(lon[0]*lon[1])+":"+str(lon[2])+str(lon[3]) 

loc = (lat, lon)

prms= { 'loc' : loc, 'antpos_ideal' : bls }

ofilename = ofile_base+"_antpos.pkl"

with open(ofilename, 'wb') as ofile:
	pickle.dump(prms, ofile)



ap_file = ofile_base+"_antpos.pkl"
prms = pickle.load(open(ap_file,'r'))

with open(generic_cal_filepath,'r') as cfile:
	data = cfile.readlines()

ofile_base=ofile_base.replace('-','_')
ofile_base = ofile_base.replace('.','')   #Ensure a valid cal file name
data[data.index('INSERT_PRMS_HERE\n')] = "prms = "+ repr(prms)  #Convert to a string
with open(ofile_base+"_cal.py", 'w') as cfile:
	cfile.writelines(data)

#p1 = subprocess.Popen('cp '+generic_cal_filepath+' ./'+file_base+"_cal.py", shell=True)
#p1.wait()



#(stdout, stderr) = out.communicate()

