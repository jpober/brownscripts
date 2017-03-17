#!/bin/env python

#SBATCH -J grow_miriad
#SBATCH --mem=30G
#SBATCH --time=1:00:00

# Take a MIRIAD uvdata file and extract the antenna positions 

import aipy as a, numpy as np
import sys, optparse, shutil
import pickle, os, subprocess
from astropy import constants as const
from astropy.coordinates import Angle
from astropy import units as u

generic_cal_filepath = '/users/alanman/GrowMIRIAD/generic_cal.py'

o = optparse.OptionParser()
o.set_usage('extract_miriad_antpos.py [options] <input file>')  #Miriad files only!

o.add_option('-o', dest='ofile', help='Output file')
o.add_option('--no-clobber', dest='clobber', action='store_false', help='Overwrite existing files (default = yes)', default=True)
#o.add_option('-s', dest='suffix', help='Optional suffix for output files', default=None)
#o.add_option('--new', action='store_true', dest='newflag', help='Create a new file of purely noise.', default=True)

opts,args = o.parse_args(sys.argv[1:])

if opts.ofile is None:
	opts.ofile=args[0]+"_new"


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
uvi.select('auto',-1,-1,include=False)
ntimes=uvi['ntimes']
npol=uvi['npol']
nchan=uvi['nchan']
nants=uvi['nants']

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


### Make a cal file.

##Convert the uvws dictionary to a numpy array:
#uvw_orig = np.array(data_accumulator[pol_list[0]][0][0])*const.c.to('m/ns').value
#uvws = np.array([uvw_orig for n in range(nants)])
#uvws = np.insert(uvws,0,[0.,0.,0.],axis=0)

da = np.array(data_accumulator[pol_list[0]])

uvw_all = da[:,0]
uvws = np.zeros((nants,3))
for j in range(1,nants):
       try:
            uvws[j,:] = (-1)*uvw_all[np.where((da[:,2] == 0) & (da[:,3] == j))][0] * const.c.to('m/ns').value
#       except IndexError:
#            pass   #Skip failures
	## -1 ---> Flip from the sky to the ground


bls = {}
top_names = ['top_x', 'top_y', 'top_z']
for i,a in enumerate(uvi['antnums']):
	bls[a] = {}
	for j in range(3):
		bls[int(a)][top_names[j]] = uvws[i,j]

lat = uvi['latitud']  # in units of radians
lon = uvi['longitu']
tloc = (lat, lon)

del(uvi)

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

