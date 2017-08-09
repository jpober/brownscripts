#!/bin/env python

#SBATCH -J extract_antpos
#SBATCH --mem=30G
#SBATCH --time=3:00:00


### Given a uvfits or miriad file, extracts the relative positions of antennas from the baseline array using pyuvdata

from pyuvdata import UVData
import sys, numpy as np, os
import pickle
import pylab as p
from astropy.coordinates import Angle
from astropy import units as u

uvd = UVData()

generic_cal_filepath = '/users/alanman/brownscripts/ael/generic_cal.py'

if os.path.isfile(sys.argv[1]):
	uvd.read_uvfits(sys.argv[1])
elif os.path.isdir(sys.argv[1]):
	uvd.read_miriad(sys.argv[1])
else:
	print "Invalid file type"
	sys.exit()

if uvd.phase_type == 'phased': uvd.unphase_to_drift()

Nants = uvd.Nants_telescope

bls={}

ant1 = 0
for ant2 in np.arange(Nants):
        bl_ind = 2048 * (ant1+1) + ant2+1 + 2**16
	loc=np.where(uvd.baseline_array==bl_ind)
	if len(loc[0]) == 0:
	     bl_ind = 2048 * (ant2+1) + ant1+1 + 2**16
	     loc=np.where(uvd.baseline_array==bl_ind)
        if len(loc[0]) == 0:
            import IPython; IPython.embed() 
	bls[ant2] = uvd.uvw_array[loc][0]

for ant in bls:
	bls[ant] = bls[ant].tolist()

if '.' in sys.argv[1]:
	ofile_base = '.'.join(sys.argv[1].split('.')[:-1])
else:
	ofile_base = sys.argv[1]

ofilename = ofile_base+"_antpos.pkl"

tloc = uvd.telescope_location_lat_lon_alt

lat, lon = Angle(tloc[0],u.radian).signed_dms, Angle(tloc[1], u.radian).signed_dms
lat = str(lat[0]*lat[1])+":"+str(lat[2])+str(lat[3])   # signed D:M:S format
lon = str(lon[0]*lon[1])+":"+str(lon[2])+str(lon[3]) 

loc = (lat, lon)

antpos={}
#Convert to the format used in cal files.
top_names = ['top_x', 'top_y', 'top_z']
for ant in bls:
	tmp = bls[ant]
	antpos[ant] = {}
	for j in range(3): antpos[ant][top_names[j]] = tmp[j]

prms= { 'loc' : loc, 'antpos_ideal' : antpos }

#with open(ofilename, 'wb') as ofile:
#	pickle.dump(prms, ofile)

with open(generic_cal_filepath,'r') as cfile:
	data = cfile.readlines()

ofile_base=ofile_base.replace('-','_')
ofile_base = ofile_base.replace('.','')   #Ensure a valid cal file name
data[data.index('INSERT_PRMS_HERE\n')] = "prms = "+ repr(prms)  #Convert to a string
with open(ofile_base+"_cal.py", 'w') as cfile:
	cfile.writelines(data)


antpos = np.tile([0.,0.,0.],len(bls)).reshape((len(bls),3))

for i,k in enumerate(bls):
	antpos[i] = bls[k]

np.savez(ofile_base+"_antpos",antpos=antpos)

#p.scatter(antpos[:,0], antpos[:,1])
#p.show()
