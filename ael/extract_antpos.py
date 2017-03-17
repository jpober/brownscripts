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

generic_cal_filepath = '/users/alanman/GrowMIRIAD/generic_cal.py'

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

#uvd.baseline_array = uvd.antnums_to_baseline(uvd.ant_1_array, uvd.ant_2_array, attempt256=True)
#uvd.baseline_array = uvd.antnums_to_baseline(uvd.ant_2_array, uvd.ant_1_array)

#flip_inds = False
#bl_ind = 2048 * (0+1) + 1+1 + 2**16
#loc=np.where(uvd.baseline_array==bl_ind)
#if len(loc[0]) == 0:
#	flip_inds = True

#for bl_ind in uvd.baseline_array:
#	print ((bl_ind - 2**16) - 1)/(2048.) - 1    #, ' : ', ((bl_ind - 2**16) - 2048 - 1)
#sys.exit()

#print len(uvd.baseline_array)/654

ant1 = 0
for ant2 in np.arange(Nants):
	# ant1 = 0 --- (positions relative to ant1)
#	if flip_inds:
#	     bl_ind = 2048 * (ant2+1) + 0+1 + 2**16
#	else:
#	     bl_ind = 2048 * (0+1) + ant2+1 + 2**16
        bl_ind = 2048 * (ant1+1) + ant2+1 + 2**16
	loc=np.where(uvd.baseline_array==bl_ind)
	if len(loc[0]) == 0:
	     bl_ind = 2048 * (ant2+1) + ant1+1 + 2**16
	     loc=np.where(uvd.baseline_array==bl_ind)
#	print len(loc[0]),  ant2
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

#Convert to the format used in cal files.
top_names = ['top_x', 'top_y', 'top_z']
for ant in bls:
	tmp = bls[ant]
	bls[ant] = {}
	for j in range(3): bls[ant][top_names[j]] = tmp
for i,a in enumerate(uvd.antenna_numbers):
	antpos[a] = {}
	for j in range(3):
		antpos[int(a)][top_names[j]] = bls[i,:][j]

prms= { 'loc' : loc, 'antpos' : bls }

with open(ofilename, 'wb') as ofile:
	pickle.dump(prms, ofile)

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
	

p.scatter(antpos[:,0], antpos[:,1])
p.show()
