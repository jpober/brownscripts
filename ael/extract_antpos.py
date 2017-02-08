#!/bin/env python

#SBATCH -J extract_antpos
#SBATCH --mem=10G
#SBATCH --time=1:00:00


### Given a uvfits file, extracts the relative positions of antennas from the baseline array using pyuvdata

from uvdata import UVData
import sys, numpy as np
import pickle


uvd = UVData()

uvd.read_uvfits(sys.argv[1])

uvd.unphase_to_drift()

Nants = uvd.Nants_telescope


bls={}
for ant2 in range(Nants):
	# ant1 = 0 --- (positions relative to ant1)
	bl_ind = 2048 * (0+1) + ant2+1 + 2**16
	print bl_ind
	loc=np.where(uvd.baseline_array==bl_ind)
	print loc
	bls[ant2] = uvd.uvw_array[loc][0]

for ant in bls.keys():
	bls[ant] = bls[ant].tolist()

ofilename = sys.argv[1].split('.')[:-1]
ofilename = ".".join(ofilename)+"_antpos.npz"

with open(ofilename, 'wb') as ofile:
	pickle.dump(bls, ofile)

print type(bls)
