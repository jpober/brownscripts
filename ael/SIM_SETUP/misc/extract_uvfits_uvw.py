#!/bin/env python

import sys
import glob
from astropy.io import fits
import numpy as np


files = glob.glob('*.uvfits')

hdus = [fits.open(f)[0] for f in files]

uvw = []
uus, vvs, wws = [],[],[]
for uv in hdus:
	lines = zip(uv.data['UU'].tolist(), uv.data['VV'].tolist(), uv.data['WW'].tolist())
	lines = ["(%.2e , %.2e , %.2e) " % l for l in lines]
	uvw.append(lines)

uvw = np.array(uvw)
uvw = np.swapaxes(uvw,0,1)
uvw = [' '.join(l)+"\n" for l in uvw]

print np.array(uvw).shape

#print uvw

with open('UVWs.txt','w+') as f:
	f.writelines(uvw)

