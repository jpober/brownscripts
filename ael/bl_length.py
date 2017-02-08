#!/bin/env python

import aipy as a
import numpy as np
import os, sys
from shutil import copy2

herapath="/users/alanman/SIM_SETUP/HERA/hsa7458_v000_HH.py"
hera37path="/users/alanman/SIM_SETUP/HERA/hera37_cal.py"
paperpath="/users/alanman/SIM_SETUP/PAPER/psa6622_v003_paper128.py"
mwapath="/users/alanman/SIM_SETUP/MWA/mwa_128_cal.py"

if len(sys.argv) < 3:
	print "Usage: bl_length.py <instrument> <baseline(s)>" 

instr = sys.argv[1]
bls=sys.argv[2:]


if ',' in ','.join(bls): bls=','.join(bls).split(',')  #In case input is in the form of a comma-separated list

if instr.lower()=='hera':
	cal=os.path.basename(herapath)
	cal=cal.split('.')[0]
	copy2(herapath, '.')
elif instr.lower()=='hera37':
	cal=os.path.basename(hera37path)
	cal=cal.split('.')[0]
	copy2(hera37path, '.')
elif instr.lower()=='paper':
	cal=os.path.basename(paperpath)
	cal=cal.split('.')[0]
	copy2(paperpath, '.')
elif instr.lower()=='mwa':
	cal=os.path.basename(mwapath)
	cal=cal.split('.')[0]
	copy2(mwapath, '.')


aa = a.cal.get_aa(cal,.1,.1,.1)
for bl in bls:
	i,j = map(int,bl.split("_"))
	bl=aa.get_baseline(i,j,src='z')
	bl=bl*a.const.len_ns / 100.
	print "_".join(map(str,[i,j]))," : ",np.linalg.norm(bl)

os.remove(cal+".py")
os.remove(cal+".pyc")
