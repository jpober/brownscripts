#!/bin/env python
import os, glob, sys, re
import numpy as np 
from astropy.io import fits
import optparse


o = optparse.OptionParser()
o.set_usage('ExtractData.py [options] param1 param2...')
o.add_option('-s', '--search', dest='search', default=None)

opts, keys = o.parse_args(sys.argv[1:])

#(TODO -- add search options to select subset)

def stringSplitByNumbers(x):
    r = re.compile('(\d+)')
    l = r.split(x)
    return [int(y) if y.isdigit() else y for y in l]

fail=False
files = glob.glob('*.uvfits')

files = sorted(files, key= stringSplitByNumbers)

for f in files:
     uv = fits.open(f)[0]
     try:
       response = []
       for k in keys:
	     if k.lower() == 'jdate': response.append(str(k) + ' = '+ str(uv.header['PZERO5']))
	     if k.lower() == 'tzero': response.append(str(k) + ' = '+ str(uv.data['DATE'][0]))
	     elif k.lower() == 'nblts': response.append(str(k) + ' = '+ str(len(uv.data['DATE'])))
	     elif k.lower() == 'dfreq': response.append(str(k) + ' = '+ str(uv.header['CDELT4']))
	     elif k.lower() == 'sfreq': response.append(str(k) + ' = '+ str(uv.header['CRPIX4']))
             else: response.append(str(k) + " = " + str(uv.header[k]))
     except KeyError:
        fail=True
	break
     print f, ", ".join(response)


if fail:
    uv = fits.open(files[0])[0]
#    print uv.header
    avail = uv.header.keys()
    avail.append('JDate'); avail.append('Ntimes')
    print ','.join(avail)
