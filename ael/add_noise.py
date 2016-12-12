#!/bin/env python

import numpy as np, ephem, uvdata.uvdata as uvd, sys
#import glob, matplotlib.pyplot as plt
from aipy.const import k   #Boltzmann constant in erg/K
import optparse,os


o = optparse.OptionParser()
o.set_usage('paper_add_noise.py [options] input_file(s).uvfits')

o.add_option('-o', dest='opath', help='Destination directory',default='./noise_files')
o.add_option('-s', dest='suffix', help='Optional suffix for output files', default=None)
o.add_option('--new', action='store_true', dest='newflag', help='Create a new file of purely noise.', default=False)
o.add_option('-N', dest='nfiles', help='Number of output files per input.', default=1)   #To create two or more noise realizations per file.
o.add_option('-A', dest='Area', help='Effective Area (cm^2)', default=55000.)  #cm^2
o.add_option('--Trec', dest='Trec', help='Receiver temperature (K)', default=100. ) #K
o.add_option('--Tref', dest='Tref', help='Sky temperature at 150 MHz (K)', default=400 ) #K
#o.add_option('-df', dest='df', help='Channel width (Hz)', default=492610.837 )  #Hz
#o.add_option('-dt', dest='dt', help='Integration length (seconds)', default=30.0 ) #sec

opts,args = o.parse_args(sys.argv[1:])

uv = uvd.UVData()

if not os.path.exists(opts.opath):
    os.makedirs(opts.opath)

#print opts.opath
#print args
#print opts.nfiles


for f in args:    #Loop over list of files
 uv.read_miriad(f)
 dt = uv.integration_time
 df = uv.channel_width
 obsid = f.split('.')[0]

 for n in range(int(opts.nfiles)):
   if not opts.suffix is None: obsid = obsid + opts.suffix

   if opts.nfiles > 1:
     if opts.nfiles == 2:
	suff = 'even' if n == 0 else 'odd'
	ofile = obsid+"_"+suff
     else:
	ofile = obsid + "_n" + str(n)
   else:
     ofile = obsid + str("_noise")

   uv.uvw_array = uv.uvw_array.astype(np.float64)   #Needed for MIRIAD's type requirements... temporary fix

   #uv.write_miriad(opts.opath+'/'+obsid)

   opath = opts.opath + "/" + ofile


   Tsys = float(opts.Tref)*np.power(uv.freq_array[0]/(150e6),-2.6) + float(opts.Trec)*np.ones(uv.Nfreqs)
   sigs =  k*Tsys/(opts.Area*np.sqrt(df*dt))*1e23/np.sqrt(2) 	# Jy
   for blt in range(uv.Nblts):
       for ns, sig in enumerate(sigs):
#TODO -- Change option here -- Add noise, or replace with noise?
	  if opts.newflag:
	   	uv.data_array[blt,0,ns,:] = np.random.normal(0,sig,(uv.Npols)) + 1j*np.random.normal(0,sig,(uv.Npols))
	  else:
	   	uv.data_array[blt,0,ns,:] += np.random.normal(0,sig,(uv.Npols)) + 1j*np.random.normal(0,sig,(uv.Npols))

#   uv.data_array[:,:

#   print 'Sigma = ', sig
#   uv.data_array = np.random.normal(0, sig, (uv.Nbls * uv.Ntimes, uv.Nspws, uv.Nfreqs,uv.Npols)) \
#		  + 1j*np.random.normal(0, sig, (uv.Nbls * uv.Ntimes, uv.Nspws, uv.Nfreqs,uv.Npols))
#   print uv.data_array.shape

#   uv.write_uvfits(opath+".uvfits")

   uv.write_miriad(opath, clobber=True)
