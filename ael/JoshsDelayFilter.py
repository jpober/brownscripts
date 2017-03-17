#!/bin/env python

#SBATCH -J delay_filter
#SBATCH --mem=50G
#SBATCH --time=8:00:00

import numpy as n
from pyuvdata import UVData
from numpy import fft
import aipy as a
import optparse, sys, os
import pylab as pl
from numpy.random import randn
import scipy.signal
import capo as C

o = optparse.OptionParser()
o.add_option('-c','--cal', dest='cal', help='Calibration file')

opts,args = o.parse_args(sys.argv[1:])
#window='blackman-harris'
#dir='/users/jkerriga/data/jkerriga/8DayLST/even/Pzen.2456242.30605.uvcRREcACOTUcHPA'
uv = a.miriad.UV(args[0])
chans = uv['nchan']
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], chans)
filters = C.dspec.wedge_width_by_bl(aa, uv['sdf'], chans, offset=0.0)
#print filters.keys()
del(uv)
print args
#print filters[(41,49)]
#uthresh,lthresh = filters[(41,49)]
#mir = uvdata.miriad.Miriad()
for files in args:
    mir = UVData()
    try:
        mir.read_miriad(files)
    except:
        pass

    for i in filters.keys():
        print i
        bsl = mir.antnums_to_baseline(i[0],i[1])
        bidx = mir.baseline_array==bsl
        d1 = mir.data_array[bidx,0,:,0] #*n.logical_not(mir.flag_array[bidx,0,:,0]).astype(float)
        d2 = mir.data_array[bidx,0,:,1] #*n.logical_not(mir.flag_array[bidx,0,:,0]).astype(float)

#        bh = n.ones(151,dtype=complex)
	bh = a.dsp.gen_window(chans,window='blackman-harris')
#	bh = n.ones(chans,dtype=complex)
#        bh1 = n.zeros(chans,dtype=complex)
#        bh1[15:166] = bh
#        d2 = n.zeros((d1.shape[0],chans),dtype=complex)
#        d2[:,0:203] = d1[:,:]
#        d2 = d2*bh
        d1 = d1*bh
        d2 = d2*bh

        # Transform Data
        #D1 = n.fft.fft(d2,axis=1)
        D1 = n.fft.fft(d1,axis=1)
        D1_ = n.fft.fftshift(D1,axes=1)    
        D2 = n.fft.fft(d2,axis=1)
        D2_ = n.fft.fftshift(D2,axes=1)

        # Find proper delays to filter
        # Design delay filter
        w1 = n.ones(n.array(D1_.shape),dtype=complex)#*10**(-14)
        uthresh,lthresh = filters[i]
        w1[:,uthresh:lthresh] = 0.
        w1 = n.fft.fftshift(w1,axes=1)
        #Filter delay data
        convD1 = w1*D1_
        convD2 = w1*D2_
        DD1 = n.fft.ifftshift(convD1,axes=1)
        dd1 = n.fft.ifft(DD1,axis=1)
        DD2 = n.fft.ifftshift(convD2,axes=1)
        dd2 = n.fft.ifft(DD2,axis=1)
        mir.data_array[bidx,0,:,0] = dd1[:,0:chans]
        mir.data_array[bidx,0,:,1] = dd2[:,0:chans]
    
    mir.flag_array[n.isnan(mir.data_array)] = True
    mir.flag_array[mir.data_array == 0j] = True
    mir.phase_center_epoch=2000.0
    mir.vis_units='Jy'
#    mir.antenna_positions = n.zeros((63,3))
    mir.write_miriad(files+'F')

    del(mir)


