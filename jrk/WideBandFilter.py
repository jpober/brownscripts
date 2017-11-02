import numpy as n
import pyuvdata
from numpy import fft
import aipy as a
import optparse, sys, os
import pylab as pl
from numpy.random import randn
import scipy.signal
import capo as C

o = optparse.OptionParser()

opts,args = o.parse_args(sys.argv[1:])
#window='blackman-harris'
#dir='/users/jkerriga/data/jkerriga/8DayLST/even/Pzen.2456242.30605.uvcRREcACOTUcHPA'
chans = 203
uv = a.miriad.UV(args[0])
aa = a.cal.get_aa('psa6240_FHD', uv['sdf'], uv['sfreq'], chans)
filters = C.dspec.wedge_width_by_bl(aa, uv['sdf'], chans, offset=15.0)
print filters.keys()
del(uv)
print args
#print filters[(41,49)]
#uthresh,lthresh = filters[(41,49)]
#mir = uvdata.miriad.Miriad()
for files in args:
    mir = pyuvdata.miriad.Miriad()
    try:
        mir.read_miriad(files)
    except:
        pass

    for i in filters.keys():
        print i
        bsl = mir.antnums_to_baseline(i[0],i[1])
        bidx = mir.baseline_array==bsl
        d1 = mir.data_array[bidx,0,:,0]*n.logical_not(mir.flag_array[bidx,0,:,0]).astype(float)
        if d1.shape[0] == 0:
            continue
        #bh1 = n.ones(203,dtype=complex)
        bh1 = a.dsp.gen_window(chans,window='blackman-harris')#n.zeros(chans,dtype=complex)
        #bh1[15:166] = bh
        #d2 = n.zeros((d1.shape[0],chans),dtype=complex)
        #d2[:,0:203] = d1[:,:]
        #d1 = d1*bh1            
        # Find proper delays to filter
        # Design delay filter
        w1 = n.zeros(n.array(d1.shape),dtype=complex)#*10**(-14)
        uthresh,lthresh = filters[i]
        print uthresh,lthresh
        w1[:,uthresh:lthresh] = 1. 
        w1 = n.fft.fftshift(w1)
        zeroLen = (w1[0,:]== 0).sum()
        w1[:,w1[0,:]== 0 ] = 1 -scipy.signal.tukey(zeroLen,alpha=0.3)
        D1 = n.fft.fft(d1,axis=1)
        D1_ = n.fft.fftshift(D1,axes=1)
        #w1 = n.fft.fftshift(w1)
        #pl.plot(n.abs(w1[0,:]))
        #pl.show()
        #W1 = n.fft.fft(w1,axis=1)
        #W1_ = n.fft.fftshift(W1,axes=1)
        Fdata = n.fft.ifft(n.fft.ifftshift(D1_*w1))
        #Fdata = []
        #for i in range(d1.shape[0]):
        #    Fdata.append(n.convolve(d1[i,:]*bh1,W1_[i,:],mode='same'))
        #print n.shape(Fdata)
        try:
            mir.data_array[bidx,0,:,0] = n.array(Fdata)
        except:
            pass
    mir.flag_array[n.isnan(mir.data_array)] = True
    mir.flag_array[mir.data_array == 0j] = True
    mir.phase_center_epoch=2000.0
    mir.vis_units='Jy'
    #mir.antenna_positions = n.zeros((63,3))
    mir.write_miriad(files+'F')

    del(mir)


