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
chans = 203.
uv = a.miriad.UV(args[0])
aa = a.cal.get_aa('psa6240_FHD', uv['sdf'], uv['sfreq'], chans)
filters = C.dspec.wedge_width_by_bl(aa, uv['sdf'], chans, offset=0.0)
print filters.keys()
del(uv)
print args
#print filters[(41,49)]
#uthresh,lthresh = filters[(41,49)]
#mir = uvdata.miriad.Miriad()
def Fzero(flags,direction):
    ct=0
    flags = n.logical_not(flags)
    try:
        flags = flags[0,:]
    except:
        pass
    while flags[ct] != 1:
        if direction == 1:
            ct+=1
        else: 
            ct-=1
        if ct == 203 or ct == -203:
            break
        if n.sum(flags) == 0:
            break
    zidx = ct
    return zidx

for files in args:
    mir = pyuvdata.miriad.Miriad()
    try:
        mir.read_miriad(files)
    except:
        pass

    for i in filters.keys():
        bsl = mir.antnums_to_baseline(i[0],i[1])
        bidx = mir.baseline_array==bsl
        f = mir.flag_array[bidx,0,:,0]
        d1 = mir.data_array[bidx,0,:,0]#*n.logical_not(mir.flag_array[bidx,0,:,0]).astype(float)
        try:
            idx1 = Fzero(f,1)
            idx2 = Fzero(f,-1)
            idx2 = n.abs(idx2+204)
        except:
            continue
        bh1 = n.ones(idx2-idx1,dtype=complex)
        #bh1 = a.dsp.gen_window(idx2-idx1,window='blackman-harris')#n.zeros(chans,dtype=complex)
        #bh1[15:166] = bh
        #d2 = n.zeros((d1.shape[0],chans),dtype=complex)
        #d2[:,0:203] = d1[:,:]
        d1[:,idx1:idx2] = d1[:,idx1:idx2]*bh1

        # Transform Data
        D1 = n.fft.fft(d1,axis=1)
        D1_ = n.fft.fftshift(D1,axes=1)
    
        # Find proper delays to filter
        # Design delay filter
        w1 = n.zeros(n.array(D1_.shape),dtype=complex)#*10**(-14)
        #ants = mir.baseline_to_antnums(i)
        uthresh,lthresh = filters[i]
        w1[:,uthresh:lthresh] = 1.
        w1 = n.fft.fftshift(w1,axes=1)
        #Filter delay data
        convD1 = w1*D1_
        DD1 = n.fft.ifftshift(convD1,axes=1)
        dd1 = n.fft.ifft(DD1,axis=1)
        mir.data_array[bidx,0,idx1:idx2,0] = dd1[:,idx1:idx2]/bh1
    
    mir.flag_array[n.isnan(mir.data_array)] = True
    mir.flag_array[mir.data_array == 0j] = True
    mir.phase_center_epoch=2000.0
    mir.vis_units='Jy'
    #mir.antenna_positions = n.zeros((63,3))
    mir.write_miriad(files+'F')

    del(mir)


