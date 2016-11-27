import numpy as n
import uvdata
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
filters = C.dspec.wedge_width_by_bl(aa, uv['sdf'], chans, offset=-50.0)
del(uv)
print args
print filters[(41,49)]
uthresh,lthresh = filters[(41,49)]
#mir = uvdata.miriad.Miriad()
for files in args:
    mir = uvdata.miriad.Miriad()
    try:
        mir.read_miriad(files)
    except:
        pass
    #mir.flag_array[mir.data_array==0j] = True
    d1 = mir.data_array[:,0,:,0]*n.logical_not(mir.flag_array[:,0,:,0]).astype(float)
    #bh = scipy.signal.chebwin(203,400,sym=True)
    #bh = scipy.signal.tukey(203,0.1,sym=True)
    #bh = a.dsp.gen_window(151,window='blackman-harris').astype(complex)
    bh = n.ones(151,dtype=complex)
    bh1 = n.zeros(chans,dtype=complex)
    bh1[15:166] = bh
    d2 = n.zeros((d1.shape[0],chans),dtype=complex)
    d2[:,0:203] = d1[:,:]
    #d1[:,14:167] = d1[:,14:167]*bh
    d2 = d2*bh1
    #flg = n.logical_not(mir.flag_array[:,0,:,0]).astype(float)
    #flg = n.zeros((d1.shape[0],1024),dtype=complex)
    #flg[:,0:203] = f1
    # Transform Data
    D1 = n.fft.fft(d2,axis=1)
    D1_ = n.fft.fftshift(D1,axes=1)
    shp = D1_.shape
    print shp
    # Find proper delays to filter
    freqs = n.linspace(mir.freq_array[0][0],(D1.shape[1]*(mir.freq_array[0][1]-mir.freq_array[0][0])+mir.freq_array[0][0]),D1.shape[1])
    delays = n.fft.fftfreq(freqs.size,freqs[1]-freqs[0])
    delays = n.fft.fftshift(delays)
    delay_max = 30.0/(2.99*10**8)
    ### Moves in 15ns, which puts the sidelobes down to 10^0 for the +15ns so any structure added is *somewhat* negligible
    #delay_bin = (n.abs(delays)<delay_max+45.0*10**(-9)).sum()/2
    #print delay_bin
    # Design delay filter
    w1 = n.ones(n.array(D1_.shape),dtype=complex)#*10**(-14)
    w1[:,uthresh:lthresh] = 0.
    w1 = n.fft.fftshift(w1,axes=1)
    #pl.plot(w1[0,:])
    #pl.show()
    print w1.shape
    #w1[:,fb[1]:fb[0]] = n.zeros((shp[0],fb[0]-fb[1]))
    #w1[:,shp[1]/2 -delay_bin:shp[1]/2+delay_bin] = 0#n.zeros((shp[0],2*delay_bin+1),dtype=complex)
    #W1 = n.fft.fftshift(w1)
    #W1_ = n.fft.fft(W1)
    #Filter delay data
    convD1 = w1*D1_
    #dd1  = n.zeros(n.array(W1_.shape),dtype=complex)
    #for i in range(d1.shape[0]):
    #    dd1[i,:] = n.convolve(bh*flg[i,:],W1_[i,:],mode='same')
    #amp = n.convolve(bh1,W1_[0,:],mode='same')
    DD1 = n.fft.ifftshift(convD1,axes=1)
    dd1 = n.fft.ifft(DD1,axis=1)
    #dd1[:,14:167] = dd1[:,14:167]/bh
    #dd1 = n.where(flg>0,dd1/bh1,0)
    
    mir.data_array[:,0,:,0] = dd1[:,0:203]
    #mir.data_array[:,0,:,0] = mir.data_array[:,0,:,0] - dd1[:,0:203]
    
    mir.flag_array[n.isnan(mir.data_array)] = True
    mir.flag_array[mir.data_array == 0j] = True
    mir.phase_center_epoch=2000.0
    mir.vis_units='Jy'
    mir.antenna_positions = n.zeros((63,3))
    mir.write_miriad(files+'F_')

    del(mir)


