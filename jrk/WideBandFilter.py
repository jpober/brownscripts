import numpy as n
import uvdata
from numpy import fft
import aipy as a
import optparse, sys, os
import pylab as pl
from numpy.random import randn
import scipy.signal

o = optparse.OptionParser()

opts,args = o.parse_args(sys.argv[1:])
#window='blackman-harris'
#dir='/users/jkerriga/data/jkerriga/8DayLST/even/Pzen.2456242.30605.uvcRREcACOTUcHPA'
print args
#mir = uvdata.miriad.Miriad()
for files in args:
    mir = uvdata.miriad.Miriad()
    try:
        mir.read_miriad(files)
    except:
        pass
    mir.flag_array[mir.data_array==0j] = True
    d1 = mir.data_array[:,0,:,0]*n.logical_not(mir.flag_array[:,0,:,0]).astype(float)
    #bh = scipy.signal.chebwin(203,300,sym=True)
    #bh = scipy.signal.tukey(151,0.1,sym=True)
    #bh = a.dsp.gen_window(203,window='blackman-harris').astype(complex)
    bh = n.ones(203,dtype=complex)
    bh1 = n.zeros(203,dtype=complex)
    bh1[0:203] = bh
    d2 = n.zeros((d1.shape[0],203),dtype=complex)
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
    # Find proper delays to filter
    freqs = n.linspace(mir.freq_array[0][0],(D1.shape[1]*(mir.freq_array[0][1]-mir.freq_array[0][0])+mir.freq_array[0][0]),D1.shape[1])
    delays = n.fft.fftfreq(freqs.size,freqs[1]-freqs[0])
    delays = n.fft.fftshift(delays)
    delay_max = 30.0/(2.99*10**8)
    ### Moves in 15ns, which puts the sidelobes down to 10^0 for the +15ns so any structure added is *somewhat* negligible
    delay_bin = (n.abs(delays)<delay_max-10.0*10**(-9)).sum()/2
    #print delay_bin
    # Design delay filter
    w1 = n.zeros(n.array(D1_.shape),dtype=complex)#*10**(-14)
    w1[:,shp[1]/2-delay_bin:shp[1]/2+delay_bin+1] = n.ones((shp[0],2*delay_bin+1))

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
    
    #mir.flag_array[n.isnan(mir.data_array)] = True
    mir.epoch=2000.0
    mir.phase_center_epoch=2000.0
    mir.vis_units='Jy'
    mir.antenna_positions = n.zeros((63,3))
    mir.write_miriad(files+'F')

    del(mir)


