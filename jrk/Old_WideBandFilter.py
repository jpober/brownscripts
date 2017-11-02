import numpy as n
import uvdata
from numpy import fft
import aipy as a
import optparse, sys, os
import pylab as pl
from numpy.random import randn

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
    #15:166
    mir.flag_array[mir.data_array==0j] = True
    d1 = mir.data_array[:,0,:,0]*n.logical_not(mir.flag_array[:,0,:,0]).astype(float)
    #nfreqs = (mir.flag_array[0,0,:,0]==False).sum()
    #nblts = (mir.flag_array[:,0,0,0]==False).sum()
    #print nfreqs,nblts
    #d1 = d1.reshape(nblts,nfreqs)
    #d1 = d1[:,0,:,0]
    bh = a.dsp.gen_window(203,window='blackman-harris')
    #bh = n.ones(151)
    d2 = n.zeros((d1.shape[0],d1.shape[1]+821),dtype=complex)
    d2[:,0:203] = d1*bh
    #BH = n.fft.fft(bh)
    #BH_ = n.fft.fftshift(BH)
    #f1 = n.logical_not(mir.flag_array[:,0,:,0]).astype(float)
    #F1 = n.fft.fft(f1,axis=1)
    #F1_ = n.fft.fftshift(F1,axes=1)
    #window = a.dsp.gen_window(141,window='blackman-harris')
    D1 = n.fft.fft(d2,axis=1)
    D1_ = n.fft.fftshift(D1,axes=1)
    shp = D1_.shape
    w1 = n.ones(n.array(D1_.shape))
    #print w1.shape
    #f = 1-mir.flag_array[:,0,:,0].astype(float)
    #f = window*f
    #F = n.fft.fft(f,axis=1)
    #F_ = n.fft.fftshift(F,axes=1)
    #pl.subplot(111)
    #pl.imshow(n.log10(n.abs(D1_)),aspect='auto')
    #pl.show()
    #mid =  a.dsp.gen_window(40,window='blackman-harris')
    #window = a.dsp.gen_window(203,window='blackman-harris')
    #print d2.shape
    freqs = n.linspace(mir.freq_array[0][0],(d2.shape[1]*(mir.freq_array[0][1]-mir.freq_array[0][0])+mir.freq_array[0][0]),d2.shape[1])
    #print freqs
    delays = n.fft.fftfreq(freqs.size,freqs[1]-freqs[0])
    delays = n.fft.fftshift(delays)
    delay_max = 30.0/(2.99*10**8)
    delay_bin = (n.abs(delays)<delay_max+15.0*10**(-9)).sum()/2
    print delay_bin
    w1[:,shp[1]/2-delay_bin:shp[1]/2+delay_bin+1] = n.zeros((shp[0],2*delay_bin+1))
    #w1[:,shp[1]/2-delay_bin] = 0.5
    #w1[:,shp[1]/2+delay_bin+1] = 0.5
    
    convD1 = D1_*w1
    #for i in range(d1.shape[0]):
    #    D1_[i,:] = D1_[i,:]*w1
    #D_clean,info = a.deconv.clean(D1,F,tol=10**(-9),maxiter=100)
    #D1_MDL = D1 - D_clean
    #F1 = n.fft.ifftshift(F1_)
    #f1_ = n.fft.ifft(F1)
    #BH = n.fft.ifftshift(BH_)
    #bh_ = n.fft.ifft(BH)
    DD1 = n.fft.ifftshift(convD1,axes=1)
    dd1 = n.fft.ifft(DD1,axis=1)
    dd1 = dd1[:,0:203]/bh
    #n.where(f1>0,dd1/(n.fft.ifft(n.fft.ifftshift(F1_))),0)
    #pl.imshow(n.log10(n.abs(D_)),aspect='auto')
    #pl.show()
    mir.data_array[:,0,:,0] = dd1
    #mir.flag_array[n.log10(n.abs(mir.data_array))>5] = True
    mir.flag_array[n.isnan(mir.data_array)] = True
    mir.epoch=2000.0
    mir.phase_center_epoch=2000.0
    mir.vis_units='Jy'
#mir.data_array = mir.data_array[:,:,54:128,:]
        #mir.vis_units='J'
        #mir.zenith_ra = mir.phase_center_ra*n.ones(len(mir.baseline_array))
        #mir.zenith_dec = mir.phase_center_dec*n.ones(len(mir.baseline_array))
        #mir.is_phased=False
    mir.write_miriad(files+'F')
#        del(mir)
    #except KeyError or TypeError:
    #for ant in range(64):
    #bsl = mir.antnums_to_baseline(41,49)
    #d = mir.data_array[bsl==mir.baseline_array,0,:,0]
    #    d = mir.data_array[:,0,:,0]
    #    D = n.fft.ifft(d)
    #    D_ = n.fft.ifftshift(D)
    #w = a.dsp.gen_window(141,window='blackman-harris')
    #w = 1-w
    #    w = n.ones(141)
    #    w[len(w)/2-10:len(w)/2+10] = 0.0 #00000000001
    #    for i in range(d.shape[0]):
    #        D_[i,:] = D_[i,:]*w
    #    DD = n.fft.fftshift(D_)
    #    dd = n.fft.fft(DD)
    #mir.data_array[bsl==mir.baseline_array,0,:,0] = d
    #    mir.data_array[:,0,:,0] = dd
        #mir.vis_units='J'
        #mir.is_phased=False
        #mir.zenith_dec = mir.phase_center_dec*n.ones(len(mir.baseline_array))
        #mir.zenith_ra = mir.phase_center_ra*n.ones(len(mir.baseline_array))
    #    mir.write_miriad(files+'F')
    del(mir)


