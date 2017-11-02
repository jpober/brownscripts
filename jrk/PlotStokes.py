import os
import numpy as n
from scipy import signal
from sklearn.preprocessing import normalize
from sklearn.decomposition import PCA
from glob import glob
import pyuvdata
import pylab as pl
import optparse, sys, os
import aipy as a
from sklearn.externals import joblib

def loadFullDay(day,pol):
#    HERAlist = glob('/Users/josh/Desktop/HERA/data/zen.2457458.*.xx.HH.uvcUA')
    HERAlist = glob('/users/jkerriga/data/jkerriga/PaperData/even_rect/sep0,1/lst.2456'+str(day)+'.*.*.*HBAL')
    HERAlist = n.sort(HERAlist)
    HERAdata = []
    times = []
    for k in HERAlist:
        uvHERA = a.miriad.UV(k)
        a.scripting.uv_selector(uvHERA, '41_49', pol)
        for p,d,f in uvHERA.all(raw=True):
            HERAdata.append(d)
            times.append(uvHERA['lst'])
    HERAdata = n.array(HERAdata)
    times = n.array(times)
    return HERAdata,times

def featArray(data):
    sh = n.shape(data)
    freqs = n.linspace(0,sh[1],sh[1])
    NNvar = n.zeros_like(data)
    dvar = n.var(n.abs(data))
    for i in range(sh[0]):
        for j in range(sh[1]):
            #samples = []
            #for p in range(100):
            k = n.random.randint(-1,1,size=1000)
            l = n.random.randint(-1,1,size=1000)
            #    try:
            #samples = n.abs(data[k+i,l+j])
            #    except:
            #        pass
            NNvar[i,j] = n.var(n.abs(data[k+i,l+j]))
    X1 = n.zeros((sh[0]*sh[1],3))
    X1[:,0] = (n.real(data)).reshape(sh[0]*sh[1])
    X1[:,1] = (n.imag(data)).reshape(sh[0]*sh[1])
    #X1[:,2] = (n.log10(n.abs(NNvar)) - n.median(n.log10(n.abs(NNvar)))).reshape(sh[0]*sh[1])
    NNvar = NNvar - n.median(NNvar)
    X1[:,2] = (n.log10(n.abs(NNvar))).reshape(sh[0]*sh[1])
    #X1[:,3] = (n.array([freqs]*sh[0])).reshape(sh[0]*sh[1])
    #X1[:,4] = (n.array([times]*sh[1])).reshape(sh[0]*sh[1])
    X1[n.abs(X1)>10**100] = 0
    for m in range(X1.shape[1]):
        X1[:,m] = X1[:,m]/n.abs(X1[:,m]).max()
    X1 = n.nan_to_num(X1)
    return X1

o = optparse.OptionParser()
opts,obs = o.parse_args(sys.argv[1:])

for o in range(242,243):
    print o
    data,times = loadFullDay(o,'I')
    #datay,times = loadFullDay(o,'yy')
    sh = n.shape(data)

#dataQ = datax-datay

def FFT(data):
    win = a.dsp.gen_window(203,window='blackman-harris')
    DATA = n.fft.fft(data)
    DATA_ = n.fft.fftshift(DATA)
    return DATA_
DATA_ = FFT(data)
from scipy.signal import csd
CSD = n.zeros((100,203)).astype(complex)
#pl.imshow(n.log10(n.abs(DATA_)),aspect='auto',interpolation='none')
for i in range(100):
    #CSD = csd(data[i,:],data[i+1,:],window='blackmanharris')
    D1 = FFT(data[i,:])
    D2 = FFT(data[i+1,:])
    CSD[i,:] = D1*n.conj(D2)
    #pl.plot(n.log10(n.abs(CSD)))
pl.subplot(211)
pl.imshow(n.log10(n.abs(DATA_[0:100,:]*n.conj(DATA_[0:100,:]))),aspect='auto',interpolation='none',vmax=4)
pl.subplot(212)
pl.imshow(n.log10(n.abs(CSD)),aspect='auto',interpolation='none',vmax=4)
pl.show()

