#!/usr/bin/python
import numpy as n
import uvdata
from numpy import fft
import aipy as a
import optparse, sys, os
from scipy import sparse
from numpy.random import randn
import os
import pylab as pl
import capo as C
from sklearn import mixture
from scipy.stats import skew
from scipy.stats import kurtosis
import matplotlib as mpl
from scipy import linalg
import itertools

o = optparse.OptionParser()
sep='(0,1)'
opts,args = o.parse_args(sys.argv[1:])
color_iter = itertools.cycle(['navy', 'c', 'cornflowerblue', 'gold',
                              'darkorange'])

BLS=os.popen("python ./capo/pspec_pipeline/getbls.py --sep='0,1' -C psa6240_FHD /Users/Josh/Desktop/PSA64Noise/DATA/zen.2456242.30605.uvcRREcACOTUcA").read()
BLS = BLS.split(',')
base = uvdata.miriad.Miriad()
data = n.zeros((14,203,len(BLS)))
for files in args:
#    print 'Adding..',files
    if files == args[0]:
        try:
            base.read_miriad(files)
        except:
            pass
        for b in range(len(BLS)):
#            print BLS[b]
            b1,b2 = BLS[b].split('_')
            bls = base.baseline_array==base.antnums_to_baseline(b1,b2)
            if base.data_array[bls,0,:,0].shape[0] == 0:
                continue
            else:
                data[:,:,b] = base.data_array[bls,0,:,0]
    else:
        print 'Error Will Robinson...Error'
        del(mir)
data = data[:,:,:]
mus = n.mean(data,2)
MAX = n.max(data,2)
MED = n.median(data,2)
MIN = n.min(data,2)
SKU = skew(data,axis=2)
KURT = kurtosis(data,axis=2)
stds = n.var(data,2)
phs = n.mean(n.angle(data),2)
phsKurt = kurtosis(n.angle(data),2)
idx1 = (n.abs(phsKurt)> 2.0).astype(int)
idx2 = (n.nan_to_num(n.log10(stds))>3.7).astype(int)
idx3 = (n.nan_to_num(n.log10(stds))>1.4).astype(int) 
#IDX = ((idx1+idx2+idx3)/3.).astype(bool)
IDX = ((idx1+idx1)/2).astype(bool)

shps = data.shape
print shps
X = n.zeros((shps[0]*shps[1],7))
X[:,0] = mus.reshape(shps[0]*shps[1])
X[:,1] = MAX.reshape(shps[0]*shps[1])
X[:,2] = phsKurt.reshape(shps[0]*shps[1])
X[:,3] = MIN.reshape(shps[0]*shps[1])
X[:,4] = SKU.reshape(shps[0]*shps[1])
X[:,5] = phs.reshape(shps[0]*shps[1])
X[:,6] = stds.reshape(shps[0]*shps[1])

dpgmm = mixture.BayesianGaussianMixture(n_components=2,covariance_type='full').fit(X)

def plot_results(X, Y_, means, covariances, index, title):
    splot = pl.subplot(1, 1, 1 + index)
    for i, (mean, covar, color) in enumerate(zip(
            means, covariances, color_iter)):
        if not n.any(Y_ == i):
            continue
        pl.scatter(X[Y_ == i, 0], X[Y_ == i, 6], .8, color=color,label=str(i))
    pl.legend()
    pl.xticks(())
    pl.yticks(())
    pl.title(title)

plot_results(X, dpgmm.predict(X), dpgmm.means_, dpgmm.covariances_, 0,'BGMM')
pl.show()

Y = dpgmm.predict(X)
print Y.shape
print Y.reshape(14,203)
_Y_ = Y.reshape(14,203)


MASK = n.zeros((14,203))
MASK[_Y_==1] = 1.

#MASK[_Y_==1] = 1.
#pl.figure()
#pl.subplot(311)
#pl.imshow(n.log10(n.abs(data[:,:,10])),aspect='auto',interpolation=None)
#pl.subplot(312)
#pl.imshow(MASK,aspect='auto',interpolation=None)
#pl.subplot(313)
#pl.imshow(n.log10(n.abs(data[:,:,10]*n.logical_not(MASK))),aspect='auto',interpolation=None)
#pl.show()
#print IDX.sum()
#pl.subplot(211)
#pl.imshow(n.log10(n.abs(data[:,10])),aspect='auto',interpolation=None)
#pl.subplot(212)
#pl.imshow(IDX.astype(float),aspect='auto',interpolation=None)
#pl.show()
for b in range(len(BLS)):
    b1,b2 = BLS[b].split('_')
    bls = base.baseline_array==base.antnums_to_baseline(b1,b2)
    if base.flag_array[bls,0,:,0].shape[0] == 0:
        continue
    else:
        print base.flag_array[bls,0,:,0][_Y_==1].shape
        base.flag_array[bls,0,:,0] += MASK.astype(bool)
base.write_miriad(args[0]+'r')
#pl.figure()
#pl.subplot(411)
#pl.scatter(mus,phs)
#pl.subplot(412)
#pl.scatter(n.abs(MAX),n.abs(MIN))
#pl.subplot(413)
#pl.scatter(n.abs(SKU),n.abs(mus))
#pl.subplot(414)
#pl.scatter(n.abs(SKU),n.abs(KURT))
#pl.show()
#pl.scatter(mus[IDX],n.log10(stds[IDX]),color='r')
#pl.scatter((mus),n.log10(stds))

