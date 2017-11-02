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
sep=('0,1','-1,1','1,1')
opts,args = o.parse_args(sys.argv[1:])
color_iter = itertools.cycle(['navy', 'c', 'cornflowerblue', 'gold',
                              'darkorange'])

#BLS=os.popen("python ./capo/pspec_pipeline/getbls.py --sep='0,1' -C psa6240_FHD /Users/Josh/Desktop/PSA64Noise/DATA/zen.2456242.30605.uvcRREcACOTUcA").read()
#BLS = BLS.split(',')
base = uvdata.miriad.Miriad()
#data = n.zeros((14,203,len(BLS)))
try:
    base.read_miriad(args[0])
except:
    pass

def StackData(BLS):
    data = n.zeros((14,203,len(BLS)))
    for b in range(len(BLS)):
        b1,b2 = BLS[b].split('_')
        bls = base.baseline_array==base.antnums_to_baseline(b1,b2)
        try:
            if len(base.data_array[bls,0,:,0]) == 0:
                continue
            else:
                data[:,:,b] = base.data_array[bls,0,:,0]
        except:
            pass
    return data

def DefFeat(data):
    mus = n.mean(data,2)
    MAX = n.max(data,2)
    MED = n.median(data,2)
    MIN = n.min(data,2)
    SKU = skew(data,axis=2)
    KURT = kurtosis(data,axis=2)
    stds = n.var(data,2)
    phs = n.mean(n.angle(data),2)
    phsKurt = kurtosis(n.angle(data),2)
    shps = data.shape
    X = n.zeros((shps[0]*shps[1],7))
    X[:,0] = mus.reshape(shps[0]*shps[1])
    X[:,1] = MAX.reshape(shps[0]*shps[1])
    X[:,2] = phsKurt.reshape(shps[0]*shps[1])
    X[:,3] = MIN.reshape(shps[0]*shps[1])
    X[:,4] = SKU.reshape(shps[0]*shps[1])
    X[:,5] = phs.reshape(shps[0]*shps[1])
    X[:,6] = stds.reshape(shps[0]*shps[1])
    return X

def GMMask(X):
    dpgmm = mixture.BayesianGaussianMixture(n_components=2,covariance_type='full').fit(X)
    Y = dpgmm.predict(X)
    _Y_ = Y.reshape(14,203)
    MASK = n.zeros((14,203))
    MASK[_Y_==1] = 1.
    return MASK


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

#plot_results(X, dpgmm.predict(X), dpgmm.means_, dpgmm.covariances_, 0,'BGMM')
#pl.show()

def AddFlags(MASK,BLS):
    for b in range(len(BLS)):
        b1,b2 = BLS[b].split('_')
        bls = base.baseline_array==base.antnums_to_baseline(b1,b2)
        try:
            if base.flag_array[bls,0,:,0].shape[0] == 0:
                continue
            else:
                base.flag_array[bls,0,:,0] += MASK.astype(bool)
        except:
            pass


for s in sep:
    print s
    BLS=os.popen("python ~/capo/pspec_pipeline/getbls.py --sep="+s+" -C psa6240_FHD /users/jkerriga/data/jkerriga/PGBH/even/Pzen.2456242.30605.uvcRREcACOTUcHPAB").read()
    print BLS
    BLS = BLS.split(',')
    D = StackData(BLS)
    X = DefFeat(D)
    M = GMMask(X)
    AddFlags(M,BLS)
base.write_miriad(args[0]+'r')
