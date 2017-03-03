#!/usr/bin/python
import numpy as n
import uvdata
import aipy as a
import optparse, sys, os
from sklearn.preprocessing import normalize
import os
from scipy.signal import cwt
import pylab as pl
import capo as C
from sklearn import mixture
from scipy.stats import skew
from scipy.stats import kurtosis
from scipy import linalg

o = optparse.OptionParser()
sep=('0,1','-1,1','1,1')
opts,args = o.parse_args(sys.argv[1:])

#BLS=os.popen("python ./capo/pspec_pipeline/getbls.py --sep='0,1' -C psa6240_FHD /Users/Josh/Desktop/PSA64Noise/DATA/zen.2456242.30605.uvcRREcACOTUcA").read()
#BLS = BLS.split(',')
base = uvdata.miriad.Miriad()

try:
    base.read_miriad(args[0])
except:
    pass

def haar(points, a):
    haarf = n.zeros(points)
    haarf[0:a/2] = 1.
    haarf[a/2:a] = -1.
    return haarf

def featArray(data,timesObs):
    freqs1 = n.linspace(100,200,203)
    sh = n.shape(data)
    CW1mean = n.zeros_like(data)
    for i in range(CW1mean.shape[1]):
        CW1 = cwt(n.abs(data[:,i]),haar,n.arange(1,10,1))
        CW1 = n.ma.masked_where(CW1==0,CW1)
        CW1mean[:,i] = n.ma.mean(n.abs(CW1),0)

    CT1mean = n.zeros_like(data)
    for j in range(CW1mean.shape[0]):
        CT1 = cwt(n.abs(data[j,:]),haar,n.arange(1,10,1))
        CT1 = n.ma.masked_where(CT1==0,CT1)
        CT1mean[j,:] = n.mean(n.abs(CT1),0)
    X1 = n.zeros((sh[0]*sh[1],4))
    X1[:,0] = (n.real(data)).reshape(sh[0]*sh[1])
    X1[:,1] = (n.imag(data)).reshape(sh[0]*sh[1])
    #X1[:,2] = n.abs(CW1mean).reshape(sh[0]*sh[1])
    #X1[:,3] = n.abs(CT1mean).reshape(sh[0]*sh[1])
    X1[:,2] = (n.array([timesObs]*sh[1])).reshape(sh[0]*sh[1])
    X1[:,3] = (n.array([freqs1]*sh[0])).reshape(sh[0]*sh[1])
    X1 = n.nan_to_num(X1)
    X1 = normalize(X1,norm='l2',axis=0)
    return X1

def LabelMaker(newlabels,sh):
    newlabels = newlabels.reshape(sh[0],sh[1])
    ml1 = 0
    ml1num = 0
    for t in n.unique(newlabels):
        if (newlabels==t).sum()>ml1num:
            ml1 = t
            ml1num = (newlabels==t).sum()
    return newlabels,ml1


def AddFlags(maskedArray,baseline):
    for b in range(len(BLS)):
        b1,b2 = BLS[b].split('_')
        bls = base.baseline_array==base.antnums_to_baseline(b1,b2)
        if base.flag_array[bls,0,:,0].shape[0] == 0:
            continue
        else:
            base.flag_array[bls,0,:,0] += MASK.astype(bool)

for s in sep:
    print s
    BLS=os.popen("python ~/capo/pspec_pipeline/getbls.py --sep="+s+" -C psa6240_FHD /users/jkerriga/data/jkerriga/PGBH/even/Pzen.2456242.30605.uvcRREcACOTUcHPAB").read()
    print BLS
    BLS = BLS.split(',')
    for b in BLS:
        a1,a2 = b.split('_')
        bsl = base.baseline_array==base.antnums_to_baseline(a1,a2)
        data = base.data_array[bsl,0,:,0]
        if data.shape[0] == 0:
            continue
        times = base.time_array[bsl]
        Xarr = featArray(data,times)
        gmm = mixture.GaussianMixture(n_components=3,covariance_type='diag',n_init=10,max_iter=1000).fit(Xarr)
        Nlabel = gmm.predict(Xarr)
        sh = n.shape(data)
        Nlabel_,maxlabel = LabelMaker(Nlabel,sh)
        MASK = n.zeros(n.shape(data)).astype(bool)
        MASK[Nlabel_!=maxlabel] = True
        base.flag_array[bsl,0,:,0] += MASK
base.write_miriad(args[0]+'r')
