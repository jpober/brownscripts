#!/usr/bin/python
import numpy as n
import pyuvdata
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
from scipy import signal
from sklearn import cluster
from glob import glob
from scipy import ndimage
o = optparse.OptionParser()
opts,args = o.parse_args(sys.argv[1:])

#BLS=os.popen("python ./capo/pspec_pipeline/getbls.py --sep='0,1' -C psa6240_FHD /Users/Josh/Desktop/PSA64Noise/DATA/zen.2456242.30605.uvcRREcACOTUcA").read()
#BLS = BLS.split(',')
base = pyuvdata.miriad.Miriad()
#day = pyuvdata.miriad.Miriad()
dlist = glob('/users/jkerriga/data/jkerriga/HERA/zen.2457458.*.xx.HH.uvcUA')

print dlist
try:
    base.read_miriad(args[0],run_check=False)
#    for i in dlist:
#        print i
#        day.read_miriad(i)
#        if i == dlist[0]:
#            dayData = day.data_array[:,0,:,0]
#            dayTimes = day.lst_array
#        else:
#            dayData.vstack(day.data_array[:,0,:,0])
#            dayTimes.vstack(day.lst_array)
except:
    pass
dayData = []
dayTimes = []
for k in dlist:
    uv = a.miriad.UV(k)
    a.scripting.uv_selector(uv,'9_10','xx')
    for p,d,f in uv.all(raw=True):
        dayData.append(d)
        dayTimes.append(uv['lst'])
dayData = n.array(dayData)
dayTimes = n.array(dayTimes)

print n.shape(dayData),n.shape(dayTimes)
def haar(points, a):
    haarf = n.zeros(points)
    haarf[0:a/2] = 1.
    haarf[a/2:a] = -1.
    return haarf

def featArray(data,timesObs):
    freqs1 = n.linspace(100,200,n.shape(data)[1])
    sh = n.shape(data)
    CW1mean = n.zeros_like(data)
    for i in range(CW1mean.shape[1]):
        CW1 = cwt(n.abs(data[:,i]),haar,n.arange(1,10,1))
        CW1 = n.ma.masked_where(CW1==0,CW1)
        CW1mean[:,i] = n.ma.mean(n.abs(CW1),0)

    CT1mean = n.zeros_like(data)
    for j in range(CW1mean.shape[0]):
        CT1 = cwt(data[j,:],signal.morlet,n.arange(1,3,1))
        CT1 = n.ma.masked_where(CT1==0,CT1)
        CT1mean[j,:] = n.mean(n.abs(CT1),0)
    processed = ndimage.sobel(n.abs(data))
    X1 = n.zeros((sh[0]*sh[1],5))
    X1[:,0] = (n.real(data)).reshape(sh[0]*sh[1])
    X1[:,1] = (n.imag(data)).reshape(sh[0]*sh[1])
    #X1[:,2] = n.abs(CW1mean).reshape(sh[0]*sh[1])
    X1[:,2] = n.abs(CT1mean).reshape(sh[0]*sh[1])
#    X1[:,2] = n.log10(n.abs(processed)).reshape(sh[0]*sh[1])
    X1[:,3] = (n.array([timesObs]*sh[1])).reshape(sh[0]*sh[1])
    X1[:,4] = (n.array([freqs1]*sh[0])).reshape(sh[0]*sh[1])
    X1 = n.nan_to_num(X1)
#    for m in range(n.shape(X1)[1]):
#        X1[:,m] = X1[:,m]/X1[:,m].max()
    X1 = normalize(X1,norm='l2',axis=1)
    X1 = n.nan_to_num(X1)
    return X1


def featTargetArray(data,timesObs,maskedLST):
    freqs1 = n.linspace(100,200,1024)
    sh = n.shape(data)
    #CW1mean = n.zeros_like(data)
    #for i in range(CW1mean.shape[1]):
    #    CW1 = cwt(n.abs(data[:,i]),haar,n.arange(1,10,1))
    #    CW1 = n.ma.masked_where(CW1==0,CW1)
    #    CW1mean[:,i] = n.ma.mean(n.abs(CW1),0)

    CT1mean = n.zeros_like(data)
    for j in range(CW1mean.shape[0]):
        CT1 = cwt(n.abs(data[j,:]),haar,n.arange(1,3,1))
        CT1 = n.ma.masked_where(CT1==0,CT1)
        CT1mean[j,:] = n.mean(n.abs(CT1),0)
    X1 = n.zeros((sh[0]*sh[1],6))
    X1[:,0] = (n.real(data)).reshape(sh[0]*sh[1])
    X1[:,1] = (n.imag(data)).reshape(sh[0]*sh[1])
    #X1[:,2] = n.abs(CW1mean).reshape(sh[0]*sh[1])
    X1[:,2] = n.abs(CT1mean).reshape(sh[0]*sh[1])
    X1[:,3] = (n.array([timesObs]*sh[1])).reshape(sh[0]*sh[1])
    X1[:,4] = (n.array([freqs1]*sh[0])).reshape(sh[0]*sh[1])
    X1[:,5] = (maskedLST).reshape(sh[0]*sh[1])
    X1 = n.nan_to_num(X1)
   # for m in range(n.shape(X1)[1]):
   #     X1[:,m] = X1[:,m]/X1[:,m].max()
    X1 = normalize(X1,norm='l2',axis=0)
    X1 = n.nan_to_num(X1)
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

a1 = 9
a2 = 10
bsl_data = base.baseline_array==base.antnums_to_baseline(a1,a2)
data = base.data_array[bsl_data,0,:,0]
#if data.shape[0] == 0:
#    break
times = base.lst_array[bsl_data]
Xtrain = featArray(dayData,dayTimes)
#for i in range(1,15):
fit_param = 3
#    print i
gmm = mixture.BayesianGaussianMixture(n_components=fit_param,covariance_type='diag',n_init=1,max_iter=1000,init_params='random').fit(Xtrain)
del(Xtrain)

Xtarget = featArray(data,times)
#gmm = mixture.BayesianGaussianMixture(n_components=fit_param,covariance_type='full',n_init=1000,max_iter=10000,init_params='random').fit(Xtarget)
Nlabel = gmm.predict(Xtarget)
sh = n.shape(data)
Nlabel_,maxlabel = LabelMaker(Nlabel,sh)
#   score.append(gmm.score(Xtrain))
#db = cluster.DBSCAN(n_jobs=2,min_samples=10).fit(Xtrain)

MASK = n.zeros(n.shape(data)).astype(bool)
MASK[Nlabel_!=maxlabel] = True
base.flag_array[bsl_data,0,:,0] += MASK
base.write_miriad(args[0]+'r')
