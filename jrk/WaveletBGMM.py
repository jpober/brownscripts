#!/usr/bin/python
import numpy as n
import pyuvdata
import aipy as a
from scipy import signal
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
#sep=('0,1','-1,1','1,1')
## because its HERA
sep=('0,1')
BLS = ['9_10']
opts,args = o.parse_args(sys.argv[1:])

#BLS=os.popen("python ./capo/pspec_pipeline/getbls.py --sep='0,1' -C psa6240_FHD /Users/Josh/Desktop/PSA64Noise/DATA/zen.2456242.30605.uvcRREcACOTUcA").read()
#BLS = BLS.split(',')
#try:
base = pyuvdata.miriad.Miriad()
base.read_miriad(args[0])
#except:
#    pass
print args[1]
#base.read_miriad(args[0])
#binned = uvdata.miriad.Miriad()
#binned.read_miriad(args[1],run_check=False)
#except:
#    pass

def haar(points, a):
    haarf = n.zeros(points)
    haarf[0:a/2] = 1.
    haarf[a/2:a] = -1.
    return haarf

def aipyImport(mirFile,ants,pol):
    uv = a.miriad.UV(mirFile)
    a.scripting.uv_selector(uv,ants,pol)
    dayData = []
    dayTimes = []
    for p,d,f in uv.all(raw=True):
        dayData.append(d)
        dayTimes.append(uv['lst'])
    return n.array(dayData),n.array(dayTimes)


def featArray(data,timesObs):
    freqs1 = n.linspace(100,200,data.shape[1])
    sh = n.shape(data)
    CW1mean = n.zeros_like(data)
    for i in range(CW1mean.shape[1]):
        CW1 = cwt(n.abs(data[:,i]),haar,n.arange(1,10,1))
        CW1 = n.ma.masked_where(CW1==0,CW1)
        CW1mean[:,i] = n.ma.mean(n.abs(CW1),0)

    CT1mean = n.zeros_like(data)
    for j in range(CW1mean.shape[0]):
        CT1 = cwt(n.abs(data[j,:]),signal.morlet,n.arange(1,3,1))
        CT1 = n.ma.masked_where(CT1==0,CT1)
        CT1mean[j,:] = n.mean(n.abs(CT1),0)
    X1 = n.zeros((sh[0]*sh[1],5))
    X1[:,0] = (n.real(data)).reshape(sh[0]*sh[1])
    X1[:,1] = (n.imag(data)).reshape(sh[0]*sh[1])
    #X1[:,2] = n.abs(CW1mean).reshape(sh[0]*sh[1])
    X1[:,2] = n.log10(n.abs(CT1mean)).reshape(sh[0]*sh[1])
    X1[:,3] = (n.array([timesObs]*sh[1])).reshape(sh[0]*sh[1])
    X1[:,4] = (n.array([freqs1]*sh[0])).reshape(sh[0]*sh[1])
    X1 = n.nan_to_num(X1)
    for m in range(n.shape(X1)[1]):
        X1[:,m] = X1[:,m]/n.abs(X1[:,m]).max()
    X1[n.abs(X1)>10**100] = 0
    print X1.max(),X1.min()
#    X1 = normalize(X1,norm='l2',axis=1)
    X1 = n.nan_to_num(X1)
    return X1


def featTargetArray(data,timesObs,maskedLST):
    freqs1 = n.linspace(100,200,data.shape[1])
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
    X1 = n.zeros((sh[0]*sh[1],6))
    X1[:,0] = (n.real(data)).reshape(sh[0]*sh[1])
    X1[:,1] = (n.imag(data)).reshape(sh[0]*sh[1])
    #X1[:,2] = n.abs(CW1mean).reshape(sh[0]*sh[1])
    X1[:,2] = n.log10(n.abs(CT1mean)).reshape(sh[0]*sh[1])
    X1[:,3] = (n.array([timesObs]*sh[1])).reshape(sh[0]*sh[1])
    X1[:,4] = (n.array([freqs1]*sh[0])).reshape(sh[0]*sh[1])
    X1[:,5] = (maskedLST[0:sh[0],:]).reshape(sh[0]*sh[1])
    X1 = n.nan_to_num(X1)
    for m in range(n.shape(X1)[1]):
        X1[:,m] = X1[:,m]/n.abs(X1[:,m]).max()
    X1[n.abs(X1)>10**100] = 0
#    X1 = normalize(X1,norm='l2',axis=1)
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

for s in sep:
    print s
#    BLS=os.popen("python ~/capo/pspec_pipeline/getbls.py --sep="+s+" -C psa6240_FHD /users/jkerriga/data/jkerriga/PGBH/even/Pzen.2456242.30605.uvcRREcACOTUcHPAB").read()
    print BLS
#    BLS = BLS.split(',')
    for b in BLS:
        print b
        a1,a2 = b.split('_')
        bsl_data = base.baseline_array==base.antnums_to_baseline(a1,a2)
        #bsl_train = binned.baseline_array==binned.antnums_to_baseline(a1,a2)
        data = base.data_array[bsl_data,0,:,0]
        #train_data = binned.data_array[bsl_train,0,:,0]
        #try:
        train_data,train_times = aipyImport(args[1],b,'xx')
        #except:
        #    continue
        if data.shape[0] == 0:
            continue
        if train_data == ():
            continue
        times = base.time_array[bsl_data]
        #train_times = binned.time_array[bsl_train]
        print train_data.shape,train_times.shape
        Xtrain = featArray(train_data,train_times)
        pipo = 1.
        fit_param = 3
        pi = n.mean(n.abs(train_data))
        #while pipo >= 1.:
#        gmm = mixture.BayesianGaussianMixture(n_components=fit_param,covariance_type='full',n_init=10,max_iter=10000).fit(Xtrain)
#        Nlabel = gmm.predict(Xtrain)
#        sh = n.shape(train_data)
#        Nlabel_,maxlabel = LabelMaker(Nlabel,sh)

#        po = n.mean(n.abs(train_data[Nlabel_==maxlabel]))
        #fit_param += 1
#        pipo = po/pi
         #   print pipo,fit_param
#        maskedLST = n.zeros(n.shape(train_data))
#        maskedLST[Nlabel_ == maxlabel] = 1.
        Xarr = featArray(data,times)
        gmm2 = mixture.GaussianMixture(n_components=fit_param,covariance_type='full',n_init=10,max_iter=10000).fit(Xarr)
        Nlabel_target = gmm2.predict(Xarr)
        print n.sum(n.unique(Nlabel_target))
        sh_targ = n.shape(data)
        Nlabel_target,maxlabel_target = LabelMaker(Nlabel_target,sh_targ)
        MASK = n.zeros(n.shape(data)).astype(bool)
        MASK[Nlabel_target!=maxlabel_target] = True
        base.flag_array[bsl_data,0,:,0] += MASK
base.write_miriad(args[0]+'r')
