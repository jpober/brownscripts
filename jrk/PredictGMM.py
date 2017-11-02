import os
import numpy as n
from scipy import signal
from sklearn import mixture
from sklearn import cluster
from scipy.signal import cwt
from scipy.stats import skewtest,kurtosistest
from sklearn.preprocessing import normalize
from sklearn.decomposition import PCA
from glob import glob
import pyuvdata
import pylab as pl
import optparse, sys, os
import aipy as a
from sklearn.externals import joblib

def loadFullDay(day):
#    HERAlist = glob('/Users/josh/Desktop/HERA/data/zen.2457458.*.xx.HH.uvcUA')
    HERAlist = glob('/users/jkerriga/data/jkerriga/PaperData/Pzen.2456'+str(day)+'.*.*SPA')
    HERAlist = n.sort(HERAlist)
    HERAdata = []
    times = []
    for k in HERAlist:
        uvHERA = a.miriad.UV(k)
        a.scripting.uv_selector(uvHERA, '41_49', 'I')
        for p,d,f in uvHERA.all(raw=True):
            HERAdata.append(d)
            times.append(uvHERA['lst'])
    HERAdata = n.array(HERAdata)
    times = n.array(times)
    return HERAdata,times

def localStats(data,k,l):
    samples = []
    for p in range(100):
        i = n.random.randint(-1,1)
        j = n.random.randint(-1,1)
        try:
            samples.append(n.abs(data[k+i,l+j]))
        except:
            pass
    return n.var(samples)

def injectFRB(data):
    data = n.array(data)
    for i in range(50,1,-1):
        j = 100./i #*n.log10(1000.*i)
        j = int(j)
        pos = 80#n.random.randint(1,90)
        data[i,j+pos] = data[i,j+pos] + 100000.*(n.random.rand()+1j*n.random.rand())/i
    return data

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

def featArrayPredict(data,times,NNarr):
    sh = n.shape(data)
    freqs = n.linspace(100,200,sh[1])
    X1 = n.zeros((sh[0]*sh[1],4))
    X1[:,0] = n.real(data).reshape(sh[0]*sh[1])
    X1[:,1] = n.imag(data).reshape(sh[0]*sh[1])
    X1[:,2] = NNarr #n.log10(n.abs(NNvar)).reshape(sh[0]*sh[1])
    X1[:,3] = (n.array([freqs]*sh[0])).reshape(sh[0]*sh[1])
    #X1[:,4] = (n.array([times]*sh[1])).reshape(sh[0]*sh[1])
    X1[n.abs(X1)>10**100] = 0
    for m in range(X1.shape[1]):
        X1[:,m] = X1[:,m]/n.abs(X1[:,m]).max()
    X1 = n.nan_to_num(X1)
    return X1

def normalize(X):
    normX = (X-n.mean(X))/n.std(X)
    return normX

def doStatTests(X,labels,ml):
    pca = PCA(n_components=1)
    pca.fit(X[labels==ml])
    XpcaML = pca.transform(X[labels==ml])
    labelsOut = labels
    normXpcaML = (XpcaML-n.mean(XpcaML))/n.std(XpcaML)
    #maxKurt = kurtosistest(normXpcaML)[1]
    #maxSkew = skewtest(normXpcaML)[1]
    for i in n.unique(labels):
        if len(X[labels==i])==0:
            continue
        else:
            Xpca = pca.transform(X[labels==i])
            Xpca = (Xpca-n.mean(Xpca))/n.std(Xpca)
            if len(Xpca) < 9:
                labelsOut[labels==i] = -1
                continue
            if False:
                if len(Xpca) < 9:
                    labelsOut[labels==i] = -1
                    continue
                pl.figure()
                if skewtest(Xpca)[1] > 0.5 or kurtosistest(Xpca)[1] > 0.5:
                    tag = 'RFI'
                else:
                    tag = 'Not RFI'
                sk = skewtest(Xpca)[1]
                kt = kurtosistest(Xpca)[1]
                sk1 = skewtest(XpcaML)[1]
                kt1 = kurtosistest(XpcaML)[1]
                pl.subplot(211)
                pl.hist(Xpca,50,label=tag+':'+str(sk)+':'+str(kt))
                pl.legend()
                pl.subplot(212)
                pl.hist(XpcaML,50,label=tag+':'+str(sk1)+':'+str(kt1))
                pl.legend()
                pl.show()
            if i == ml:
                continue
            if skewtest(Xpca)[1] > 0.01: #or kurtosistest(Xpca)[1] > 1.:
                labelsOut[labels==i] = -1
            #else:
            #    labelsOut[labels==i] = ml
    return labelsOut

def LabelMaker(newlabels):
    ml1 = 0
    ml1num = 0
    for t in n.unique(newlabels):
        if (newlabels==t).sum()>ml1num:
            ml1 = t
            ml1num = (newlabels==t).sum()
    newlabels[newlabels==ml1] = -1
    newlabels[newlabels!=-1] += 10
    return newlabels

def findMaxLabel(labels):
    mlabel = ()
    maxCt = 0
    for i in n.unique(labels):
        if (labels==i).sum() > maxCt:
            maxCt = (labels==i).sum()
            mlabel = i
        else:
            continue
    return mlabel

def findMinLabel(labels):
    mlabel = ()
    minCt = 10**100
    for i in n.unique(labels):
        if (labels==i).sum() < minCt:
            minCt = (labels==i).sum()
            mlabel = i
        else:
            continue
    return mlabel



o = optparse.OptionParser()
opts,obs = o.parse_args(sys.argv[1:])

dpgmm = joblib.load('HERAPersistDay.pkl')
CL1 = []
CL2 = []
CL3 = []

for o in obs:
#for o in range(242,365):
    print o
    #try:
    uv = pyuvdata.miriad.Miriad()
    uv.read_miriad(o)
    #data,times = loadFullDay(o)
#    print data.shape
    #except:
    #    pass
    for b in n.unique(uv.baseline_array):
        idx = uv.baseline_array==b
        #idx = uv.baseline_array==uv.antnums_to_baseline(41,49)
        if b == n.unique(uv.baseline_array)[0]:
            data = uv.data_array[idx,0,:,0]
        else:
            data = (data+uv.data_array[idx,0,:,0])/2.0
    #inject = n.random.randint(1,5)
#    print inject
#    if inject == 3:
#    data = injectFRB(data)
    sh = n.shape(data)
#        times = uv.lst_array[idx]
    X = featArray(data)
    labels = dpgmm.predict(X)
    labels = labels.reshape(sh[0],sh[1])
    if (labels==1000).sum() > 0:
        print 'Found and FRB!!!'
        print o
    pl.imshow(n.log10(n.abs(data)),aspect='auto')
    pl.show()
#        break
#        ml = findMaxLabel(labels)
#        mask = n.zeros_like(data).astype(bool)
#        mask[labels!=ml] = True
#        uv.flag_array[idx,0,:,0] = mask
    del(labels)
    del(X)
#    uv.write_miriad(o+'r')
#    del(uv)


