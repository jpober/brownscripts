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

def loadFullDay():
    HERAlist = glob('/Users/josh/Desktop/HERA/data/zen.2457458.*.xx.HH.uvcUA')
    HERAdata = []
    times = []
    for k in HERAlist:
        uvHERA = a.miriad.UV(k)
        a.scripting.uv_selector(uvHERA, '9_10', 'xx')
        for p,d,f in uvHERA.all(raw=True):
            HERAdata.append(d)
            times.append(uvHERA['lst'])
    HERAdata = n.array(HERAdata)
    times = n.array(times)
    return HERAdata,times

def localStats(data,k,l):
    samples = []
    for p in range(100):
        i = n.random.randint(-2,2)
        j = n.random.randint(-2,2)
        try:
            samples.append(n.abs(data[k+i,l+j]))
        except:
            pass
    return 1/n.var(samples)
#    if dvar > n.var(samples):
#        return 1
#    else:
#        return 0

def featArray(data,times):
    sh = n.shape(data)
    freqs = n.linspace(100,200,sh[1])
#    CT1mean = n.zeros_like(data)
    #for j in range(CT1mean.shape[0]):
    #    CT1 = cwt(data[j,:],signal.morlet,n.arange(1,2,1))
    #    CT1mean[j,:] = n.mean(n.abs(CT1),0)
    NNvar = n.zeros_like(data)
    dvar = n.var(n.abs(data))
    for i in range(sh[0]):
        for j in range(sh[1]):
            NNvar[i,j] = localStats(data,i,j)
    X1 = n.zeros((sh[0]*sh[1],5))
    X1[:,0] = n.real(data).reshape(sh[0]*sh[1])
    X1[:,1] = n.imag(data).reshape(sh[0]*sh[1])
    X1[:,2] = n.log10(n.abs(NNvar)).reshape(sh[0]*sh[1])
    X1[:,3] = (n.array([freqs]*sh[0])).reshape(sh[0]*sh[1])
    X1[:,4] = (n.array([times]*sh[1])).reshape(sh[0]*sh[1])
    X1[n.abs(X1)>10**100] = 0
    for m in range(X1.shape[1]):
        X1[:,m] = X1[:,m]/n.abs(X1[:,m]).max()
    X1 = n.nan_to_num(X1)
    #X1 = normalize(X1,norm='l2',axis=0)
    return X1,n.log10(n.abs(NNvar)).reshape(sh[0]*sh[1])

def featArrayPredict(data,times,NNarr):
    sh = n.shape(data)
    freqs = n.linspace(100,200,sh[1])
    X1 = n.zeros((sh[0]*sh[1],5))
    X1[:,0] = n.real(data).reshape(sh[0]*sh[1])
    X1[:,1] = n.imag(data).reshape(sh[0]*sh[1])
    X1[:,2] = NNarr #n.log10(n.abs(NNvar)).reshape(sh[0]*sh[1])
    X1[:,3] = (n.array([freqs]*sh[0])).reshape(sh[0]*sh[1])
    X1[:,4] = (n.array([times]*sh[1])).reshape(sh[0]*sh[1])
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

def Cluster(X,ncomps):
    dpgmm = mixture.BayesianGaussianMixture(n_components=ncomps,covariance_type='full',n_init=20,max_iter=10000,init_params='kmeans',weight_concentration_prior_type='dirichlet_process').fit(X)
    labels = dpgmm.predict(X)
    return labels

def ClusterVanilla(X,ncomps):
    dpgmm = mixture.GaussianMixture(n_components=ncomps,covariance_type='full',n_init=1,max_iter=10000,init_params='random').fit(X)
    labels = dpgmm.predict(X)
    return labels

def kMeans(X,nclust):
    km = cluster.KMeans(n_clusters=nclust).fit(X)
    return km.labels_


o = optparse.OptionParser()
opts,obs = o.parse_args(sys.argv[1:])


uv = pyuvdata.miriad.Miriad()

#fullDay,fullDayTimes = loadFullDay()
#Xfull = featArray(fullDay,fullDayTimes)

#km = cluster.KMeans(n_clusters=5).fit(Xfull)
#labelsfull = km.labels_

#dpgmm = mixture.BayesianGaussianMixture(n_components=3,covariance_type='full',n_init=1,max_iter=10000,init_params='random',weight_concentration_prior_type='dirichlet_distribution').fit(Xfull)
#labelsfull = dpgmm.predict(Xfull)
#ml = findMaxLabel(labelsfull)
CL1 = []
CL2 = []
CL3 = []
for o in obs:
    print o
    #try:
    uv.read_miriad(o)
    #except:
    #    pass
    idx = uv.baseline_array==uv.antnums_to_baseline(9,10)
    data = uv.data_array[idx,0,:,0]
    X,NNArr = featArray(data,uv.lst_array[idx])
    dpgmm = mixture.BayesianGaussianMixture(n_components=3,covariance_type='full',n_init=20,max_iter=10000,init_params='kmeans',weight_concentration_prior_type='dirichlet_process').fit(X)
    labels = dpgmm.predict(X)
    for b in n.unique(uv.baseline_array):
        idx = (b==uv.baseline_array)
        if uv.baseline_to_antnums(b) == (9,10):
            continue
        data = uv.data_array[idx,0,:,0]
        #X,NNArr = featArray(data,uv.lst_array[idx])
        sh = data.shape
#        print n.size(X)
        X = featArrayPredict(data,uv.lst_array[idx],NNArr)
        print 'Predicting...'
        labels = dpgmm.predict(X)
        #labels = kMeans(X,100)
            #labels = Cluster(X,3)
        
        print 'CL1:',(labels==0).sum(),'CL2:',(labels==1).sum(),'CL3:',(labels==2).sum()
        #ml = findMaxLabel(labels)
        CL1.append((labels==0).sum())
        CL2.append((labels==1).sum())
        CL3.append((labels==2).sum())
        
        #reLabel = LabelMaker(labels)
        #else:
        #    labels = dpgmm.predict(X)
        #    reLabel = LabelMaker(labels)
        
#        ml = findMaxLabel(labels)
        #for m in range(1):
        #    ml = findMaxLabel(labels)
        #    labels2 = kMeans(X[labels==ml])
        #    labels[labels==ml] = labels2+3
        
        #ml = findMaxLabel(labels)
        labels = labels.reshape(sh[0],sh[1])
        ml = findMaxLabel(labels)
        mask = n.zeros_like(data).astype(bool)
        mask[labels!=ml] = True
        del(labels)
        del(X)
        ### Try re-mixture modeling with mask as 4th parameter
        #X2 = featArray2(data,mask.astype(float))
        #labels2 = Cluster(X2)
        #ml = findMaxLabel(labels2)
        #labels2 = labels2.reshape(sh[0],sh[1])
        #mask2 = n.zeros_like(data).astype(bool)
        #mask2[labels2!=ml] = True

        uv.flag_array[idx,0,:,0] = mask
        #uv.data_array[idx,0,:,0] = uv.data_array[idx,0,:,0]*n.logical_not(mask)

    uv.write_miriad(o+'r')
pl.hist(CL1,50)
pl.hist(CL2,50)
pl.hist(CL3,50)
pl.show()
