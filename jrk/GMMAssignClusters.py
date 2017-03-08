#import matplotlib.pyplot as plt
import pylab as plt
from sklearn import mixture
import aipy as a
import numpy as n
from glob import glob
from sklearn.externals import joblib
def loadAipyData():
    HERAlist = glob('/Users/josh/Desktop/HERA/data/zen.2457458.40*.xx.HH.uvcUA')
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

def featArray(data):
    sh = n.shape(data)
    freqs = n.linspace(100,200,sh[1])
    NNvar = n.zeros_like(data)
    dvar = n.var(n.abs(data))
    for i in range(sh[0]):
        for j in range(sh[1]):
            samples = []
            for p in range(500):
                k = n.random.randint(-1,1)
                l = n.random.randint(-1,1)
                try:
                    samples.append(n.abs(data[k+i,l+j]))
                except:
                    pass
            NNvar[i,j] = n.var(samples)
    X1 = n.zeros((sh[0]*sh[1],4))
    X1[:,0] = n.real(data).reshape(sh[0]*sh[1])
    X1[:,1] = n.imag(data).reshape(sh[0]*sh[1])
    X1[:,2] = (n.log10(n.abs(NNvar)) - n.median(n.log10(n.abs(NNvar)))).reshape(sh[0]*sh[1])
    X1[:,3] = (n.array([freqs]*sh[0])).reshape(sh[0]*sh[1])
    #X1[:,4] = (n.array([times]*sh[1])).reshape(sh[0]*sh[1])
    X1[n.abs(X1)>10**100] = 0
    for m in range(X1.shape[1]):
        X1[:,m] = X1[:,m]/n.abs(X1[:,m]).max()
    X1 = n.nan_to_num(X1)
    return X1

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

data,times = loadAipyData()
print 'Data loaded. Preparing feature array...'
X = featArray(data)
print 'X feature array loaded.'
dpgmm = mixture.BayesianGaussianMixture(n_components=10,covariance_type='full',n_init=1,max_iter=1000,init_params='kmeans',weight_concentration_prior_type='dirichlet_process').fit(X)
labels = dpgmm.predict(X)
sh = n.shape(data)
ml = findMaxLabel(labels)
labels = labels.reshape(sh[0],sh[1])
mask = n.ones_like(data)
#mask[labels!=ml] = 0
fig = plt.figure()
plt.ion()

options={'y':-1,
    'n':100}
labelsC = n.copy(labels)
RFIlabels = []
for c in n.unique(labels):
    instaMask = n.ones_like(data)
    instaMask[labels != c] = 0 
    plt.imshow(n.log10(n.abs(data*instaMask)),aspect='auto',cmap='jet')
    plt.show()
    text=raw_input(str(options))
    #if text in options.keys():
    labelsC[labels==c] = options[text]
    if text == 'y':
        RFIlabels.append(c)
    plt.clf()


mask[labelsC!=100] = 0
fig2 = plt.figure()
plt.imshow(n.log10(n.abs(data*mask)),aspect='auto',cmap='jet')
plt.pause(1000000)
plt.show()
plt.clf()

n.savetxt('RFIlabels.txt',RFIlabels)
joblib.dump(dpgmm,'HERAPersistDay.pkl')
