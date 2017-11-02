import numpy as n
from glob import glob
import pylab as pl
from numpy import random

pspecs = glob('/users/jkerriga/brownscripts/jrk/PGProcessedUF_RatioWideFRF/95_115/I/*/pspec_boot*npz')
print pspecs
NBOOT = 300
pk = []
for i in pspecs:
    x = n.load(i)
    pk.append(x['nocov_vs_t'])
#x = n.load(pspecs[0])
pk = n.array(n.abs(pk))
ntimes = pk.shape[-1]
pk_boots = []
print pk.shape
pkrs = pk.reshape(pk.shape[0]*pk.shape[2],pk.shape[1])
print pkrs.shape
pl.hist(pkrs[:,2],10000)
num,bins = n.histogram(pkrs[:,2],10000)
pl.axvline(n.median(pkrs[:,3]),0,n.max(num),color='r',label='median')
pl.axvline(n.mean(pkrs[:,3]),0,n.max(num),color='g',label='mean')
pl.xlim(0,2.0)
pl.ylim(0,n.max(num)+0.2*n.max(num))
pl.legend()
pl.show()
for boot in xrange(NBOOT):
    for t in xrange(ntimes):
        t = random.choice(range(pk.shape[-1]))
        b = random.choice(range(pk.shape[0]))
        pk_boots.append(pk[b,:,t])
pk_boots = n.array(pk_boots)
print pk_boots.shape
pk_boots = n.sort(pk_boots,axis=0)
up_thresh = int(n.around(0.975 * pk_boots.shape[0])) #2 sigma, single tail
print up_thresh
#pl.hist(pk_boots[:,10],1000)
#pl.show()

pk = n.mean(n.abs(pk),axis=0)
pk = n.mean(pk,axis=1)
print pk.shape
print pk
err = n.std(pk_boots,axis=0)#pk-1.96*n.std(pk_boots,axis=0)/n.sqrt(pk_boots.shape[0])#
err = (pk_boots[up_thresh,:] - pk) / 2
print err.shape
k = x['kpl']
#pk = n.mean(pk,axis=0)
#err = n.std(pk,axis=0)
#print err.shape
#pk = n.mean(pk,axis=1)
#pk = n.abs(pk)
#err = n.mean(err,axis=1)
fig = pl.figure()
ax = fig.add_subplot(111)
pkfold = n.zeros(11)
kfold = n.zeros(11)
errfold = n.zeros(11)
pkfold[1:11] = n.mean((pk[::-1][1:11],pk[1:11]),axis=0)
kfold[1:11] = n.mean((k[::-1][1:11],k[1:11]),axis=0)

errfold[1:11] = n.mean((err[::-1][1:11],err[1:11]),axis=0)
ax.errorbar(k,pk,yerr=err,fmt='+')
#pl.plot(k,pk)
pl.hlines(1,-0.5,0.5,linestyles='--')
pl.ylim(0,1.05)
pl.show()

