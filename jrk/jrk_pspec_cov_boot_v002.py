#! /usr/bin/env python

import aipy as a, numpy as n
from matplotlib import pylab as p
import capo as C
import sys, optparse, re, os, random

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True, cal=True)
o.add_option('--plot', action='store_true',
    help='Generate plots')
o.add_option('--identity', action='store_true', default=False,
    help='Use the identity matrix instead of C^-1')
opts,args = o.parse_args(sys.argv[1:])

NBOOT = 300
MEDIAN = True
CLIP = False #True #trim time axis
LO,HI = 40,600
PLOT = opts.plot

#Save npz file info
pk_vs_t = {}
err_vs_t = {}
temp_noise_var = {}
nocov_vs_t = {}
chans = []
afreqs = []

for filename in args:
    print 'Reading', filename
    f = n.load(filename)
    kpl,cmd = f['kpl'], f['cmd']
    path = os.path.dirname(filename)
    if not pk_vs_t.has_key(path):
        pk_vs_t[path] = []
        err_vs_t[path] = []
        temp_noise_var[path] = []
        nocov_vs_t[path] = []
    pk_vs_t[path].append(f['pk_vs_t'])
    scalar = f['scalar']
    afreqs = f['afreqs']
    chans = f['chans']
    nocov_vs_t[path].append(f['nocov_vs_t'])

paths = pk_vs_t.keys() #different path depending on what sep (bltype) the boots are from 
k0 = n.abs(kpl).argmin() #where k=0 is

pk_2d = n.array([pk_vs_t[path] for path in paths]) #(bltype,bootstraps,kpls,times)
nocov_2d = n.array([nocov_vs_t[path] for path in paths]) #(bltype,bootstraps,kpls,times), T averaged over all bls
avg_pk_2d = n.average(pk_2d, axis=1) #(bltype,kpls,times), best estimate of pk in each integration
wgts = n.ones_like(avg_pk_2d)

if opts.identity: #override power spectrum with the version w/o covariance diagonalization
    print 'Overriding power spectrum with non-covariance diagonalized version'
    pk_2d = nocov_2d
    pk_2d = n.abs(pk_2d)
if CLIP: #trim time axis
    pk_2d = pk_2d[...,LO:HI]
    avg_pk_2d = avg_pk_2d[...,LO:HI]
    wgts = wgts[...,LO:HI]
else:
    pass

print n.abs(pk_2d).shape
#print 'Average:',n.max(n.mean(pk_2d))
newpk_2d = n.mean(n.mean(n.mean(pk_2d,axis=0),axis=1),axis=1)
print newpk_2d.shape
pk_2d = pk_2d.transpose([1,2,3,0]).copy() #(bootstraps,kpls,times,bltype)
pk_2d.shape = pk_2d.shape[:-2] + (pk_2d.shape[-2] * pk_2d.shape[-1],) #(bootstraps,kpls,timebls)
wgts = wgts.transpose([1,2,0]).copy() #(kpls,times,bltypes)
wgts.shape = wgts.shape[:-2] + (wgts.shape[-2] * wgts.shape[-1],) #(kpls,timebls)

npaths = pk_2d.shape[0]
ntimes = pk_2d.shape[-1]

pk_boot = []
pk_fold_boot = []

print 'Bootstrapping over boot files and times...'
for boot in xrange(NBOOT):
    if boot % 10 == 0: print '   ',boot,'/',NBOOT
    dsum,dwgt = [],[]
    for t in xrange(ntimes):
        t = random.choice(range(pk_2d.shape[-1])) #random time
        b = random.choice(range(pk_2d.shape[0])) #random boot npz file
        dsum += [pk_2d[b,:,t]]# * wgts[:,t]] #!!!!!!!!!removed wgts
        #dsum.append(pk_2d[b,:,t])
        dwgt += [wgts[:,t]]
    if MEDIAN: dsum,dwgt = n.median(dsum, axis=0), n.median(dwgt, axis=0) #median over time
    else: dsum,dwgt = n.average(dsum,axis=0), n.average(dwgt,axis=0)
    pk_boot.append(dsum)#/dwgt)
    dsum_fold = dsum[k0:].copy()
    dwgt_fold = dwgt[k0:].copy()
    dsum_fold[1:] = 0
    dwgt_fold[1:] = 0
    dsum_pos,dwgt_pos = dsum[k0+1:].copy(), dwgt[k0+1:].copy() 
    dsum_neg,dwgt_neg = dsum[k0-1::-1].copy(), dwgt[k0-1::-1].copy() #for odd # of channels    
    for h in xrange(2): #bootstrap over which half of the spectrum (or both) are used
        h = random.randint(0,1)
        dsum_fold[1:] += [dsum_pos, dsum_neg][h]
        dwgt_fold[1:] += [dwgt_pos, dwgt_neg][h]
    pk_fold_boot.append(dsum_fold)# / dwgt_fold)# !!!!!!removes wgts
pk_boot = n.array(pk_boot).T
pk_fold_boot = n.array(n.abs(pk_fold_boot)).T
print 'Sorting bootstraps...'
pk = n.abs(n.mean(pk_boot, axis=1)) #average over all boots
pk_fold = n.abs(n.mean(pk_fold_boot, axis=1))
#this is excluding imag component in noise estimate `
#pk_boot = n.sort(pk_boot.real, axis=1) #dropping imag component here
#pk_fold_boot = n.sort(pk_fold_boot.real, axis=1) #dropping imag component here
pk_boot = n.sort(pk_boot, axis=1) #dropping imag component here
pk_fold_boot = n.sort(pk_fold_boot, axis=1) 
if True:
    print 'Deriving errors from histogram...'
    up_thresh = int(n.around(0.975 * pk_boot.shape[1])) #2 sigma, single tail
    #important to only include real component in estimation of error
    err = (pk_boot[:,up_thresh] - pk) / 2 #effective "1 sigma" derived from actual 2 sigma
    up_thresh_fold = int(n.around(0.975 * pk_fold_boot.shape[1])) #2 sigma, single tail
    err_fold = (pk_fold_boot[:,up_thresh_fold] - pk_fold) / 2 #effective "1 sigma" derived from actual 2 sigma

    print 'Plotting Pk...'
    print '%4s %11s %18s' % ('k','Pk+err','Pk +/- err')
    for k_,pk_,err_ in zip(kpl, pk.real/1e6, 2*err/1e6):
        print '%6.3f, %9.5f (%9.5f +/- %9.5f)' % (k_, pk_+err_,pk_,err_)
    print '-'*70
    print 'Plotting Pk_fold...'
    print '%4s %11s %18s' % ('k','Pk_fold+err','Pk_fold +/- err')
    for k_,pk_,err_ in zip(kpl[k0:], pk_fold.real/1e6, 2*err_fold/1e6):
        print '%6.3f, %9.5f (%9.5f +/- %9.5f)' % (k_, pk_+err_,pk_,err_)
else:
    err = n.std(pk_boot, axis=1)
    err_fold = n.std(pk_fold_boot, axis=1)

print 'Writing pspec.npz...'
n.savez('pspec.npz', kpl=kpl, pk=pk, err=err, pk_fold=pk_fold, err_fold=err_fold, cmd=cmd,afreqs=afreqs,chans=chans)
    

