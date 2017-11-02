import numpy as n
import uvdata
import pylab as pl
import optparse, sys, os
import aipy as AP
from glob import glob
import capo as C
from astropy.utils.data import clear_download_cache
clear_download_cache()
o = optparse.OptionParser()
o.add_option("--dirty", action='store', dest='dirty')
o.add_option("--residual", action='store', dest='residual')
opts,args = o.parse_args(sys.argv[1:])
chans=20
d = n.zeros((2,chans),dtype=complex)
#model = n.zeros((1,chans),dtype=complex)#[]
#modelF = n.zeros((1,chans),dtype=complex) #[]
model = {}
modelF = {}
#blk = AP.dsp.gen_window(203,window='blackman-harris')
blk = n.ones(20)
dirtlist = glob(opts.dirty)
reslist = glob(opts.residual)
uv = AP.miriad.UV(dirtlist[0])
aa = AP.cal.get_aa('psa6240_FHD', uv['sdf'], uv['sfreq'], uv['nchan'])
filters = C.dspec.wedge_width_by_bl(aa, uv['sdf'], chans, offset=15.0)
ch = ['30_50','95_115','115_135','135_155']
for c in ch:
    model[c] = n.zeros((1,chans),dtype=complex)
    modelF[c] = n.zeros((1,chans),dtype=complex)
#bins = n.add(filters[(41,49)],chans/2)
#bin = n.zeros(203)
#uthresh,lthresh = filters[(41,49)]
#bin[uthresh] = 1.
#bin[lthresh] = 1.
del(uv)
print dirtlist,reslist
for i in dirtlist:
    mir = uvdata.miriad.Miriad()
    try:
        mir.read_miriad(i)
    except:
        pass
    bsl = mir.antnums_to_baseline(41,49)==mir.baseline_array
    for c in ch:
        chan = c.split('_')
        intmodel = mir.data_array[bsl,0,int(chan[0]):int(chan[1]),0]
    #intmodel = n.mean(intmodel,0)
    #model = n.mean((model,intmodel),0)
        model[c] = n.vstack((model[c],n.abs(intmodel)))
    del(bsl)
    del(mir)

for j in reslist:
    print j
    mir = uvdata.miriad.Miriad()
    try:
        mir.read_miriad(j)
    except:
        pass
    bsl = mir.antnums_to_baseline(41,49)==mir.baseline_array
    for c in ch:
        chan = c.split('_')
        intmodelF = mir.data_array[bsl,0,int(chan[0]):int(chan[1]),0]
    #intmodelF = n.mean(intmodelF,0)
    #modelF = n.mean((modelF,intmodelF),0)
        modelF[c] = n.vstack((modelF[c],n.abs(intmodelF)))
    del(bsl)
    del(mir)

### Find and average the PSPECs 
dir = 'PGProcessedUF_RatioWideFRF'
npz = glob('/users/jkerriga/brownscripts/jrk/'+dir+'/*/I/*.npz')
pspecs = {}
pk = []
for i in npz:
    print i
    x = n.load(i)
    #pspecs[chans+'_pk'] = x['pk']
    pk.append(n.mean(x['pk']))
print len(pk)
visPrat = []
for c in ch:
    visPrat.append(n.mean(modelF[c])/n.mean(model[c]))
print len(visPrat)
freqs = n.linspace(100*10**6,200*10**6,203)[95:115]
delays = n.fft.fftfreq(freqs.size,freqs[1]-freqs[0])
delays = n.fft.fftshift(delays)
delaymax = 30.0/(2.99*10**8)

#for c in ch:
#    mu_model = n.mean(n.abs(model[c]),0)
#    mu_modelF = n.mean(n.abs(modelF[c]),0)
#    pl.semilogy(delays,mu_model,'b')
#    pl.semilogy(delays,mu_modelF,'r')
pl.plot(pk[::-1],visPrat)
pl.show()

