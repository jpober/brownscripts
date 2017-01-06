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
model = n.zeros((1,chans),dtype=complex)#[]
modelF = n.zeros((1,chans),dtype=complex) #[]
#model = []
#modelF = []
#blk = AP.dsp.gen_window(203,window='blackman-harris')
blk = n.ones(20)
dirtlist = glob(opts.dirty)
reslist = glob(opts.residual)
uv = AP.miriad.UV(dirtlist[0])
aa = AP.cal.get_aa('psa6240_FHD', uv['sdf'], uv['sfreq'], uv['nchan'])
filters = C.dspec.wedge_width_by_bl(aa, uv['sdf'], chans, offset=15.0)
bins = n.add(filters[(41,49)],chans/2)
bin = n.zeros(203)
uthresh,lthresh = filters[(41,49)]
bin[uthresh] = 1.
bin[lthresh] = 1.
del(uv)
print dirtlist,reslist
for i in dirtlist:
    print i
    mir = uvdata.miriad.Miriad()
    try:
        mir.read_miriad(i)
    except:
        pass
    bsl = mir.antnums_to_baseline(41,49)==mir.baseline_array
    intmodel = n.fft.fftshift(n.fft.fft(mir.data_array[bsl,0,95:115,0]*blk))
    #intmodel = n.mean(intmodel,0)
    #model = n.mean((model,intmodel),0)
    model = n.vstack((model,n.abs(intmodel)))
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
    intmodelF = n.fft.fftshift(n.fft.fft(mir.data_array[bsl,0,95:115,0]*blk))
    #intmodelF = n.mean(intmodelF,0)
    #modelF = n.mean((modelF,intmodelF),0)
    modelF = n.vstack((modelF,n.abs(intmodelF)))
    del(bsl)
    del(mir)


freqs = n.linspace(100*10**6,200*10**6,203)[95:115]
delays = n.fft.fftfreq(freqs.size,freqs[1]-freqs[0])
delays = n.fft.fftshift(delays)
delaymax = 30.0/(2.99*10**8)
print model.shape

mu_model = n.mean(n.abs(model),0)
mu_modelF = n.mean(n.abs(modelF),0)
#mu_model = #n.abs(model[0,:])
#mu_modelF = #n.abs(modelF[0,:])
pl.semilogy(delays,mu_model,'b')
pl.semilogy(delays,mu_modelF,'r')
#pl.vlines(delaymax,10**-7,10**6)

#pl.vlines(delaymax+15*10**-9,10**-8,10**6,linestyles='dashed',color='r')
#pl.vlines(delays[bins[0]],10**-8,10**6,linestyles='dashed')
#pl.vlines(-1*delaymax,10**-7,10**6)
#pl.vlines(-1*delaymax-15*10**-9,10**-7,10**6,linestyles='dashed',color='r')
#pl.vlines(delays[bins[1]],10**-8,10**6,linestyles='dashed')
#pl.plot(20*n.log10(mu_model/n.max(mu_model)),'b')
#pl.plot(20*n.log10(mu_modelF/n.max(mu_modelF)),'r')
pl.show()

