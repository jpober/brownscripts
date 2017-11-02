import numpy as n
import pyuvdata
import pylab as pl
import optparse, sys, os
import aipy as a
from glob import glob
import capo as C
import scipy.stats
from astropy.utils.data import clear_download_cache
from scipy.signal import medfilt
clear_download_cache()
c1 = 85
c2 = 125
o = optparse.OptionParser()
o.add_option("--dirty", action='store', dest='dirty')
o.add_option("--residual", action='store', dest='residual')
opts,args = o.parse_args(sys.argv[1:])


def AipyImport(fileList):
    PAPERdata = []
    uv = a.miriad.UV(fileList[0])
    aa = a.cal.get_aa('psa6240_FHD', uv['sdf'], uv['sfreq'], chans)
    filters = C.dspec.wedge_width_by_bl(aa, uv['sdf'], chans, offset=0.0)
    #ants = uv.antenna_names
    #times = []
    #flags = []
    ants = [(41,49)]
    for files in fileList:
        for an in ants:#filters.keys():
            print an
        #for an1 in ants:
        #    for an2 in ants:
            uvPAPER = a.miriad.UV(files)
            a.scripting.uv_selector(uvPAPER, str(an[0])+'_'+str(an[1]),'I')
            for p,d,f in uvPAPER.all(raw=True):
                PAPERdata.append(n.abs(n.fft.fftshift(n.fft.fft(d[c1:c2]*blk))))
            #times.append(uvPAPER['lst'])
            #flags.append(f)
    PAPERdata = n.array(PAPERdata)
    #times = n.array(times)
    #flags = n.array(flags)
    return PAPERdata

chans=c2-c1
d = n.zeros((2,chans),dtype=complex)
model = n.zeros((1,chans),dtype=complex)#[]
modelF = n.zeros((1,chans),dtype=complex) #[]
#model = []
#modelF = []
blk = a.dsp.gen_window(c2-c1,window='blackman-harris')
#blk = n.ones(c2-c1)
dirtlist = glob(opts.dirty)
reslist = glob(opts.residual)
print dirtlist
uv = a.miriad.UV(dirtlist[0])
aa = a.cal.get_aa('psa6240_FHD', uv['sdf'], uv['sfreq'], uv['nchan'])
filters = C.dspec.wedge_width_by_bl(aa, uv['sdf'], chans, offset=15.0)
bins = n.add(filters[(41,49)],chans/2)
bin = n.zeros(c2-c1)
uthresh,lthresh = filters[(41,49)]
bin[uthresh] = 1.
bin[lthresh] = 1.
bin_ = n.fft.fftshift(bin)
del(uv)
print dirtlist,reslist

dirtd = AipyImport(dirtlist)
resd = AipyImport(reslist)

freqs = n.linspace(100*10**6,200*10**6,203)[c1:c2]
delays = n.fft.fftfreq(freqs.size,freqs[1]-freqs[0])
delays = n.fft.fftshift(delays)
delaymax = 30.0/(2.99*10**8)

err = n.std(n.abs(resd)/(dirtd),0)/n.sqrt(len(resd))
mu_model = n.mean(n.abs(dirtd),0)
mu_modelF = n.mean(n.abs(resd),0)
dist = n.abs(resd)/n.abs(dirtd)
cleaned_dist = dist.reshape(-1,chans)
#cleaned_dist = dist[n.isnan(dist)==False].reshape(-1,chans)
#cleaned_dist = cleaned_dist[n.isinf(cleaned_dist)==False].reshape(-1,chans)

cleaneddist = []
#idxs = n.argwhere(cleaned_dist[:,15]-n.mean(cleaned_dist[:,15]) < 2.*n.std(cleaned_dist[:,15]))

#print n.shape(idxs),'index'
print n.shape(cleaned_dist),'clean dist'
#for i in range(21):
#    cleaneddist.append(cleaned_dist[cleaned_dist[:,i]-n.mean(cleaned_dist[:,i]) < 3.*n.std(cleaned_dist[:,i]),i]) #medfilt(cleaned_dist) ### each kmode has a different number of good measurements
#cleaneddist = n.array(cleaneddist).reshape(-1,21)
#cleaneddist = cleaned_dist[idxs,:].reshape(-1,21)

rat = n.median(cleaned_dist,0)#n.array(scipy.stats.mode(cleaned_dist,0))[0,0,:]
print n.shape(rat)
rat_fold = []
delay_fold = []
mid = (chans-1)/2
#for i in range(mid):
#    rat_fold.append(n.mean((rat[mid + i],rat[mid - i])))
#    delay_fold.append(n.mean((delays[mid+i],delays[mid-i])))
print n.shape(rat),'rat'
print n.shape(delays),'delays'
rat_err = n.std(n.abs(resd)/n.abs(dirtd),0)
#mu_model = #n.abs(model[0,:])
#mu_modelF = #n.abs(modelF[0,:])
#pl.semilogy(delays,mu_model,'b.')
#pl.semilogy(delays,mu_modelF,'r.')
#pl.ylim(0,1)
print n.shape(delays),n.shape(dirtd)
#pl.errorbar(delays,rat,yerr=rat_err)
pl.semilogy(delays,mu_model)
pl.semilogy(delays,mu_modelF)
#pl.vlines(delays[uthresh+chans/2],10**-1,10**0)
#pl.vlines(delays[lthresh+chans/2],10**-1,10**0)
pl.vlines(delays[bin_==1],10**-1,10**0)

#pl.vlines(delaymax+15*10**-9,10**-8,10**6,linestyles='dashed',color='r')
#pl.vlines(delays[bins[0]],10**-8,10**6,linestyles='dashed')
#pl.vlines(-1*delaymax,10**-7,10**6)
#pl.vlines(-1*delaymax-15*10**-9,10**-7,10**6,linestyles='dashed',color='r')
#pl.vlines(delays[bins[1]],10**-8,10**6,linestyles='dashed')
#pl.plot(20*n.log10(mu_model/n.max(mu_model)),'b')
#pl.plot(20*n.log10(mu_modelF/n.max(mu_modelF)),'r')
pl.show()


#dist[n.isnan(dist)] = 0
#dist[n.isinf(dist)] = 0
for i in range(chans):
    print '% below 1: ',(cleaned_dist[:,i]<=1.).sum()/(1.*len(cleaned_dist[:,i]))
pl.figure()
mode = n.array(scipy.stats.mode(cleaned_dist[:,14],0))
print mode[0]
for i in range(chans):
    pl.subplot(6,4,i+1)
    cd = cleaned_dist[:,i]
    cd = cd[~n.isnan(cd)]
    cd = cd[~n.isinf(cd)]
    pl.hist(cd,1000)
    pl.axvline(n.median(cd),0,10**3)
    pl.axvline(n.mean(cd),0,10**3)
    pl.axvline(mode[0],0,10**3)
    pl.xlim(0,1.5)
pl.show()
