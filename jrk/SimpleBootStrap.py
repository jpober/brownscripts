import numpy as n
import aipy as a
import pylab as pl
from glob import glob


ch1 = 85
ch2 = 125

def aipyImport(mirFile,ants,pol):
    dayData = []
    dayTimes = []
    for i in mirFile:
        print i
        uv = a.miriad.UV(i)
        a.scripting.uv_selector(uv,ants,pol)
        for p,d,f in uv.all(raw=True):
            dayData.append(d)
            dayTimes.append(uv['lst'])
    return n.array(dayData),n.array(dayTimes)

def delay(data):
    DATA = n.fft.fft(data,axis=1)
    DATA_ = n.fft.fftshift(DATA,axes=1)
    return DATA_

lsts = glob('/users/jkerriga/data/jkerriga/PaperData/even15/sep*/*uvHBA')
lstsS = glob('/users/jkerriga/data/jkerriga/PaperData/even15/sep*/*uvSBA')
data,times = aipyImport(lsts,'41_49','I')
dataS,timesS = aipyImport(lstsS,'41_49','I')
zslice = data[:,ch1:ch2]
zsliceS = dataS[:,ch1:ch2]
win = a.dsp.gen_window(ch2-ch1,window='blackman-harris')

ZSLICE = n.mean(n.abs(delay(zslice*win)),0)
std = n.std(n.abs(delay(zslice*win)),0)
ZSLICES = n.mean(n.abs(delay(zsliceS*win)),0)
stdS = n.std(n.abs(delay(zsliceS*win)),0)

diff = n.abs(delay(zslice*win)) - n.abs(delay(zsliceS*win))
diffstd = n.std(diff)
diffmu = n.mean(diff,0)

x = n.linspace(-20,20,40)
#pl.errorbar(x,ZSLICE,yerr=std,color='red')
#pl.errorbar(x,ZSLICES,yerr=stdS)
pl.errorbar(x,diffmu,yerr=diffstd)
#pl.plot(ZSLICE,'.')
#pl.plot(ZSLICES,'+')
pl.show()
