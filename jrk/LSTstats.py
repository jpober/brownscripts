import numpy as n
import aipy as a
from glob import glob
import pylab as pl




def loadAipyData():
    PAPERlist = glob('/gpfs/data/jpober/jkerriga/PaperData/odd/Pzen.2456*.*HPABR')
    PAPERlist = n.sort(PAPERlist)
    print PAPERlist
    PAPERdata = []
    times = []
    for k in PAPERlist:
        uvPAPER = a.miriad.UV(k)
        a.scripting.uv_selector(uvPAPER, '41_49', 'I')
        for p,d,f in uvPAPER.all(raw=True):
            if uvPAPER['lst'] > 1. and uvPAPER['lst'] < 1.01:
                PAPERdata.append(d[100])
                times.append(uvPAPER['lst'])
                break
            #print times
    PAPERdata = n.array(PAPERdata)
    times = n.array(times)
    return PAPERdata,times


data,times = loadAipyData()
pl.hist(n.abs(data))
pl.show()
