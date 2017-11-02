import numpy as n
import uvdata
import pylab as pl
import optparse, sys, os
import aipy as a
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
vislist = sorted(glob('/gpfs/data/jpober/jkerriga/ProvingGrounds/odd/Pzen.2456243.*.uvcRREcACOTUcHPA'))
FSlist = sorted(glob('/gpfs/data/jpober/jkerriga/ProvingGrounds/odd/Pzen.2456243.*.uvcRREcACOTUcSPA'))
FAlist = sorted(glob('/gpfs/data/jpober/jkerriga/ProvingGrounds/odd/Pzen.2456243.*.uvcRREcACOTUcSPABR'))
lstBinlist = sorted(glob('/gpfs/data/jpober/jkerriga/PaperData/even15/sep0,1/lst.2456242.*.*.uvSBA'))
frflist = sorted(glob('/gpfs/data/jpober/jkerriga/PaperData/even15/sep0,1/lst.2456242.*.*.uvSBMxA'))


def AipyImport(fileList):
    PAPERdata = []
    times = []
    flags = []
    for f in fileList:
        uvPAPER = a.miriad.UV(f)
        a.scripting.uv_selector(uvPAPER, '41_49','I')
        for p,d,f in uvPAPER.all(raw=True):
            PAPERdata.append(d)
            times.append(uvPAPER['lst'])
            flags.append(f)
    PAPERdata = n.array(PAPERdata)
    times = n.array(times)
    flags = n.array(flags)
    return PAPERdata,(times*180/n.pi)/15,flags

visData,visTimes,visflags = AipyImport(vislist)
FSData,FSTimes,FSflags = AipyImport(FSlist)
FAData,FATimes,FAflags = AipyImport(FAlist)
lstBinData,lstBinTimes,binflags = AipyImport(lstBinlist)
frfData,frfTimes,frfflags = AipyImport(frflist)
FSData = FSData*n.logical_not(FSflags)


f = pl.figure()
#f, axes = pl.subplots(5,1,sharey=True)
f.text(0.06,0.5,'LST (Hrs)',va='center',rotation='vertical')
f.text(0.4,0.13,'Freq. (MHz)',va='center')
#pl.setp(axes.flat,aspect=0.1,adjustable='datalim')
freqs = n.linspace(100,200,203)
f1 = 30
f2 = 160
t1 = 0
t2 = 644
sp = 20
print visTimes
rat = 30#0.1*freqs[f2]/visTimes.max()
ax1 = pl.subplot(151)
pl.imshow(n.log10(n.abs(visData[:,f1:f2])),aspect=rat,interpolation='none',vmin=-1,vmax=4,extent=[freqs[f1],freqs[f2],visTimes.max(),visTimes.min()])
ax1.annotate('A',xy=(120,1),xytext=(120,1),size=20)
ax1.tick_params(axis='x',labelsize=10)
ax1.set_xticklabels(n.linspace(freqs[f1+sp],freqs[f2-sp],3).astype(int))
ax1.set_xticks(n.linspace(freqs[f1+sp],freqs[f2-sp],3))

ax2 = pl.subplot(152)
ax2.set_yticklabels([])
pl.imshow(n.log10(n.abs(FSData[:,f1:f2])),aspect=rat,interpolation='none',vmin=-1,vmax=4,extent=[freqs[f1],freqs[f2],visTimes.max(),visTimes.min()])
ax2.annotate('B',xy=(120,1),xytext=(120,1),size=20)
ax2.set_xticklabels(n.linspace(freqs[f1+sp],freqs[f2-sp],3).astype(int))
ax2.set_xticks(n.linspace(freqs[f1+sp],freqs[f2-sp],3))
ax2.tick_params(axis='x',labelsize=10)

ax3 = pl.subplot(153)
ax3.set_yticklabels([])
pl.imshow(n.log10(n.abs(FAData[:,f1:f2])),aspect=rat,interpolation='none',vmin=-1,vmax=4,extent=[freqs[f1],freqs[f2],visTimes.max(),visTimes.min()])
ax3.annotate('C',xy=(120,1),xytext=(120,1),size=20)
ax3.set_xticklabels(n.linspace(freqs[f1+sp],freqs[f2-sp],3).astype(int))
ax3.set_xticks(n.linspace(freqs[f1+sp],freqs[f2-sp],3))
ax3.tick_params(axis='x',labelsize=10)

ax4 = pl.subplot(154)
ax4.set_yticklabels([])
pl.imshow(n.log10(n.abs(lstBinData[:,f1:f2])),aspect=rat,interpolation='none',vmin=-1,vmax=4,extent=[freqs[f1],freqs[f2],visTimes.max(),visTimes.min()])
ax4.annotate('D',xy=(120,1),xytext=(120,1),size=20)
ax4.set_xticklabels(n.linspace(freqs[f1+sp],freqs[f2-sp],3).astype(int))
ax4.set_xticks(n.linspace(freqs[f1+sp],freqs[f2-sp],3))
ax4.tick_params(axis='x',labelsize=10)

ax5 = pl.subplot(155)
ax5.set_yticklabels([])
pl.imshow(n.log10(n.abs(frfData[:,f1:f2])),aspect=rat,interpolation='none',vmin=-1,vmax=4,extent=[freqs[f1],freqs[f2],visTimes.max(),visTimes.min()])
ax5.annotate('E',xy=(120,1),xytext=(120,1),size=20)
ax5.set_xticklabels(n.linspace(freqs[f1+sp],freqs[f2-sp],3).astype(int))
ax5.set_xticks(n.linspace(freqs[f1+sp],freqs[f2-sp],3))
ax5.tick_params(axis='x',labelsize=10)
#pl.xlim(freqs[f1],freqs[f2])
#pl.xlabel('Freq. (Ghz)')

f.subplots_adjust(right=0.8)
#cbar_ax = f.add_axes([0.85, 0.25, 0.03, 0.5])
cbar_ax = f.add_axes([0.85, 0.19, .03, .62])
pl.colorbar(cax=cbar_ax)
f.text(0.83,0.15,'$log_{10}(Jy)$',va='center')
#pl.tight_layout()
#cbar_ax.xlabel('$log_{10}$(Jy)')
pl.savefig('PSPECFlow.pdf')
pl.show()

