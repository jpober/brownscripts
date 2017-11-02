import numpy as n
import uvdata
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
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
o.add_option("--bsl",action='store',dest='bsl')
o.add_option("--dt", dest='dt', action='store_true')
opts,args = o.parse_args(sys.argv[1:])

#d = n.zeros((2,chans),dtype=complex)
dirtlist = sorted(glob(opts.dirty))
reslist = sorted(glob(opts.residual))
dirty = n.zeros((1,203))
res = n.zeros((1,203))
bslnum = opts.bsl.split('_')
dt = opts.dt

for i in dirtlist:
    print i
    mir = uvdata.miriad.Miriad()
    try:
        mir.read_miriad(i)
    except:
        pass
    bsl = mir.antnums_to_baseline(int(bslnum[0]),int(bslnum[1]))==mir.baseline_array
    dirty = n.vstack((dirty,mir.data_array[bsl][:,0,:,0]))
    del(bsl)
    del(mir)

for j in reslist:
    print j
    mir = uvdata.miriad.Miriad()
    try:
        mir.read_miriad(j)
    except:
        pass
    bsl = mir.antnums_to_baseline(int(bslnum[0]),int(bslnum[1]))==mir.baseline_array
    res = n.vstack((res,mir.data_array[bsl][:,0,:,0]))
    del(bsl)
    del(mir)
dirty = dirty[1:dirty.shape[0],:]
dirty = dirty[:,20:156]
res = res[1:res.shape[0],:]
res = res[:,20:156]
#f, (ax1,ax2,ax3) = pl.subplots(3,1,sharex=True)
win = AP.dsp.gen_window(136,window='blackman-harris')
#dirty = n.abs(dirty)
#res = n.abs(res)
freqs = n.linspace(100,200,203)
delaybins = n.linspace(-101,101,203)
if opts.dt == True:
    print 'TRUE'
    fig, axes = plt.subplots(nrows=2,ncols=2)
    axes = axes.flatten()
    DIRTY_ = n.fft.ifftshift(n.fft.ifft(dirty*win,axis=1),axes=1)
    RES_ = n.fft.ifftshift(n.fft.ifft(res*win,axis=1),axes=1)
    print axes.shape
else:
    fig, axes = plt.subplots(nrows=1,ncols=2)
dirty = n.abs(dirty)
res = n.abs(res)

im=axes[0].imshow(n.log10(dirty),aspect='auto',interpolation='none')
axes[0].set_title('Dirty')
im=axes[1].imshow(n.log10(res),aspect='auto',interpolation='none')
axes[1].set_title('Residual')
if opts.dt == True:
    im2=axes[2].imshow(n.log10(n.abs(DIRTY_)),aspect='auto',interpolation='none')
    axes[2].set_title('Dirty Delay Trans.')
    im2=axes[3].imshow(n.log10(n.abs(RES_)),aspect='auto',interpolation='none')
    axes[3].set_title('Residual Delay Trans.')
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im, cax=cbar_ax)
#fig.colorbar(im2,cax=cbar_
plt.savefig('/users/jkerriga/Plots/WaterfallResDirty.png')
#pl.show()
