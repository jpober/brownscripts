import numpy as n
import pyuvdata
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
dirtlist = sorted(glob(opts.dirty))
reslist = sorted(glob(opts.residual))
dNsamp = n.zeros((1,203))
rNsamp = n.zeros((1,203))
for i in dirtlist:
    print i
    mir = pyuvdata.miriad.Miriad()
    try:
        mir.read_miriad(i)
    except:
        pass
    bsl = mir.antnums_to_baseline(41,49)==mir.baseline_array
    dNsamp = n.vstack((dNsamp,mir.nsample_array[bsl][:,0,:,0]))
    del(bsl)
    del(mir)

for j in reslist:
    print j
    mir = pyuvdata.miriad.Miriad()
    try:
        mir.read_miriad(j)
    except:
        pass
    bsl = mir.antnums_to_baseline(41,49)==mir.baseline_array
    rNsamp = n.vstack((rNsamp,mir.nsample_array[bsl][:,0,:,0]))
    del(bsl)
    del(mir)
diff = rNsamp-dNsamp
print (diff[dNsamp==0]>0).sum()

for i in range(203):
    if (rNsamp[:,i]!=0).sum() == len(rNsamp[:,i]):
        print i
        break
    else:
        continue

for j in range(203):
    if (dNsamp[:,j]!=0).sum() == len(dNsamp[:,j]):
        print j
        break
    else:
        continue


print float((diff[dNsamp==0]>0).sum())/(dNsamp==0).sum()
diff[diff<0] = 0
f, (ax1,ax2,ax3) = pl.subplots(3,1,sharex=True)
im=ax1.imshow(dNsamp,aspect='auto',interpolation='none')
im=ax2.imshow(rNsamp,aspect='auto',interpolation='none')
im=ax3.imshow(diff,aspect='auto',interpolation='none')
f.subplots_adjust(right=0.8)
cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
f.colorbar(im, cax=cbar_ax)
pl.show()


pl.plot(n.sum(rNsamp,0),'.')
pl.plot(n.sum(dNsamp,0),'+')
pl.show()
