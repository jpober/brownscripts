#!/bin/env python

import matplotlib
matplotlib.use('Agg')
import numpy as np, pylab as p, healpy as hp
from scipy.io import readsav
import sys, optparse, time

op = optparse.OptionParser()

### If freq is specified, args = filenames and loop over those.
### If obs is specified, args = frequency channels and loop over these

op.add_option("-f", "--freq", help='Frequency channel (defaults to 0)')
op.add_option("-o", "--obs", help='Filename')
opts, args = op.parse_args(sys.argv[1:])

fil = None
vmin, vmax= None, None
for f in args:
     start = time.time()
     if opts.freq is None and opts.obs is None:
	print "Please specify a channel or a filename"
        sys.exit()
     if (not opts.freq is None) and (not opts.obs is None):
	print "Please use either -f or -o, not both."
	sys.exit()
     if not opts.freq is None:
	freq  = int(opts.freq)
	fname = f
     elif not opts.obs is None:
	fname = opts.obs
	freq = int(f)

     ofname = ".".join(fname.split('.')[:-1])+"_"+str(freq)+".png"
     print fname
     if (not opts.obs is None) and (fil is None):
	     fil = readsav(fname)
     if not opts.freq is None:
	     fil = readsav(fname)
     fig = p.figure()
     
     inds = fil['hpx_inds']
     nside = int(fil['nside'])
     vecs = hp.pixelfunc.pix2vec(nside, inds)
     mean_vec = (np.mean(vecs[0]), np.mean(vecs[1]), np.mean(vecs[2]))
     ### Generate a rotation matrix 
     dt, dp = hp.rotator.vec2dir(mean_vec,lonlat=True)

     mwp = hp.projector.MollweideProj(xsize=2000, rot=(dt,dp,0))
     i,j   = mwp.xy2ij(mwp.vec2xy(vecs[0],vecs[1], vecs[2]))
     imin, imax = min(i), max(i)
     jmin, jmax = min(j), max(j)
     
     map = np.zeros(hp.nside2npix(nside))
     map[inds] = fil['res_cube'][freq,:]
     if vmin is None or vmax is None: vmin, vmax = min(map[inds]), max(map[inds])
     fun = lambda x,y,z: hp.pixelfunc.vec2pix(nside,x,y,z,nest=False)
     arr = mwp.projmap(map, fun)
     arr = arr[imin:imax, jmin:jmax]
     
     p.imshow(arr, vmin= vmin, vmax = vmax)
     if not opts.obs is None:
	p.title("Chan="+str(freq))
     else:
	p.title(fname)
#     p.show()
     p.colorbar()
     p.savefig(ofname,bbox_inches='tight')
     print 'Elapsed: ', time.time() - start
