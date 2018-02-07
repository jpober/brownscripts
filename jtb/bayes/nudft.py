import numpy as np
import matplotlib.gridspec as gridspec
import optparse, sys, os

from scipy.integrate import dblquad
from matplotlib.pyplot import *
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Option parser
o = optparse.OptionParser()
o.add_option('--uvdata',
    type = str,
    help = 'Filename for input uvw binary .npy file.')
o.add_option('--pixel_count',
    type = int,
    help = 'Number of pixels along one axis in sky image.  Total number of pixels is pixel_count**2.',
    default = 50)
opts,args = o.parse_args(sys.argv[1:])

## ----------------- Construct sky ----------------- ##
ls = np.cos(np.linspace(-np.pi, np.pi, opts.pixel_count))
ms = np.sin(np.linspace(-np.pi, np.pi, opts.pixel_count))
N_im = ls.size*ms.size
extent_lm = [ls.min(), ls.max(), ms.min(), ms.max()]

L, M = np.meshgrid(ls, ms)
I = 0*L
if len(ls) % 2 == 0:
    mid_l = len(ls)/2
else:
    mid_l = int(len(ls)/2. - 0.5)
if len(ms) % 2 == 0:
    mid_m = len(ms)/2
else:
    mid_m = int(len(ms)/2. - 0.5)
I[mid_l, mid_m] = 1.0
N_freq = 1
pixel_area = (ls[1]-ls[0])*(ms[1]-ms[0])

# Construct position vectors for N_im pixels
I_vec = np.reshape(I, (N_im*N_freq,), order='F')
ls_vec, ms_vec = np.zeros(0), np.zeros(0)
for l in ls:
    for m in ms:
        ls_vec = np.append(ls_vec, l)
        ms_vec = np.append(ms_vec, m)



## ----------------- Construct uv ----------------- ##

# Read in u,v sampling points
uvws = np.load(opts.uvdata)
us = np.concatenate((-uvws[0,0,:], uvws[0,0,:]))
vs = np.concatenate((-uvws[0,1,:], uvws[0,1,:]))
N_vis = us.size*vs.size
extent_uv = [us.min(), us.max(), vs.min(), vs.max()]


## ----------------- Construct nudft ----------------- ##

# Construct DFT matrix using outer product
DFT = np.exp(-1j*2*np.pi*(np.outer(us,ls_vec) + np.outer(vs, ms_vec)))


Vs = np.dot(DFT, I_vec)*pixel_area

################# LEFT OFF HERE #################

# Need to average over redundant baselines
Vs_unique = np.zeros(0)
uvs_unique, loc_inds, inverse_inds = np.unique(np.stack((us, vs), axis=1),
                                                                      axis = 0,
                                                                      return_index = True,
                                                                      return_inverse = True)
for u,v in uvs_unique:
    avg_inds = np.where(np.logical_and(us == u, vs == v))[0]
    Vs_unique = np.append(Vs_unique, np.mean(Vs[avg_inds]))

################# LEFT OFF HERE #################



## ----------------- Analytic solution comparison ----------------- ##

# Point source, Flat beam
Vs_func = lambda u,v: pixel_area*np.exp(2*np.pi*1j*(u*ls[mid_l] + v*ms[mid_m]))
Vs_analytic = Vs_func(uvs_unique[:,0],uvs_unique[:,1])




## ----------------- Plotting ----------------- ##

fig = figure(figsize=(16,4))
gs = gridspec.GridSpec(1,4)

imgs = []
aspect = 'auto'
cmap = 'gnuplot2'

ax = subplot(gs[0])
ogim = ax.imshow(I,
                            interpolation = 'nearest',
                            origin = 'center',
                            extent = extent_lm,
                            aspect = aspect,
                            cmap = cmap)
ax.set_title('sky')
ax.set_xlabel('l', size=16)
ax.set_ylabel('m', size=16, rotation=0)

visax = subplot(gs[1])
myim = visax.scatter(uvs_unique[:,0],
                                uvs_unique[:,1],
                                c = np.abs(Vs_unique),
                                cmap = cmap)
visax.set_title('my fft')

anax = subplot(gs[2])
anim = anax.scatter(uvs_unique[:,0],
                              uvs_unique[:,1],
                              c = np.abs(Vs_analytic),
                              cmap = cmap)
anax.set_title('analytic solution')

diffax = subplot(gs[3])
diffim = diffax.scatter(uvs_unique[:,0],
                                uvs_unique[:,1],
                                c = np.abs(Vs_unique) - np.abs(Vs_analytic),
                                cmap = cmap)
diffax.set_title('my fft - analytic solution')

for a in fig.axes[1:]:
    a.set_xlabel('u', size=16)
    a.set_ylabel('v', size=16, rotation=0)
    a.set_aspect('equal')

imgs = [ogim, myim, anim, diffim]

for i,axe in enumerate(fig.axes):
    divider = make_axes_locatable(axe)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    colorbar(imgs[i], cax=cax)


gs.tight_layout(fig)

show()
