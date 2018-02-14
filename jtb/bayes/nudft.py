import numpy as np
import matplotlib.gridspec as gridspec
import optparse, sys, os

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
o.add_option('--cmap',
    type = str,
    help = 'Colormap to use for plotting.',
    default = 'gnuplot')
o.add_option('--plot_phase',
    action = 'store_true',
    help = 'If passed, also plot and compare phase.')
o.add_option('--l_offset',
    type = float,
    help = 'Moves source in l-direction by l_offset*ls.max().  Must be between 0 and 1.')
o.add_option('--m_offset',
    type = float,
    help = 'Moves source in m-direction by m_offset*ms.max().  Must be between 0 and 1.')
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

if opts.l_offset:
    l_off = int(mid_l*opts.l_offset)
else:
    l_off = 0
if opts.m_offset:
    m_off = int(mid_m*opts.m_offset)
else:
    m_off = 0

I[mid_m + m_off, mid_l + l_off] = 1.0
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
Vs_func = lambda u,v: pixel_area*np.exp(-2*np.pi*1j*(u*ls[mid_l + l_off] + v*ms[mid_m + m_off]))
Vs_analytic = Vs_func(uvs_unique[:,0],uvs_unique[:,1])




## ----------------- Plotting ----------------- ##

if opts.plot_phase:
    fig = figure(figsize = (16,8))
    gs = gridspec.GridSpec(2,4)

    imax = subplot(gs[:,0])
    visax = subplot(gs[0,1])
    vispax = subplot(gs[1,1])
    anax = subplot(gs[0,2])
    anpax = subplot(gs[1,2])
    diffax = subplot(gs[0,3])
    diffpax = subplot(gs[1,3])
else:
    fig = figure(figsize = (16,4))
    gs = gridspec.GridSpec(1,4)

    imax = subplot(gs[0])
    visax = subplot(gs[1])
    anax = subplot(gs[2])
    diffax = subplot(gs[3])


aspect = 'equal'

ogim = imax.imshow(I,
                                interpolation = 'nearest',
                                origin = 'lower',
                                extent = extent_lm,
                                aspect = aspect,
                                cmap = opts.cmap)
imax.set_title('sky')
imax.set_xlabel('l', size=16)
imax.set_ylabel('m', size=16, rotation=0)


myim = visax.scatter(uvs_unique[:,0],
                                uvs_unique[:,1],
                                c = np.abs(Vs_unique),
                                cmap = opts.cmap)
visax.set_title('my fft')


anim = anax.scatter(uvs_unique[:,0],
                              uvs_unique[:,1],
                              c = np.abs(Vs_analytic),
                              cmap = opts.cmap)
anax.set_title('analytic solution')


diffim = diffax.scatter(uvs_unique[:,0],
                                uvs_unique[:,1],
                                c = np.abs(Vs_unique) - np.abs(Vs_analytic),
                                cmap = opts.cmap)
diffax.set_title('my fft - analytic solution')

if opts.plot_phase:
    mypim = vispax.scatter(uvs_unique[:,0],
                                       uvs_unique[:,1],
                                       c = np.imag(Vs_unique),
                                       cmap = opts.cmap)

    anpim = anpax.scatter(uvs_unique[:,0],
                                     uvs_unique[:,1],
                                     c = np.imag(Vs_analytic),
                                     cmap = opts.cmap)

    diffpim = diffpax.scatter(uvs_unique[:,0],
                                        uvs_unique[:,1],
                                        c = np.imag(Vs_unique) - np.imag(Vs_analytic),
                                        cmap = opts.cmap)

    imgs = [ogim, myim, mypim, anim, anpim, diffim, diffpim]

else:
    imgs = [ogim, myim, anim, diffim]


for a in fig.axes[1:]:
    a.set_xlabel('u', size = 16)
    a.set_ylabel('v', size = 16, rotation = 0)
    a.set_aspect('equal')

for i,axe in enumerate(fig.axes):
    divider = make_axes_locatable(axe)
    cax = divider.append_axes("right", size = "5%", pad = 0.05)
    colorbar(imgs[i], cax = cax)

gs.tight_layout(fig)

show()

sys.exit()
