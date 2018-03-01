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
o.add_option('--ps_catalog',
    type = str,
    default = None,
    help = 'Filename for binary numpy point source catalog.')
opts,args = o.parse_args(sys.argv[1:])

## ----------------- Construct sky ----------------- ##
pixel_num = 50
ls = np.cos(np.linspace(-np.pi, np.pi, pixel_num))
ms = np.sin(np.linspace(-np.pi, np.pi, pixel_num))
pixel_side_length = np.diff(ls)[0]
N_im = ls.size*ms.size
extent_lm = [ls.min(), ls.max(), ms.min(), ms.max()]

L, M = np.meshgrid(ls, ms)

#Make source catalog
nsources = 50
grid_pos = np.zeros((nsources,2), dtype = int)
true_pos = np.zeros((nsources,2))
for i in range(grid_pos.shape[0]):
    grid_pos[i, 0] = np.random.randint(0, ls.shape[0])
    grid_pos[i, 1] = np.random.randint(0, ms.shape[0])
    true_pos[i, 0] = ls[grid_pos[i, 0]] + np.random.uniform(low=-pixel_side_length, high=pixel_side_length)
    true_pos[i, 1] = ms[grid_pos[i, 1]] + np.random.uniform(low=-pixel_side_length, high=pixel_side_length)

I = 0*L
I[grid_pos[:,0], grid_pos[:,1]] = 1.0
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
DFT = np.exp(-1j*2*np.pi*(np.outer(us, ls_vec) + np.outer(vs, ms_vec)))

# Compute visibilities from DFT
Vs = np.dot(DFT, I_vec)

# Average over redundant baselines
Vs_unique = np.zeros(0)
uvs_unique = np.unique(np.stack((us, vs), axis=1), axis = 0)

for u,v in uvs_unique:
    avg_inds = np.where(np.logical_and(us == u, vs == v))[0]
    Vs_unique = np.append(Vs_unique, np.mean(Vs[avg_inds]))


## ----------------- Compute Power Spectrum ----------------- ##
# P(k_u,k_v) = sum_r(Visibilities(u,v)**2) s.t. sqrt(u**2+v**2) = r
# Averaged over circular annuli

rs = np.round(np.sqrt(uvs_unique[:,0]**2 + uvs_unique[:,1]**2), decimals=14)
rs_unique = np.unique(rs)
pspec = np.zeros_like(rs_unique)*1.0j

for i,r in enumerate(rs_unique):
    inds = np.where(rs == r)[0]
    pspec[i] = np.mean(np.abs(Vs_unique[inds])**2)


## ----------------- Analytic solution comparison ----------------- ##

# Point source, Flat beam
Vs_func = lambda u,l,v,m: np.exp(-2*np.pi*1j*(u*l + v*m))
Vs_analytic = np.zeros(uvs_unique.shape[0], dtype=complex)
for i,(u,v) in enumerate(uvs_unique):
    Vs_analytic[i] = np.sum(Vs_func(u, true_pos[:,0], v, true_pos[:,1]))
pspec_analytic = np.zeros_like(rs_unique)*1.0j
for i,r in enumerate(rs_unique):
    inds = np.where(rs == r)[0]
    pspec_analytic[i] = np.mean(np.abs(Vs_analytic[inds])**2)


## ----------------- Plotting ----------------- ##
fig = figure(figsize = (15,5))
gs = gridspec.GridSpec(1,3)

gridax = subplot(gs[0,0])
# gridax.scatter(L, M, marker='.', alpha=0.25)
gridax.hlines(ms, ls.min(), ls.max(), alpha=0.25)
gridax.vlines(ls, ms.min(), ms.max(), alpha=0.25)
gridax.scatter(true_pos[:,0], true_pos[:,1], marker='o', c='r', s=4, edgecolor='r')
gridax.set_xlabel('l', size=16)
gridax.set_ylabel('m', size=16)
gridax.set_xlim([ls.min(), ls.max()])
gridax.set_ylim([ms.min(), ms.max()])

pspecax = subplot(gs[0,1])
pspecax.plot(rs_unique, pspec, 'o-', label = 'Numerical')
pspecax.plot(rs_unique, pspec_analytic, 'o-', label = 'Analytic')
pspecax.set_xlabel(r'"k$_\perp$"', size = 16)
pspecax.set_ylabel('Power', size = 16)
pspecax.legend(loc='upper left', ncol=2)

diffax = subplot(gs[0,2])
diffax.plot(rs_unique, pspec - pspec_analytic, 'o-')
diffax.set_xlabel(r'"k$_\perp$"', size = 16)
diffax.set_ylabel('Power [Numerical - Analytic]', size = 16)

tight_layout()

show()
