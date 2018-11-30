import numpy as np
import healpy as hp
import time, optparse, sys, os, re

from astropy.io import fits
from scipy import interpolate


## ----------------- Option Parser ----------------- ##
o = optparse.OptionParser()

o.add_option('--npix_side',
    type = int,
    default = 31,
    help = 'Number of pixels on a side for the lm and uv grids.')

o.add_option('--fov',
    type = float,
    default = 10,
    help = 'Field of view in degrees.')

opts,args = o.parse_args(sys.argv[1:])

# Read in cst files
filenames = args
print 'Reading in %d files...' %len(filenames)
freqs = np.array([float(re.findall(r'\d+', f.split('_')[-1])[0]) for f in filenames])
order = np.argsort(freqs)
freqs = freqs[order]
nfreqs = len(freqs)
filenames = np.array(filenames)[order]
nfiles = len(filenames)

# Get sky coordinates from cst output
data = np.loadtxt(filenames[0], skiprows=2, usecols=[0, 1])
thetas, phis = np.deg2rad(data[:, 0]), np.deg2rad(data[:, 1])
npix = thetas.size
coords = np.stack((thetas, phis)).T
# Convert to l,m coordinates
# ls_cst = np.cos(phis)*np.sin(thetas)
# ms_cst = np.sin(phis)*np.sin(thetas)
ls_cst = thetas*np.cos(phis)
ms_cst = thetas*np.sin(phis)
# sort_inds = np.lexsort((coords[:, 1], coords[:, 0]))

# try interpolating to thetas, phis and then using indices
# for l, m from converted thetas, phis

# Construct sky coordinates to interpolate to
if opts.fov != 10.:
    npix_side = int(opts.fov*31/10)
    if npix_side%2 == 0:
        npix_side += 1
else:
    npix_side = opts.npix_side
FOV = np.round(np.deg2rad(opts.fov), decimals=14)
ls = np.linspace(-FOV/2, FOV/2, npix_side)
ms = np.copy(ls)
nlm = ls.size*ms.size
lm_pixel_half = np.diff(ls)[0]/2.
L, M = np.meshgrid(ls, ms)
ls_vec, ms_vec = L.flatten(), M.flatten()
# l_inds = np.logical_and(ls_cst >= ls.min(), ls_cst <= ls.max())
# m_inds = np.logical_and(ms_cst >= ms.min(), ms_cst <= ms.max())
# lm_inds = l_inds*m_inds*(thetas <= np.pi/2)
interp_thetas = np.sqrt(ls_vec**2 + ms_vec**2)
interp_phis = np.arctan2(ms_vec,ls_vec)
interp_phis[interp_phis < 0.0] += 2*np.pi

# Get thetas, phis that correspond to ls, ms
# interp_coords = np.stack((ls_cst, ms_cst)).T

# Per frequency interpolation
cst_beam = np.zeros((npix, nfiles))
interp_beam = np.zeros((nlm, nfreqs))

# # griddata stuff just in case
# theta_vals = np.tile(np.arange(0, 181), 360)
# phi_vals = np.arange(0, 360).reshape((-1, 1))
# phi_vals = np.tile(phi_vals, 181).flatten()

for fi, f in enumerate(filenames):
   cst_beam[:, fi] = np.loadtxt(f, skiprows=2, usecols=2)
   # interp_func = interpolate.NearestNDInterpolator(interp_coords, cst_beam[:, fi])
   #
   # interp2d is not bad, but there should be some asymmetry in the beam that's no there in imshow
   # interp_func = interpolate.interp2d(ls_cst,
   #                                    ms_cst,
   #                                    cst_beam[:, fi],
   #                                    kind = 'linear')
   #
   # This is the one that moves the beam maximum around when first going to healpix
   # interp_func = interpolate.RectBivariateSpline(np.unique(coords[sort_inds, 0]),
   #                                               np.unique(coords[sort_inds, 1]),
   #                                               cst_beam[sort_inds, fi].reshape([181,360], order='F'))
   # for i in range(nlm):
   #    interp_beam[i, fi] = interp_func(ls_vec[i], ms_vec[i])
   #
   interp_beam[:, fi] = interpolate.griddata((thetas, phis),
                                             cst_beam[:, fi],
                                             (interp_thetas, interp_phis),
                                             method = 'linear')

# ---------------------------------- Plotting ---------------------------------- #
import matplotlib.gridspec as gridspec

from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib.pyplot import *

fig = figure(figsize=(9, 4))
gs = GridSpec(1, 2)
imgs = []
extent_lm = [ls.min() - lm_pixel_half, ls.max() + lm_pixel_half,
             ms.min() - lm_pixel_half, ms.max() + lm_pixel_half]
extent_lm = np.rad2deg(extent_lm)
freq_ind = 0

interp_ax = fig.add_subplot(gs[0])
im = interp_ax.imshow(interp_beam[:, freq_ind].reshape([npix_side]*2),
                      extent = extent_lm,
                      origin = 'lower')
interp_ax.set_xlabel('l [deg]', size=16)
interp_ax.set_ylabel('m [deg]', size=16)
interp_ax.set_title('Interpolated Beam')
imgs.append(im)

# Cuts through the beams plot
l0_inds = np.where(ls_vec == 0.0)[0]
m0_inds = np.where(ms_vec == 0.0)[0]
phi0_inds = np.where(np.logical_or(phis == np.pi, phis == 0))[0]
phi90_inds = np.where(np.logical_or(phis == np.pi/2, phis == 3*np.pi/2))[0]

cut_ax = fig.add_subplot(gs[1])
ymin = np.min((interp_beam[l0_inds, freq_ind].min(), interp_beam[m0_inds, 0].min()))
ymax = np.max((interp_beam[l0_inds, freq_ind].max(), interp_beam[m0_inds, 0].max()))
cut_ax.plot(ms_vec[l0_inds], interp_beam[l0_inds, freq_ind], 'r-', label='l = 0')
cut_ax.plot(ls_vec[m0_inds], interp_beam[m0_inds, freq_ind], 'b-', label='m = 0')
cut_ax.plot(ls_cst[phi0_inds], cst_beam[phi0_inds, freq_ind], 'bo', label='m_cst = 0')
cut_ax.plot(ms_cst[phi90_inds], cst_beam[phi90_inds, freq_ind], 'ro', label='l_cst = 0')
cut_ax.set_xlim([ls.min() - lm_pixel_half, ls.max() + lm_pixel_half])
cut_ax.set_ylim([ymin - 10, ymax + 10])
cut_ax.legend(loc='lower center', frameon=False)

for i,ax in enumerate(fig.axes[:-1]):
    ax_divider = make_axes_locatable(ax)
    cbar_ax = ax_divider.append_axes("right", size="5%", pad="2%")
    cb = fig.colorbar(imgs[i], cax=cbar_ax)

gs.tight_layout(fig)

show()
