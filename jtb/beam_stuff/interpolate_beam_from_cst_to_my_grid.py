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

o.add_option('--write',
    action = 'store_true',
    help = 'If passed, write interpolated beam to a numpy readable file.')

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
# Convert cst coordinates to l,m coordinates for plotting
ls_cst = thetas*np.cos(phis)
ms_cst = thetas*np.sin(phis)

# Construct sky coordinates to interpolate to
if opts.fov != 10. and opts.npix_side == 31:
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
interp_thetas = np.sqrt(ls_vec**2 + ms_vec**2)
interp_phis = np.arctan2(ms_vec,ls_vec)
interp_phis[interp_phis < 0.0] += 2*np.pi

# Per frequency interpolation
cst_beam = np.zeros((npix, nfiles))
interp_beam = np.zeros((nlm, nfiles))

for fi, f in enumerate(filenames):
   cst_beam[:, fi] = np.loadtxt(f, skiprows=2, usecols=2)
   interp_beam[:, fi] = interpolate.griddata((thetas, phis),
                                             cst_beam[:, fi],
                                             (interp_thetas, interp_phis),
                                             method = 'linear')

if opts.write:
    save_dir = '/users/jburba/data/jburba/beam_stuff/interp_beams/'
    filename = 'interp_beam_%dnpix-side_%ddfov_%.1f-%.1fMHz' %(npix_side,
                                                               opts.fov,
                                                               freqs.min(),
                                                               freqs.max())
    print 'Writing to %s' %(save_dir + filename + '.npy')
    np.save(save_dir + filename, interp_beam)

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
cut_ax.plot(np.rad2deg(ms_vec[l0_inds]), interp_beam[l0_inds, freq_ind], 'r-', label='l = 0')
cut_ax.plot(np.rad2deg(ms_cst[phi90_inds]), cst_beam[phi90_inds, freq_ind], 'ro', label='l_cst = 0')
cut_ax.plot(np.rad2deg(ls_vec[m0_inds]), interp_beam[m0_inds, freq_ind], 'b-', label='m = 0')
cut_ax.plot(np.rad2deg(ls_cst[phi0_inds]), cst_beam[phi0_inds, freq_ind], 'bo', label='m_cst = 0')
cut_ax.set_xlim(np.rad2deg([ls.min() - lm_pixel_half, ls.max() + lm_pixel_half]))
cut_ax.set_ylim([ymin - 10, ymax + 10])
cut_ax.set_xlabel('Distance from zenith [deg]')
cut_ax.legend(loc='lower center', frameon=False, ncol=2)

for i,ax in enumerate(fig.axes[:-1]):
    ax_divider = make_axes_locatable(ax)
    cbar_ax = ax_divider.append_axes("right", size="5%", pad="2%")
    cb = fig.colorbar(imgs[i], cax=cbar_ax)

gs.tight_layout(fig)

show()
