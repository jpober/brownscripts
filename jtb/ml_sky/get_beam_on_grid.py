import numpy as np
import healpy as hp
import matplotlib.gridspec as gridspec
import optparse, sys

from astropy.io import fits
from scipy.optimize import curve_fit
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib.pyplot import *


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

o.add_option('--fit_beam',
    action = 'store_true',
    help = 'If passed, fit beam with a 2d Gaussian and plot beam - fit.')

o.add_option('--log_scale',
    action = 'store_true',
    help = 'If passed, plot in dB.')

opts,args = o.parse_args(sys.argv[1:])

def twoD_Gaussian((x, y), amp, x0, y0, sigma_x, sigma_y):
    # See https://stackoverflow.com/questions/21566379/
    # fitting-a-2d-gaussian-function-using-scipy-optimize-curve-fit-valueerror-and-m
    # for example of 2d Guassian fitting using scipy.optimize.curve_fit
    return (amp*np.exp(-(x-x0)**2/(2*sigma_x**2) - (y-y0)**2/(2*sigma_y**2))).ravel()

# Construct l,m grid
NPIX_SIDE = opts.npix_side
FOV = np.deg2rad(opts.fov)
ls = np.linspace(-FOV/2, FOV/2, NPIX_SIDE)
ms = np.copy(ls)
extent_lm = [ls.min(), ls.max(), ms.min(), ms.max()]
extent_lm = np.rad2deg(extent_lm)
ls_vec, ms_vec = np.zeros(0), np.zeros(0)
for m in ms:
    for l in ls:
        ls_vec = np.append(ls_vec, l)
        ms_vec = np.append(ms_vec, m)

thetas = np.sqrt(ls_vec**2 + ms_vec**2)
thetas *= FOV/(2.*thetas.max())
phis = np.arctan2(np.copy(ms_vec), np.copy(ls_vec))
phis[phis < 0.0] += 2*np.pi

# Get beam on my (l, m) grid
fits_file = ('/Users/jburba/hera_things/hera-team/HERA-Beams/'
                +
                'NicolasFagnoniBeams/healpix_beam.fits')
beam_E = fits.getdata(fits_file, extname='BEAM_E')
freqs = fits.getdata(fits_file, extname='FREQS')
freq_ind = np.where(freqs == 150.)[0][0]

beam_grid = hp.get_interp_val(beam_E[:, freq_ind], thetas, phis)

if opts.fit_beam:
    L, M = np.meshgrid(ls, ms)
    fit_beam_grid = np.zeros_like(beam_grid)
    initial_guess = (beam_grid.max(), 0., 0., 1., 1.)
    popt, pcov = curve_fit(twoD_Gaussian, (L, M), beam_grid, p0=initial_guess)
    fit_beam_grid = twoD_Gaussian((L, M), *popt)

# Plotting
if opts.fit_beam:
    gs = gridspec.GridSpec(1, 3)
    fig = figure(figsize=(16, 4))
else:
    gs = gridspec.GridSpec(1, 1)
    fig = figure(figsize=(8, 8))

# thax = fig.add_subplot(gs[0])
# thim = imshow(np.rad2deg(thetas).reshape((NPIX_SIDE, NPIX_SIDE)),
#                         origin='lower',
#                         extent=extent_lm)
#
# phax = fig.add_subplot(gs[1])
# phim = imshow(np.rad2deg(phis).reshape((NPIX_SIDE, NPIX_SIDE)),
#                         origin='lower',
#                         extent=extent_lm)

beamax = fig.add_subplot(gs[0])
if opts.log_scale:
    beamim = imshow(10*np.log10(beam_grid).reshape((NPIX_SIDE, NPIX_SIDE)),
                                origin='lower',
                                extent=extent_lm)
    beamax.set_title('Beam')
else:
    beamim = imshow(beam_grid.reshape((NPIX_SIDE, NPIX_SIDE)),
                                origin='lower',
                                extent=extent_lm)
    beamax.set_title('Beam, %.1fMHz' %freqs[freq_ind])

if opts.fit_beam:
    fitax = fig.add_subplot(gs[1])
    if opts.log_scale:
        fitim = imshow(10*np.log10(fit_beam_grid).reshape([NPIX_SIDE]*2),
                               origin='lower',
                               extent=extent_lm)
        fitax.set_title('Fitted Beam')
    else:
        fitim = imshow(fit_beam_grid.reshape([NPIX_SIDE]*2),
                               origin='lower',
                               extent=extent_lm)
        fitax.set_title('Fitted Beam')

    diff_data = beam_grid - fit_beam_grid
    diffax = fig.add_subplot(gs[2])
    if opts.log_scale:
        diffim = imshow(10*np.log10(diff_data).reshape([NPIX_SIDE]*2),
                                 origin='lower',
                                 extent=extent_lm)
        diffax.set_title('Beam - Fitted Beam')
    else:
        diffim = imshow(diff_data.reshape([NPIX_SIDE]*2),
                                 origin='lower',
                                 extent=extent_lm)
        diffax.set_title('Beam - Fitted Beam')

# imgs = [thim, phim, beamim]
if opts.fit_beam:
    imgs = [beamim, fitim, diffim]
else:
    imgs = [beamim]

for i, ax in enumerate(fig.axes):
    ax_divider = make_axes_locatable(ax)
    cax = ax_divider.append_axes("right", size="5%", pad="2%")
    cb = fig.colorbar(imgs[i], cax=cax)
    if opts.fit_beam:
        if not imgs[i] == diffim:
            cb.set_label('dBi', size=16)
    ax.set_xlabel('l')
    ax.set_ylabel('m')

gs.tight_layout(fig)
show()
