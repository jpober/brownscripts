import numpy as np
import healpy as hp
import matplotlib.gridspec as gridspec
import time, optparse, sys, os

from astropy.io import fits
from scipy.optimize import curve_fit


## ----------------- Option Parser ----------------- ##
o = optparse.OptionParser()

o.add_option('--freq',
    type=str,
    default='150',
    help='Frequency(ies) of observation in MHz.')

o.add_option('--freq_res',
    type = float,
    default = 1.0,
    help = 'Channel width in MHz if --freq passed with \'-\' (1MHz for Nicolas beams).')

o.add_option('--grid_pos',
    action='store_true',
    default = 'True',
    help = 'If passed, compare with grid centers. '
               +
               'Otherwise compare with true source positions.')

o.add_option('--zenith_source',
    action = 'store_true',
    help = 'If passed, only place one source at zenith.')

o.add_option('--horizon_source',
    action = 'store_true',
    help = 'If passed, only place one source at the horizon.')

o.add_option('--uniform_sky',
    action = 'store_true',
    help = 'If passed, make sky uniform with amplitude 1.')

o.add_option('--noise_sky',
    action = 'store_true',
    help = 'If passed, draw sky from normal distribution with rms given by --rms.')

o.add_option('--l_offset',
    type = float,
    help = 'Moves source in l-direction by l_offset*ls.max().  Must be between 0 and 1.')

o.add_option('--m_offset',
    type = float,
    help = 'Moves source in m-direction by m_offset*ms.max().  Must be between 0 and 1.')

o.add_option('--nsources',
    type = int,
    default = 1,
    help= 'Number of sources to add to image')

o.add_option('--rms',
    type = float,
    default = 1.e-5,
    help = 'RMS for noise injection.')

o.add_option('--spec_index',
    type = float,
    default = '0.0',
    help = 'Spectral index for amplitude of point sources as a function of frequency,'
           + 'i.e. n = spectral index and amplitude A(v)=v^n.')

o.add_option('--log_scale',
    action = 'store_true',
    help = 'If passed, plot in log scale.')

o.add_option('--force_lim',
    action = 'store_true',
    help = 'If passed, force colorbars to same limits on all plots.')

o.add_option('--beam',
    type = 'str',
    help = 'Healpix compatible fits file containing the primary beam.')

o.add_option('--fit_beam',
    action = 'store_true',
    help = 'If passed, fit a two dimensional Gaussian to the beam for use in the maxL solution.')

o.add_option('--fractional_fit',
    action = 'store_true',
    help = 'If passed, normalize the beam at each frequency to have a peak at 1.')

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
    help = 'If passed, write fitted RMS and associated errors to a .npy file.')

o.add_option('--uvdata',
    type = str,
    help = 'Filename for input uvw binary .npy file.')

opts,args = o.parse_args(sys.argv[1:])


# Visibility function
def Vs_func(u, l, v, m):
    return np.exp(-2*np.pi*1j*(u*l+ v*m))

# Gaussian fitting function
def Gaussian(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def twoD_Gaussian((x, y), x0, y0, sigma_x, sigma_y):
    # See https://stackoverflow.com/questions/21566379/
    # fitting-a-2d-gaussian-function-using-scipy-optimize-curve-fit-valueerror-and-m
    # for example of 2d Guassian fitting using scipy.optimize.curve_fit
    return (np.exp(-(x-x0)**2/(2*sigma_x**2) - (y-y0)**2/(2*sigma_y**2))).ravel()

if opts.fractional_fit:
    # Force fit_beam = False to use the raw beam
    opts.fit_beam = False

# Constants
C = 3.e8 # meters per second
npix_side = opts.npix_side # 125 for 40deg FOV, 31 for 10deg FOV, 63 for 20deg FOV
npix = npix_side**2

# Get frequency(ies) or frequency range
if ',' in opts.freq:
    freqs = np.sort(map(float, opts.freq.split(',')))
elif '-' in opts.freq:
    freqs_arr = np.sort(map(float, opts.freq.split('-')))
    freqs = np.round(np.arange(freqs_arr[0], freqs_arr[1] + opts.freq_res, opts.freq_res),
                              decimals=3)
    freqs = freqs[np.where(freqs <= freqs_arr[1])]
elif not (',' in opts.freq and '-' in opts.freq):
    freqs = np.array([float(opts.freq)])

nfreqs = len(freqs)



## -------------------------------------- Construct sky -------------------------------------- ##
print 'Constructing sky...'

# Construct l,m grid
FOV = np.deg2rad(opts.fov)
ls = np.linspace(-FOV/2, FOV/2, npix_side)
ms = np.copy(ls)
nlm = ls.size*ms.size
lm_pixel_half = np.diff(ls)[0]/2.
ls_vec, ms_vec = np.zeros(0), np.zeros(0)
for m in ms:
    for l in ls:
        ls_vec = np.append(ls_vec, l)
        ms_vec = np.append(ms_vec, m)

#Make source catalog
nsources = opts.nsources
if len(ls) % 2 == 0:
    mid_l = len(ls)/2
else:
    mid_l = int(len(ls)/2. - 0.5)
if len(ms) % 2 == 0:
    mid_m = len(ms)/2
else:
    mid_m = int(len(ms)/2. - 0.5)

# grid_pos = np.zeros((nsources, 2), dtype = int)
# true_pos = np.zeros((nsources, 2))

Sky = np.zeros((nfreqs, npix_side, npix_side))

if opts.zenith_source:
    Sky[:, mid_m, mid_l] = 1.0

elif opts.horizon_source:
    Sky[:, 0, mid_l] = 1.0

elif opts.l_offset or opts.m_offset:
    if opts.l_offset:
        l_off = int(mid_l*opts.l_offset)
    else:
        l_off = 0

    if opts.m_offset:
        m_off = int(mid_m*opts.m_offset)
    else:
        m_off = 0

    Sky[:, mid_m + m_off, mid_l + l_off] = 1.0

elif opts.uniform_sky or opts.noise_sky:
    nsources = nlm
    if opts.noise_sky:
        Sky = np.ones_like(Sky)*np.random.normal(0, opts.rms, Sky.shape)
    else:
        Sky = np.ones_like(Sky)

else:
    for i in range(nsources):
        Sky[:, np.random.randint(0, ms.shape[0]), np.random.randint(0, ls.shape[0])] += 1.0

# if nsources == 1:
#     print 'Nsources: 1'
#     pos_vec = np.where(np.logical_and(ls_vec == ls[grid_pos[:, 1]],
#                                                         ms_vec == ms[grid_pos[:, 0]]))[0][0]
# else:
#     print 'Nsources: %d' %nsources
#     pos_vec = np.zeros(0, dtype=int)
#     for (l, m) in np.stack([ls[grid_pos[:, 1]], ms[grid_pos[:, 0]]], axis=1):
#         pos_vec = np.append(pos_vec, np.where(np.logical_and(ls_vec == l,
#                                                                                           ms_vec == m))[0][0])

# # Make sky matrix
# Sky = np.zeros((nfreqs, npix_side, npix_side))
# Sky_counts = np.zeros((npix_side, npix_side))
# for i, freq in enumerate(freqs):
#     if not opts.spec_index == 0.0:
#         Sky[i, grid_pos[:, 0], grid_pos[:, 1]] = 1./(1 + (freq - freqs.min())**opts.spec_index)
#     else:
#         Sky[i, grid_pos[:, 0], grid_pos[:, 1]] = 1.
#
# unique_grid_pos, unique_grid_pos_counts = np.unique(grid_pos,
#                                                                                  axis=0,
#                                                                                  return_counts=True)
# Sky_counts[unique_grid_pos[:, 0], unique_grid_pos[:, 1]] = unique_grid_pos_counts
#
# for i in range(freqs.size):
#     Sky[i] *= Sky_counts
#
# if opts.noise_sky:
#     Sky *= np.random.normal(0, opts.rms, Sky.shape)

# Flatten sky for matrix computation
Sky_vec = Sky.reshape((nfreqs, npix_side**2))

# Beam stuff
# if opts.beam:
#     # Get beam on sky grid
#     thetas = np.sqrt(ls_vec**2 + ms_vec**2)
#     thetas *= FOV/(2.*thetas.max())
#     phis = np.arctan2(np.copy(ms_vec), np.copy(ls_vec))
#     phis[phis < 0.0] += 2*np.pi
#
#     # Get beam on my (l, m) grid
#     beam_E = fits.getdata(opts.beam, extname='BEAM_E')
#     beam_freqs = fits.getdata(opts.beam, extname='FREQS')
#     beam_grid = np.zeros((nfreqs, thetas.size))
#
#     if opts.fit_beam:
#         L, M = np.meshgrid(ls, ms)
#         fit_beam_grid = np.zeros_like(beam_grid)
#
#     # Get beam on grid
#     ref_beam = None
#     for i, freq in enumerate(freqs):
#         freq_ind = np.where(beam_freqs == freq)[0][0]
#         beam_grid[i] = hp.get_interp_val(beam_E[:, freq_ind], thetas, phis)
#         beam_grid[i] /= beam_grid[i].max()
#
#         if opts.fractional_fit:
#             if ref_beam is None:
#                 print 'Setting reference beam at %.3fMHz' %freq
#                 ref_beam = np.copy(beam_grid[i])
#             beam_grid[i] /= ref_beam
#
#         if opts.fit_beam:
#             initial_guess = (0., 0., 1., 1.)
#             # popt, pcov = curve_fit(lambda (l, m), l0, m0, sl, sm: twoD_Gaussian((l, m), l0, m0, sl, sm), (L, M), beam_grid[i], p0=initial_guess)
#             popt, pcov = curve_fit(twoD_Gaussian, (L, M), beam_grid[i], p0=initial_guess)
#             fit_beam_grid[i] = twoD_Gaussian((L, M), *popt)

## ----------------------------------- Construct Visibilities ----------------------------------- ##
print 'Constructing visibilities...'

if opts.uvdata:
    # Can we even use a non gridded case for this scenario?
    # This would produce a maxL sky with shape (nfreqs, nuv)
    # but what (l, m) do those nuv sky pixels correspond to?
    print 'Loading uv coverage...'
    uvws_meters = np.load(opts.uvdata)
    # uvs only contains unique (u,v) assuming perfect degeneracy
    # uvs has shape (nfreqs, ntimes, nblns, 3)
    us_vec = uvws_meters[:, 0]
    vs_vec = uvws_meters[:, 1]
    # us_vec and vs_vec have shape (nfreqs, nblns)
    nvis = us_vec.shape[1]
else:
    # Construct u,v grid
    us_grid = np.fft.fftshift(np.fft.fftfreq(ls.shape[0], d=np.mean(np.diff(ls))))
    vs_grid = np.fft.fftshift(np.fft.fftfreq(ms.shape[0], d=np.mean(np.diff(ms))))
    us_vec, vs_vec = np.zeros(0), np.zeros(0)
    for v in vs_grid:
        for u in us_grid:
            us_vec = np.append(us_vec, u)
            vs_vec = np.append(vs_vec, v)
    uv_pixel_half = np.mean(np.diff(us_grid))/2.
    nvis = us_vec.shape[0]

# Use analytical solution to get visibilities using true positions
Vs = np.zeros((nfreqs, nvis), dtype=complex)
Vs_fft = np.zeros((nfreqs, nvis), dtype=complex)
Vs_mat = np.zeros((nfreqs, nvis), dtype=complex)

DFT = np.exp(-1j*2*np.pi*(np.outer(us_vec, ls_vec)
                    +
                    np.outer(vs_vec, ms_vec)))

for i in range(nfreqs):
    Vs_mat[i] = np.dot(DFT, Sky_vec[i])
    # Construct visibilities faster as FT of beam*sky
    if opts.beam:
        Vs_fft[i] = np.fft.fftshift(np.fft.fft2((beam_grid[i]*Sky_vec[i]).reshape([npix_side]*2))).flatten()
    else:
        Vs_fft[i] = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift(Sky_vec[i].reshape([npix_side]*2)))).flatten()


for i in range(nfreqs):
    for j in range(us_vec.size):
        Vs[i, j] = np.sum(Sky_vec[i]*Vs_func(us_vec[j], ls_vec, vs_vec[j], ms_vec))

## ---------------------------------- Plotting ---------------------------------- ##
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib.pyplot import *

fontsize = 14
extent_lm = [ls.min() - lm_pixel_half, ls.max() + lm_pixel_half,
                    ms.min() - lm_pixel_half, ms.max() + lm_pixel_half]
extent_lm = np.rad2deg(extent_lm)

extent_uv = [us_grid.min() - uv_pixel_half, us_grid.max() + uv_pixel_half,
                    vs_grid.min() - uv_pixel_half, vs_grid.max() + uv_pixel_half]

fig = figure(figsize=(16, 8))
gs = gridspec.GridSpec(2, 4)
imgs = []

ax = fig.add_subplot(gs[0, 0])
imgs.append(ax.imshow(Sky[0], extent=extent_lm, origin='lower'))
ax.set_title('Sky')

ax1 = fig.add_subplot(gs[0, 1])
imgs.append(ax1.imshow(np.abs(Vs[0]).reshape([npix_side]*2), extent=extent_uv, origin='lower'))
ax1.set_title('For Loop Visibilities')

ax2 = fig.add_subplot(gs[0, 2])
imgs.append(ax2.imshow(np.abs(Vs_mat[0]).reshape([npix_side]*2), extent=extent_uv, origin='lower'))
ax2.set_title('Matrix DFT Visibilities')

ax3 = fig.add_subplot(gs[0, 3])
imgs.append(ax3.imshow(np.abs(Vs_fft[0]).reshape([npix_side]*2), extent=extent_uv, origin='lower'))
ax3.set_title('Numpy Visibilities')

ax4 = fig.add_subplot(gs[1, 1])
imgs.append(ax4.imshow(np.angle(Vs[0]).reshape([npix_side]*2), extent=extent_uv, origin='lower'))
ax4.set_title('Phase')

ax5 = fig.add_subplot(gs[1, 2])
imgs.append(ax5.imshow(np.angle(Vs_mat[0]).reshape([npix_side]*2), extent=extent_uv, origin='lower'))
ax5.set_title('Phase')

ax6 = fig.add_subplot(gs[1, 3])
imgs.append(ax6.imshow(np.angle(Vs_fft[0]).reshape([npix_side]*2), extent=extent_uv, origin='lower'))
ax6.set_title('Phase')

# ax4 = fig.add_subplot(gs[1, 1:3])
# imgs.append(ax4.imshow((np.abs(Vs) - np.abs(Vs_mat)).reshape([npix_side]*2), extent=extent_uv))
# ax4.set_title('For loop - Matrix DFT')
#
# ax5 = fig.add_subplot(gs[1, 2:])
# imgs.append(ax5.imshow((np.abs(Vs) - np.abs(Vs_fft)).reshape([npix_side]*2), extent=extent_uv))
# ax5.set_title('For loop - Numpy FFT')

for i, ax in enumerate(fig.axes):
    ax_divider = make_axes_locatable(ax)
    cax = ax_divider.append_axes("right", size="5%", pad="2%")
    cb = fig.colorbar(imgs[i], cax=cax)
    if i == (len(fig.axes) - 1):
        cb.set_label('Vis Units', size=fontsize)

gs.tight_layout(fig)
# gs.update(wspace=0.05)

show()
