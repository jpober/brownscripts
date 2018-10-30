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

o.add_option('--gaussian_source',
    action = 'store_true',
    help = 'If passed, use a narrow gaussian source centered at (l_off, m_off).')

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
    # default = ('/Users/jburba/hera_things/hera-team/HERA-Beams/'
    #                 +
    #                 'NicolasFagnoniBeams/healpix_beam.fits'),
    help = 'Healpix compatible fits file containing the primary beam.')

o.add_option('--fit_beam',
    action = 'store_true',
    help = 'If passed, fit a two dimensional Gaussian to the beam for use in the maxL solution.')

o.add_option('--gof_inds',
    action = 'store_true',
    help = 'If passed, use the fractional difference between the true and fitted beams to determine '
               'what pixels to use for construction of the maximum likelihood solution for the sky via --fit_tol.')

o.add_option('--fit_tol',
    type = 'float',
    default = 0.1,
    help = 'Cutoff value for --gof_inds.  If fit_tol=0.1 (default), then all pixels where the fractional difference between '
               'the true and fitted beam is greater than 10% will be ignored in the maximum likelihood sky construction.')

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

o.add_option('--inset_fov',
    type = float,
    default = 10,
    help = 'Field of view in degrees to use when constructing the maximum likelihood solution for the sky.')

o.add_option('--write',
    action = 'store_true',
    help = 'If passed, write fitted RMS and associated errors to a .npy file.')

o.add_option('--uvdata',
    type = str,
    help = 'Filename for input uvw binary .npy file.')

o.add_option('--poly_order',
    type = str,
    # default = '2',
    help = 'If n, fit beam variation as a function of frequency with an n-th order polynomial.')

o.add_option('--data',
    type = str,
    default = '',
    help = 'File path for previously generated rms data via this script.')

opts,args = o.parse_args(sys.argv[1:])
print o.values

if opts.inset_fov > opts.fov:
    opts.inset_fov = opts.fov

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
# FOV = np.round(np.deg2rad(opts.fov), decimals=14)
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

# Build sky from opts params
Sky = np.zeros((nfreqs, npix_side, npix_side))

if opts.zenith_source:
    Sky[:, mid_m, mid_l] = 1.0

elif opts.horizon_source:
    Sky[:, 0, mid_l] = 1.0

elif (opts.l_offset or opts.m_offset) and not opts.gaussian_source:
    if opts.l_offset:
        l_off = int(mid_l*opts.l_offset)
    else:
        l_off = 0

    if opts.m_offset:
        m_off = int(mid_m*opts.m_offset)
    else:
        m_off = 0

    print '\tAdding source at l=%.1f, m=%.1f...' %(np.rad2deg(ls[mid_l + l_off]), np.rad2deg(ms[mid_m + m_off]))
    Sky[:, mid_m + m_off, mid_l + l_off] = 1.0

elif opts.uniform_sky or opts.noise_sky:
    nsources = nlm
    if opts.noise_sky:
        Sky = np.ones_like(Sky)*np.random.normal(0, opts.rms, Sky.shape)
    else:
        Sky = np.ones_like(Sky)

elif opts.gaussian_source:
    L, M = np.meshgrid(ls, ms)
    stddev = 0.03
    if opts.l_offset:
        l_off = int(mid_l*opts.l_offset)
    else:
        l_off = 0

    if opts.m_offset:
        m_off = int(mid_m*opts.m_offset)
    else:
        m_off = 0
    print '\tAdding Gaussian source at l=%.1f, m=%.1f...' %(np.rad2deg(ls[mid_l + l_off]), np.rad2deg(ms[mid_m + m_off]))
    Sky_vec = twoD_Gaussian((L, M), ls[mid_l + l_off], ms[mid_m + m_off], stddev, stddev)
    Sky_vec += twoD_Gaussian((L, M), ls[int(mid_l*(1-0.5))], 0, stddev, stddev)
    Sky = np.tile(Sky_vec, nfreqs).reshape((nfreqs, npix_side**2))

else:
    for i in range(nsources):
        Sky[:, np.random.randint(0, ms.shape[0]), np.random.randint(0, ls.shape[0])] += 1.0

# Flatten sky for matrix computation
Sky_vec = Sky.reshape((nfreqs, npix_side**2))

# Beam stuff
if opts.beam:
    # Get beam on sky grid
    thetas = np.sqrt(ls_vec**2 + ms_vec**2)
    thetas *= FOV/(2.*thetas.max())
    phis = np.arctan2(np.copy(ms_vec), np.copy(ls_vec))
    phis[phis < 0.0] += 2*np.pi

    # Get beam on my (l, m) grid
    if 'NF' in opts.beam:
        hdulist = fits.open(opts.beam)
        beam_E = np.copy(hdulist[0].data).squeeze()[0].T
        beam_freqs = np.arange(100, 201, 1)
    else:
        beam_E = fits.getdata(opts.beam, extname='BEAM_E')
        beam_freqs = fits.getdata(opts.beam, extname='FREQS')
    beam_grid = np.zeros((nfreqs, thetas.size))

    if opts.fit_beam:
        L, M = np.meshgrid(ls, ms)
        fit_beam_grid = np.zeros_like(beam_grid)

    # Get beam on grid
    ref_beam = None
    for i, freq in enumerate(freqs):
        freq_ind = np.where(beam_freqs == freq)[0][0]
        beam_grid[i] = hp.get_interp_val(beam_E[:, freq_ind], thetas, phis)
        beam_grid[i] /= beam_grid[i].max()

        if opts.fractional_fit:
            if ref_beam is None:
                print 'Setting reference beam at %.3fMHz' %freq
                ref_beam = np.copy(beam_grid[i])
            beam_grid[i] /= ref_beam

        if opts.fit_beam:
            initial_guess = (0., 0., 1., 1.)
            popt, pcov = curve_fit(twoD_Gaussian, (L, M), beam_grid[i], p0=initial_guess)
            print 'Fitting beam...'
            fit_beam_grid[i] = twoD_Gaussian((L, M), *popt)

# inset_fov = np.round(np.deg2rad(opts.inset_fov), decimals=14)
inset_fov = np.deg2rad(opts.inset_fov)
inset_fov_l_inds = np.logical_and(ls_vec >= -inset_fov/2, ls_vec <= inset_fov/2)
inset_fov_m_inds = np.logical_and(ms_vec >= -inset_fov/2, ms_vec <= inset_fov/2)
inset_fov_inds = inset_fov_l_inds*inset_fov_m_inds
ls_inset_vec = ls_vec[inset_fov_inds]
ms_inset_vec = ms_vec[inset_fov_inds]
ls_inset_grid = np.unique(ls_inset_vec)
ms_inset_grid = np.copy(ls_inset_grid)

## ----------------------------------- Construct Visibilities ----------------------------------- ##
print 'Constructing visibilities...'

# Construct u,v grid using the inset fov
us_inset_grid = np.fft.fftshift(np.fft.fftfreq(ls_inset_grid.shape[0], d=np.mean(np.diff(ls_inset_grid))))
vs_inset_grid = np.fft.fftshift(np.fft.fftfreq(ms_inset_grid.shape[0], d=np.mean(np.diff(ms_inset_grid))))
U_inset, V_inset = np.meshgrid(us_inset_grid, vs_inset_grid)
us_inset_vec, vs_inset_vec = U_inset.flatten(), V_inset.flatten()
uv_pixel_half = np.mean(np.diff(us_inset_grid))/2.
nvis = us_inset_vec.shape[0]

# Use analytical solution to get visibilities using true positions
Vs = np.zeros((nfreqs, nvis), dtype=complex)
Vs_inset = np.zeros((nfreqs, nvis), dtype=complex)

DFT = np.exp(-1j*2*np.pi*(np.outer(us_inset_vec, ls_vec)
                    +
                    np.outer(vs_inset_vec, ms_vec)))
DFT_inset = np.exp(-1j*2*np.pi*(np.outer(us_inset_vec, ls_inset_vec)
                    +
                    np.outer(vs_inset_vec, ms_inset_vec)))


for i in range(nfreqs):
    # Construct visibilities faster as FT of beam*sky
    if opts.beam:
        print 'Using beam'
        Vs[i] = np.dot(DFT, beam_grid[i]*Sky_vec[i])
        Vs_inset[i] = np.dot(DFT_inset, beam_grid[i, inset_fov_inds]*Sky_vec[i, inset_fov_inds])
    else:
        Vs[i] = np.dot(DFT, Sky_vec[i])
        Vs_inset[i] = np.dot(DFT_inset, Sky_vec[i, inset_fov_inds])


## ---------------------------------- Plotting ---------------------------------- ##
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib.pyplot import *

print 'Plotting...'

plot_size = 3
nrows = 2
ncols = 4
nplots = nrows*ncols
fig = figure(figsize=(plot_size*ncols, plot_size*nrows))

fontsize = 14
extent_lm = [ls.min() - lm_pixel_half, ls.max() + lm_pixel_half,
                    ms.min() - lm_pixel_half, ms.max() + lm_pixel_half]
extent_lm_inset_fov = [ls_inset_grid.min() - lm_pixel_half, ls_inset_grid.max() + lm_pixel_half,
                    ms_inset_grid.min() - lm_pixel_half, ms_inset_grid.max() + lm_pixel_half]
extent_lm = np.rad2deg(extent_lm)
extent_lm_inset_fov = np.rad2deg(extent_lm_inset_fov)

extent_uv_inset = [us_inset_grid.min() - uv_pixel_half, us_inset_grid.max() + uv_pixel_half,
                             vs_inset_grid.min() - uv_pixel_half, vs_inset_grid.max() + uv_pixel_half]

gs = gridspec.GridSpec(nrows, ncols, width_ratios=[1, 1, 3, 3])
imgs = []

# SKY AXIS
sky_ax = fig.add_subplot(gs[0, :2])
if opts.beam:
    if opts.log_scale:
        im = sky_ax.imshow(np.log10((Sky_vec[0]*beam_grid[0])).reshape([ls.size]*2),
                                               origin='lower',
                                               extent=extent_lm)
    else:
        im = sky_ax.imshow((Sky_vec[0]*beam_grid[0]).reshape([ls.size]*2),
                                        origin = 'lower',
                                        extent = extent_lm)
    sky_ax.set_title('Sky*Beam')
else:
    if opts.log_scale:
        im = sky_ax.imshow(np.log10(Sky_vec[0]).reshape([ls.size]*2),
                                               origin='lower',
                                               extent=extent_lm)
    else:
        im = sky_ax.imshow((Sky_vec[0]).reshape([ls.size]*2),
                                        origin = 'lower',
                                        extent = extent_lm)
    sky_ax.set_title('Sky')

if inset_fov < FOV:
    vlines(np.rad2deg([ls_inset_grid.min(), ls_inset_grid.max()]),
                np.rad2deg(ms_inset_grid.min()),
                np.rad2deg(ms_inset_grid.max()))
    hlines(np.rad2deg([ms_inset_grid.min(), ms_inset_grid.max()]),
                np.rad2deg(ls_inset_grid.min()),
                np.rad2deg(ls_inset_grid.max()))

sky_ax.set_xlabel('l [deg]', size=fontsize)
sky_ax.set_ylabel('m [deg]', size=fontsize)
imgs.append(im)

# BEAM AXIS
if opts.beam:
    beam_ax = fig.add_subplot(gs[1, :2])

    im = beam_ax.imshow(np.log10((beam_grid[0])).reshape([ls.size]*2),
                                       origin='lower',
                                       extent=extent_lm)
    beam_ax.set_title('log(Beam)')


    if inset_fov < FOV:
        vlines(np.rad2deg([ls_inset_grid.min(), ls_inset_grid.max()]),
                    np.rad2deg(ms_inset_grid.min()),
                    np.rad2deg(ms_inset_grid.max()))
        hlines(np.rad2deg([ms_inset_grid.min(), ms_inset_grid.max()]),
                    np.rad2deg(ls_inset_grid.min()),
                    np.rad2deg(ls_inset_grid.max()))

    beam_ax.set_xlabel('l [deg]', size=fontsize)
    beam_ax.set_ylabel('m [deg]', size=fontsize)
    imgs.append(im)

# VIS AXES
vs_ax = fig.add_subplot(gs[0, 2])
if opts.log_scale:
    im = vs_ax.imshow(np.log10(np.abs(Vs)[0]).reshape([us_inset_grid.size]*2),
                                           origin='lower',
                                           extent=extent_uv_inset)
    vs_ax.set_title('log|Vs| (full fov)')
else:
    im = vs_ax.imshow((np.abs(Vs)[0]).reshape([us_inset_grid.size]*2),
                                    origin = 'lower',
                                    extent = extent_uv_inset)
    vs_ax.set_title('|Vs| (full fov)')
imgs.append(im)

vs_inset_ax = fig.add_subplot(gs[0, 3])
if opts.log_scale:
    im = vs_inset_ax.imshow(np.log10(np.abs(Vs_inset)[0]).reshape([us_inset_grid.size]*2),
                                           origin='lower',
                                           extent=extent_uv_inset)
    vs_inset_ax.set_title('log|Vs| (inset fov)')
else:
    im = vs_inset_ax.imshow((np.abs(Vs_inset)[0]).reshape([us_inset_grid.size]*2),
                                           origin = 'lower',
                                           extent = extent_uv_inset)
    vs_inset_ax.set_title('|Vs| (inset fov)')
imgs.append(im)

vs_phase_ax = fig.add_subplot(gs[1, 2])
im = vs_phase_ax.imshow((np.angle(Vs)[0]).reshape([us_inset_grid.size]*2),
                                         origin = 'lower',
                                         extent = extent_uv_inset)
vs_phase_ax.set_title('Phase')
imgs.append(im)

vs_phase_inset_ax = fig.add_subplot(gs[1, 3])
im = vs_phase_inset_ax.imshow((np.angle(Vs_inset)[0]).reshape([us_inset_grid.size]*2),
                                                 origin = 'lower',
                                                 extent = extent_uv_inset)
vs_phase_inset_ax.set_title('Phase')
imgs.append(im)

for i,ax in enumerate(fig.axes):
    ax_divider = make_axes_locatable(ax)
    cax = ax_divider.append_axes("right", size="5%", pad="2%")
    cb = fig.colorbar(imgs[i], cax=cax)
    if opts.beam:
        if ax == beam_ax:
            cb.set_label('log', size=fontsize)
    if opts.log_scale:
        cb.set_label('log', size=fontsize)


gs.tight_layout(fig)
show()
