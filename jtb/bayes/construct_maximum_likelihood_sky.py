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
    help = 'If passed, compare with grid centers. '
               +
               'Otherwise compare with true source positions.')

o.add_option('--zenith_source',
    action = 'store_true',
    help = 'If passed, only place one source at zenith.')

o.add_option('--horizon_source',
    action = 'store_true',
    help = 'If passed, only place one source at the horizon.')

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

o.add_option('--beam',
    type = 'str',
    # default = ('/Users/jburba/hera_things/hera-team/HERA-Beams/'
    #                 +
    #                 'NicolasFagnoniBeams/healpix_beam.fits'),
    help = 'Healpix compatible fits file containing the primary beam.')

o.add_option('--write',
    action = 'store_true',
    help = 'If passed, write fitted RMS and associated errors to a .npy file.')

opts,args = o.parse_args(sys.argv[1:])


# Visibility function
def Vs_func(u, l, v, m):
    return np.exp(-2*np.pi*1j*(u*l+ v*m))

# Gaussian fitting function
def Gaussian(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


# Constants
C = 3.e8 # meters per second
NPIX_SIDE = 125 # 125 for 40deg FOV, 31 for 10deg FOV
NPIX = NPIX_SIDE**2

# Get frequency(ies) or frequency range
if ',' in opts.freq:
    freqs = np.sort(map(float, opts.freq.split(',')))
elif '-' in opts.freq:
    freqs_arr = np.sort(map(float, opts.freq.split('-')))
    freqs = np.round(np.arange(freqs_arr[0], freqs_arr[1] + opts.freq_res, opts.freq_res),
                              decimals=3)
    freqs = freqs[np.where(freqs <= freqs_arr[1])]
elif not (',' in opts.freq and '-' in opts.freq):
    freqs = [float(opts.freq)]

nfreqs = len(freqs)


## -------------------------------------- Construct sky -------------------------------------- ##
print 'Constructing sky...'

# Construct l,m grid
FOV = np.deg2rad(40)
ls = np.linspace(-FOV/2, FOV/2, NPIX_SIDE)
ms = np.copy(ls)
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

grid_pos = np.zeros((nsources, 2), dtype = int)
true_pos = np.zeros((nsources, 2))

if opts.zenith_source:
    for i in range(grid_pos.shape[0]):
        grid_pos[i, 0] = mid_l
        grid_pos[i, 1] = mid_m
        true_pos[i, 0] = ls[grid_pos[i, 0]] + np.random.uniform(low=-lm_pixel_half,
                                                                                          high=lm_pixel_half)
        true_pos[i, 1] = ms[grid_pos[i, 1]] + np.random.uniform(low=-lm_pixel_half,
                                                                                            high=lm_pixel_half)

elif opts.horizon_source:
    for i in range(grid_pos.shape[0]):
        grid_pos[i, 0] = mid_l
        grid_pos[i, 1] = 0
        true_pos[i, 0] = ls[grid_pos[i, 0]] + np.random.uniform(low=-lm_pixel_half,
                                                                                          high=lm_pixel_half)
        true_pos[i, 1] = ms[grid_pos[i, 1]] + np.random.uniform(low=0,
                                                                                            high=lm_pixel_half)

elif opts.l_offset or opts.m_offset:
    if opts.l_offset:
        l_off = int(mid_l*opts.l_offset)
    else:
        l_off = 0

    if opts.m_offset:
        m_off = int(mid_m*opts.m_offset)
    else:
        m_off = 0

    for i in range(grid_pos.shape[0]):
        grid_pos[i, 0] = mid_l + l_off
        grid_pos[i, 1] = mid_m + m_off
        true_pos[i, 0] = ls[grid_pos[i, 0]] + np.random.uniform(low=-lm_pixel_half,
                                                                                          high=lm_pixel_half)
        true_pos[i, 1] = ms[grid_pos[i, 1]] + np.random.uniform(low=-lm_pixel_half,
                                                                                            high=lm_pixel_half)

else:
    for i in range(grid_pos.shape[0]):
        grid_pos[i, 0] = np.random.randint(0, ls.shape[0])
        grid_pos[i, 1] = np.random.randint(0, ms.shape[0])
        true_pos[i, 0] = ls[grid_pos[i, 0]] + np.random.uniform(low=-lm_pixel_half,
                                                                                          high=lm_pixel_half)
        true_pos[i, 1] = ms[grid_pos[i, 1]] + np.random.uniform(low=-lm_pixel_half,
                                                                                            high=lm_pixel_half)

# Make sky matrix
Sky = np.zeros((nfreqs, NPIX_SIDE, NPIX_SIDE))
Sky_counts = np.zeros((NPIX_SIDE, NPIX_SIDE))
for i, freq in enumerate(freqs):
    if not opts.spec_index == 0.0:
        Sky[i, grid_pos[:, 1], grid_pos[:, 0]] += 1./(1 + (freq - freqs.min())**opts.spec_index)
    else:
        Sky[i, grid_pos[:, 1], grid_pos[:, 0]] += 1.

unique_grid_pos, unique_grid_pos_counts = np.unique(grid_pos,
                                                                                 axis=0,
                                                                                 return_counts=True)
Sky_counts[unique_grid_pos[:, 1], unique_grid_pos[:, 0]] = unique_grid_pos_counts

# Flatten sky for matrix computation
Sky_vec = Sky.flatten()

if opts.beam:
    # Get beam on sky grid
    thetas = np.sqrt(ls_vec**2 + ms_vec**2)
    thetas *= FOV/(2.*thetas.max())
    phis = np.arctan2(np.copy(ms_vec), np.copy(ls_vec))
    phis[phis < 0.0] += 2*np.pi

    # Get beam on my (l, m) grid
    beam_E = fits.getdata(opts.beam, extname='BEAM_E')
    beam_freqs = fits.getdata(opts.beam, extname='FREQS')
    beam_grid = np.zeros((nfreqs, thetas.size))

    # Get beam on grid
    for i, freq in enumerate(freqs):
        freq_ind = np.where(beam_freqs == freq)[0][0]
        beam_grid[i] = hp.get_interp_val(beam_E[:, freq_ind], thetas, phis)


## ----------------------------------- Construct Visibilities ----------------------------------- ##
print 'Constructing visibilities...'

# Construct u,v grid
us_grid = np.fft.fftshift(np.fft.fftfreq(ls.shape[0], d=np.mean(np.diff(ls))))
vs_grid = np.fft.fftshift(np.fft.fftfreq(ms.shape[0], d=np.mean(np.diff(ms))))
us_vec, vs_vec = np.zeros(0), np.zeros(0)
for v in vs_grid:
    for u in us_grid:
        us_vec = np.append(us_vec, u)
        vs_vec = np.append(vs_vec, v)
uv_pixel_half = np.mean(np.diff(us_grid))/2.

# Use analytical solution to get visibilities using true positions
Vs = np.zeros((nfreqs, NPIX), dtype=complex)
for i in range(nfreqs):
    for j in range(us_vec.size):
        if opts.grid_pos:
            if opts.beam:
                Vs[i, j] = beam_grid[i, j]*np.sum(Vs_func(us_vec[j],
                                                       ls[grid_pos[:, 0]],
                                                       vs_vec[j],
                                                       ms[grid_pos[:, 1]]))
            else:
                Vs[i, j] = np.sum(Vs_func(us_vec[j],
                                                       ls[grid_pos[:, 0]],
                                                       vs_vec[j],
                                                       ms[grid_pos[:, 1]]))
        else:
            if opts.beam:
                Vs[i, j] = beam_grid[i, j]*np.sum(Vs_func(us_vec[j],
                                                                              true_pos[:, 0],
                                                                              vs_vec[j],
                                                                              true_pos[:, 1]))
            else:
                Vs[i, j] = np.sum(Vs_func(us_vec[j],
                                                       true_pos[:, 0],
                                                       vs_vec[j],
                                                       true_pos[:, 1]))

Vs = Vs.reshape((nfreqs, NPIX_SIDE, NPIX_SIDE))

## ---------------------------------- Construct MaxL Sky ---------------------------------- ##
print 'Constructing maximum likelihood sky...'
# Construct solution using analytic solution for maximum likelihood
# Assumes a Gaussian log likelihood function
# Requires noise injection into data (visibilities above)
a = np.zeros_like(Vs)
Vs_maxL = np.zeros_like(a)
N_inv = np.eye(NPIX)/opts.rms**2

# Create data from visibilities with injected Gaussian noise
d = np.copy(Vs)

# Iterate over frequencies
print 'Entering for loop...'
for i in range(nfreqs):
    # Noise must be added so that d is still Hermitian
    d_flat = d[i].flatten()
    Vs_flat = Vs[i].flatten()

    half_ind = int(us_vec.size/2.) + 1
    for j, [u,v] in enumerate(np.stack((us_vec, vs_vec), axis=1)[:half_ind]):
        neg_ind = np.where(np.logical_and(us_vec == -u, vs_vec == -v))[0][0]
        d_flat[[j, neg_ind]] += (np.random.normal(0, opts.rms, 1)
                                          +
                                          1j*np.random.normal(0, opts.rms, 1))

    d[i] = d_flat.reshape([NPIX_SIDE]*2)

    DFT = np.exp(-1j*2*np.pi*(np.outer(us_vec, ls_vec)
                        +
                        np.outer(vs_vec, ms_vec)))

    if opts.beam:
        freq_ind = np.where(beam_freqs == freqs[i])
        P = np.diag(beam_E[:, freq_ind])
        inv_part = np.linalg.inv(np.dot(np.dot(np.dot(np.dot(P, DFT.conj().T), N_inv), DFT), P))
        right_part = np.dot(np.dot(np.dot(P, DFT.conj().T), N_inv), d[i].flatten())
    else:
        inv_part = np.linalg.inv(np.dot(np.dot(DFT.conj().T, N_inv), DFT))
        right_part = np.dot(np.dot(DFT.conj().T, N_inv), d[i].flatten())

    # Maximum likelihood solution for the sky
    a[i] = np.dot(inv_part, right_part).reshape((NPIX_SIDE, NPIX_SIDE))

    # Generate visibilities from maximum liklihood solution
    Vs_maxL[i] = np.dot(DFT, a[i].flatten()).reshape((NPIX_SIDE, NPIX_SIDE))

    if opts.beam:
        del(DFT, inv_part, right_part, P)
    else:
        del(DFT, inv_part, right_part)

print 'For loop finished...'

# # Fit for RMS of residuals
# diff_data = np.abs(Vs) - np.abs(Vs_maxL)
# counts, bins = np.histogram(diff_data[0].flatten(), bins=50)
# bin_width = np.mean(np.diff(bins))
# fit_xs = bins[:-1] + bin_width/2
# guess_params = [np.max(counts), 0.0, opts.rms]
# fit_params, fit_cov = curve_fit(Gaussian, fit_xs, counts, p0=guess_params)
# # 0: amplitude, 1: mean, 2:std dev
# gauss_fit = Gaussian(fit_xs, fit_params[0], fit_params[1], fit_params[2])

if opts.write:
    # Write fitted RMS data
    if os.path.exists('./sim_vis/'):
        if nfreqs > 1:
            filename = 'sim_vis/maxL_visdata_%sMHz_%sMHz' %(opts.freq, opts.freq_res)
        else:
            filename = 'sim_vis/maxL_visdata_%sMHz' %opts.freq
    else:
        if nfreqs > 1:
            filename = 'maxL_visdata_%sMHz_%sMHz' %(opts.freq, opts.freq_res)
        else:
            filename = 'maxL_visdata_%sMHz' %opts.freq
    print 'Writing ' + filename + '.npy ...\n'
    out_dic = {}
    out_dic['vis'] = Vs
    out_dic['maxL_sky'] = a
    out_dic['maxL_vis'] = Vs_maxL
    out_dic['input_rms'] = opts.rms
    np.save(filename + '.npy', out_dic)

    sys.exit()


## ---------------------------------- Plotting ---------------------------------- ##
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib.pyplot import *

freq_ind = 0
fontsize = 14
extent_uv = [us_grid.min() - uv_pixel_half, us_grid.max() + uv_pixel_half,
                    vs_grid.min() - uv_pixel_half, vs_grid.max() + uv_pixel_half]
extent_lm = [ls.min() - lm_pixel_half, ls.max() + lm_pixel_half,
                    ms.min() - lm_pixel_half, ms.max() + lm_pixel_half]

fig = figure(figsize=(16,8))
gs = gridspec.GridSpec(2, 3)


# Plot sky
skyax = fig.add_subplot(gs[0,0])
skyax.scatter(true_pos[:, 0], true_pos[:, 1], marker='.', c='r', alpha=0.5, s=25)
if opts.log_scale:
    skyim = skyax.imshow(np.log10(Sky[freq_ind]*Sky_counts),
                                      extent=extent_lm,
                                      origin='lower')
    skyax.set_title('Log Sky, %.1fMHz' %freqs[freq_ind], fontsize=fontsize)
else:
    skyim = skyax.imshow(Sky[freq_ind]*Sky_counts,
                                      extent=extent_lm,
                                      origin='lower')
    skyax.set_title('Sky, %.1fMHz' %freqs[freq_ind], fontsize=fontsize)
skyax.set_xlabel('l')
skyax.set_ylabel('m')


# Plot visibilities using true source position
visax = fig.add_subplot(gs[0,1])
if np.abs(Vs[freq_ind]).min() == np.abs(Vs[freq_ind]).max():
    if opts.log_scale:
        vmin = np.log10(np.abs(Vs[freq_ind]).min() - 1.e-5)
        vmax = np.log10(np.abs(Vs[freq_ind]).min() + 1.e-5)
    else:
        vmin = np.abs(Vs[freq_ind]).min() - 1.e-5
        vmax = np.abs(Vs[freq_ind]).min() + 1.e-5
else:
    if opts.log_scale:
        vmin = np.log10(np.abs(Vs[freq_ind]).min())
        vmax = np.log10(np.abs(Vs[freq_ind]).max())
    else:
        vmin = np.abs(Vs[freq_ind]).min()
        vmax = np.abs(Vs[freq_ind]).max()
if opts.log_scale:
    visim = visax.imshow(np.log10(np.abs(Vs[freq_ind])),
                                    extent=extent_uv,
                                    vmin=vmin, vmax=vmax,
                                    origin='lower')
    visax.set_title('Log |Vs|, %.1fMHz ' %freqs[freq_ind], fontsize=fontsize)
else:
    visim = visax.imshow(np.abs(Vs[freq_ind]),
                                    extent=extent_uv,
                                    vmin=vmin, vmax=vmax,
                                    origin='lower')
    visax.set_title('|Vs|, %.1fMHz ' %freqs[freq_ind], fontsize=fontsize)
visax.set_xlabel('u [m]')
visax.set_ylabel('v [m]')


# Plot sky solution for maximum likelihood
anskyax = fig.add_subplot(gs[1,0])
if opts.log_scale:
    anskyim = anskyax.imshow(np.log10(np.real(a[freq_ind])),
                                             extent=extent_lm,
                                             origin='lower')
    anskyax.set_title('Log MaxL Sky, %.1fMHz' %freqs[freq_ind], fontsize=fontsize)
else:
    anskyim = anskyax.imshow(np.real(a[freq_ind]),
                                             extent=extent_lm,
                                             origin='lower')
    anskyax.set_title('MaxL Sky, %.1fMHz' %freqs[freq_ind], fontsize=fontsize)
anskyax.set_xlabel('l')
anskyax.set_ylabel('m')


# Plot visibilities for maximum likelihood solution
anvisax = fig.add_subplot(gs[1,1])
if opts.log_scale:
    anvisim = anvisax.imshow(np.log10(np.abs(Vs_maxL[freq_ind])),
                                           extent=extent_uv,
                                           origin='lower')
    anvisax.set_title('Log|MaxL Visibilities|, %.1fMHz' %freqs[freq_ind], fontsize=fontsize)
else:
    anvisim = anvisax.imshow(np.abs(Vs_maxL[freq_ind]),
                                           extent=extent_uv,
                                           origin='lower')
    anvisax.set_title('|MaxL Visibilities|, %.1fMHz' %freqs[freq_ind], fontsize=fontsize)
anvisax.set_xlabel('u [m]')
anvisax.set_ylabel('v [m]')


# Plot |Vs| - |Vs_maxL|
diffax = fig.add_subplot(gs[:, -1])
if opts.log_scale:
    diffim = diffax.imshow(np.log10(diff_data[freq_ind]),
                                      extent=extent_uv,
                                      origin='lower')
    diffax.set_title('Log(|Vs| - |MaxL Vs|), %.1fMHz\nFitted RMS: %.2e'
                          %(freqs[freq_ind], fit_params[2]), fontsize=fontsize)
else:
    diffim = diffax.imshow(diff_data[freq_ind],
                                      extent=extent_uv,
                                      origin='lower')
    diffax.set_title('|Vs| - |MaxL Vs|, %.1fMHz\nFitted RMS: %.2e'
                          %(freqs[freq_ind], fit_params[2]), fontsize=fontsize)
diffax.set_xlabel('u [m]')
diffax.set_ylabel('v [m]')


# imgs = [skyim, visim, ftvisim, anskyim, anvisim, anftvisim]
imgs = [skyim, visim, anskyim, anvisim, diffim]

for i,ax in enumerate(fig.axes):
    ax_divider = make_axes_locatable(ax)
    cax = ax_divider.append_axes("right", size="5%", pad="2%")
    cb = fig.colorbar(imgs[i], cax=cax)

gs.tight_layout(fig)
show()
