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

o.add_option('--fit_beam',
    action = 'store_true',
    help = 'If passed, fit a two dimensional Gaussian to the beam for use in the maxL solution.')

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

def twoD_Gaussian((x, y), amp, x0, y0, sigma_x, sigma_y):
    # See https://stackoverflow.com/questions/21566379/
    # fitting-a-2d-gaussian-function-using-scipy-optimize-curve-fit-valueerror-and-m
    # for example of 2d Guassian fitting using scipy.optimize.curve_fit
    return (amp*np.exp(-(x-x0)**2/(2*sigma_x**2) - (y-y0)**2/(2*sigma_y**2))).ravel()


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
    freqs = [float(opts.freq)]

nfreqs = len(freqs)


## -------------------------------------- Construct sky -------------------------------------- ##
print 'Constructing sky...'

# Construct l,m grid
FOV = np.deg2rad(opts.fov)
ls = np.linspace(-FOV/2, FOV/2, npix_side)
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
        grid_pos[i, 0] = mid_m
        grid_pos[i, 1] = mid_l
        true_pos[i, 0] = ls[grid_pos[i, 1]] + np.random.uniform(low=-lm_pixel_half,
                                                                                          high=lm_pixel_half)
        true_pos[i, 1] = ms[grid_pos[i, 0]] + np.random.uniform(low=-lm_pixel_half,
                                                                                            high=lm_pixel_half)

elif opts.horizon_source:
    for i in range(grid_pos.shape[0]):
        grid_pos[i, 0] = 0
        grid_pos[i, 1] = mid_l
        true_pos[i, 0] = ls[grid_pos[i, 1]] + np.random.uniform(low=-lm_pixel_half,
                                                                                          high=lm_pixel_half)
        true_pos[i, 1] = ms[grid_pos[i, 0]] + np.random.uniform(low=0,
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
        grid_pos[i, 0] = mid_m + m_off
        grid_pos[i, 1] = mid_l + l_off
        true_pos[i, 0] = ls[grid_pos[i, 1]] + np.random.uniform(low=-lm_pixel_half,
                                                                                          high=lm_pixel_half)
        true_pos[i, 1] = ms[grid_pos[i, 0]] + np.random.uniform(low=-lm_pixel_half,
                                                                                            high=lm_pixel_half)

else:
    for i in range(grid_pos.shape[0]):
        grid_pos[i, 0] = np.random.randint(0, ms.shape[0])
        grid_pos[i, 1] = np.random.randint(0, ls.shape[0])
        true_pos[i, 0] = ls[grid_pos[i, 1]] + np.random.uniform(low=-lm_pixel_half,
                                                                                          high=lm_pixel_half)
        true_pos[i, 1] = ms[grid_pos[i, 0]] + np.random.uniform(low=-lm_pixel_half,
                                                                                            high=lm_pixel_half)

if opts.beam:
    if nsources == 1:
        pos_vec = np.where(np.logical_and(ls_vec == ls[grid_pos[:, 1]],
                                                            ms_vec == ms[grid_pos[:, 0]]))[0][0]
    else:
        pos_vec = np.zeros(0, dtype=int)
        for (l, m) in np.stack([ls[grid_pos[:, 1]], ms[grid_pos[:, 0]]], axis=1):
            pos_vec = np.append(pos_vec, np.where(np.logical_and(ls_vec == l,
                                                                                              ms_vec == m))[0][0])

# Make sky matrix
Sky = np.zeros((nfreqs, npix_side, npix_side))
Sky_counts = np.zeros((npix_side, npix_side))
for i, freq in enumerate(freqs):
    if not opts.spec_index == 0.0:
        Sky[i, grid_pos[:, 0], grid_pos[:, 1]] += 1./(1 + (freq - freqs.min())**opts.spec_index)
    else:
        Sky[i, grid_pos[:, 0], grid_pos[:, 1]] += 1.

unique_grid_pos, unique_grid_pos_counts = np.unique(grid_pos,
                                                                                 axis=0,
                                                                                 return_counts=True)
Sky_counts[unique_grid_pos[:, 0], unique_grid_pos[:, 1]] = unique_grid_pos_counts

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

    if opts.fit_beam:
        L, M = np.meshgrid(ls, ms)
        fit_beam_grid = np.zeros_like(beam_grid)

    # Get beam on grid
    for i, freq in enumerate(freqs):
        freq_ind = np.where(beam_freqs == freq)[0][0]
        beam_grid[i] = hp.get_interp_val(beam_E[:, freq_ind], thetas, phis)

        if opts.fit_beam:
            initial_guess = (beam_grid[i].max(), 0., 0., 1., 1.)
            popt, pcov = curve_fit(twoD_Gaussian, (L, M), beam_grid[i], p0=initial_guess)
            fit_beam_grid[i] = twoD_Gaussian((L, M), *popt)

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
for i in range(nfreqs):
    for j in range(us_vec.size):
        if opts.uvdata:
            inv_wavelength = freqs[i]*1.e6/C
        else:
            inv_wavelength = 1.
        if opts.grid_pos:
            if opts.beam:
                Vs[i, j] = np.sum(beam_grid[i, pos_vec]*Vs_func(us_vec[j]*inv_wavelength,
                                                       ls[grid_pos[:, 1]],
                                                       vs_vec[j]*inv_wavelength,
                                                       ms[grid_pos[:, 0]]))
            else:
                Vs[i, j] = np.sum(Vs_func(us_vec[j]*inv_wavelength,
                                                       ls[grid_pos[:, 1]],
                                                       vs_vec[j]*inv_wavelength,
                                                       ms[grid_pos[:, 0]]))
        else:
            if opts.beam:
                Vs[i, j] = np.sum(beam_grid[i, pos_vec]*Vs_func(us_vec[j]*inv_wavelength,
                                                                              true_pos[:, 0],
                                                                              vs_vec[j]*inv_wavelength,
                                                                              true_pos[:, 1]))
            else:
                Vs[i, j] = np.sum(Vs_func(us_vec[j]*inv_wavelength,
                                                       true_pos[:, 0],
                                                       vs_vec[j]*inv_wavelength,
                                                       true_pos[:, 1]))


## ---------------------------------- Construct MaxL Sky ---------------------------------- ##
print 'Constructing maximum likelihood sky...'
# Construct solution using analytic solution for maximum likelihood
# Assumes a Gaussian log likelihood function
# Requires noise injection into data (visibilities above)
a = np.zeros_like(Vs)
Vs_maxL = np.zeros_like(a)
N_inv = np.eye(npix)/opts.rms**2

# Create data from visibilities with injected Gaussian noise
d = np.copy(Vs)

# Iterate over frequencies
print 'Frequency [MHz]: ',
for freq_ind in range(nfreqs):
    if freq_ind == nfreqs - 1:
        print '%.0f' %freqs[freq_ind]
    else:
        print '%.0f, ' %freqs[freq_ind] ,
    # Noise must be added so that d is still Hermitian
    half_ind = int(us_vec.size/2.) + 1
    for j, [u,v] in enumerate(np.stack((us_vec, vs_vec), axis=1)[:half_ind]):
        neg_ind = np.where(np.logical_and(us_vec == -u, vs_vec == -v))[0][0]
        d[freq_ind, [j, neg_ind]] += (np.random.normal(0, opts.rms, 1)
                                                  +
                                                  1j*np.random.normal(0, opts.rms, 1))

    # Construct DFT matrix
    DFT = np.exp(-1j*2*np.pi*(np.outer(us_vec, ls_vec)
                        +
                        np.outer(vs_vec, ms_vec)))

    if opts.beam:
        if opts.fit_beam:
            P = np.diag(fit_beam_grid[freq_ind])
        else:
            P = np.diag(beam_grid[freq_ind])
        DftP = np.dot(DFT, P)
        inv_part = np.linalg.inv(np.dot(np.dot(DftP.conj().T, N_inv), DftP))
        right_part = np.dot(np.dot(DftP.conj().T, N_inv), d[freq_ind])
    else:
        inv_part = np.linalg.inv(np.dot(np.dot(DFT.conj().T, N_inv), DFT))
        right_part = np.dot(np.dot(DFT.conj().T, N_inv), d[freq_ind])

    # Maximum likelihood solution for the sky
    a[freq_ind] = np.dot(inv_part, right_part)


# Fit quadratics in frequency along frequency axis if a
fit_coeffs = np.polyfit(freqs, np.abs(a), 2) # shape (3, npix), i.e. 3 coefficients for each pixel
rms_data = np.zeros_like(a[0])
quad_fits = np.zeros_like(a)

for pix_ind in range(npix):
    quad_fits[:, pix_ind] = (fit_coeffs[0, pix_ind]*freqs**2
                                       +
                                       fit_coeffs[1, pix_ind]*freqs
                                       +
                                       fit_coeffs[2, pix_ind])

residuals = np.abs(a) - np.abs(quad_fits)
rms_data = np.std(residuals, axis=0)

if opts.write:
    # Write fitted RMS data
    if os.path.exists('./sim_vis/'):
        if nfreqs > 1:
            filename = 'sim_vis/maxL_rms_freq_fit_%sMHz_%sMHz_%.0fdfov' %(opts.freq,
                                                                                                                  opts.freq_res,
                                                                                                                  np.rad2deg(FOV))
        else:
            filename = 'sim_vis/maxL_rms_freq_fit_%sMHz_%.0fdfov' %(opts.freq,
                                                                                                      np.rad2deg(FOV))
    else:
        if nfreqs > 1:
            filename = 'maxL_rms_freq_fit_%sMHz_%sMHz_%.0fdfov' %(opts.freq,
                                                                                                      opts.freq_res,
                                                                                                      np.rad2deg(FOV))
        else:
            filename = 'maxL_rms_freq_fit_%sMHz_%.0fdfov' %(opts.freq,
                                                                                          np.rad2deg(FOV))
    print 'Writing ' + filename + '.npy ...\n'
    
    out_dic = {}
    out_dic['sky'] = Sky[freq_ind]*Sky_counts
    out_dic['vis'] = Vs
    out_dic['maxL_sky'] = a
    out_dic['maxL_vis'] = Vs_maxL
    out_dic['input_rms'] = opts.rms
    out_dic['freqs'] = freqs
    if opts.beam:
        out_dic['beam_file'] = opts.beam
    if opts.fit_beam:
        out_dic['fitted_beam'] = True
    else:
        out_dic['fitted_beam'] = False
    out_dic['rms_data'] = rms_data
    np.save(filename + '.npy', out_dic)

    sys.exit()


## ---------------------------------- Plotting ---------------------------------- ##
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib.pyplot import *

freq_ind = 0
fontsize = 14
extent_lm = [ls.min() - lm_pixel_half, ls.max() + lm_pixel_half,
                    ms.min() - lm_pixel_half, ms.max() + lm_pixel_half]
extent_lm = np.rad2deg(extent_lm)

fig = figure()
gs = gridspec.GridSpec(1, 1)


# Plot sky
rms_ax = fig.add_subplot(gs[0])
rms_im = imshow(rms_data.reshape([npix_side]*2))
rms_ax.set_xlabel('l')
rms_ax.set_ylabel('m')

imgs = [rms_im]

for i,ax in enumerate(fig.axes):
    ax_divider = make_axes_locatable(ax)
    cax = ax_divider.append_axes("right", size="5%", pad="2%")
    cb = fig.colorbar(imgs[i], cax=cax)
    cb.set_label('Fitted RMS/$10^{%.0f}$' %np.log10(opts.rms), size=16)

gs.tight_layout(fig)
show()
