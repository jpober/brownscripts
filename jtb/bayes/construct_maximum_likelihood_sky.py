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

o.add_option('--beam',
    type = 'str',
    # default = ('/Users/jburba/hera_things/hera-team/HERA-Beams/'
    #                 +
    #                 'NicolasFagnoniBeams/healpix_beam.fits'),
    help = 'Healpix compatible fits file containing the primary beam.')

o.add_option('--fit_beam',
    action = 'store_true',
    help = 'If passed, fit a two dimensional Gaussian to the beam for use in the maxL solution.')

o.add_option('--fractional_fit',
    action = 'store_true',
    help = 'If passed, normalize the beam at each frequency to have a peak at 1.')

o.add_option('--ref_freq',
    type = 'int',
    default = 150,
    help = 'Frequency to use as reference frequency for external fitting of frequency structure.')

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
print o.values

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
#
# if opts.zenith_source:
#     for i in range(grid_pos.shape[0]):
#         grid_pos[i, 0] = mid_m
#         grid_pos[i, 1] = mid_l
#         true_pos[i, 0] = ls[grid_pos[i, 1]] + np.random.uniform(low=-lm_pixel_half,
#                                                                                           high=lm_pixel_half)
#         true_pos[i, 1] = ms[grid_pos[i, 0]] + np.random.uniform(low=-lm_pixel_half,
#                                                                                             high=lm_pixel_half)
#
# elif opts.horizon_source:
#     for i in range(grid_pos.shape[0]):
#         grid_pos[i, 0] = 0
#         grid_pos[i, 1] = mid_l
#         true_pos[i, 0] = ls[grid_pos[i, 1]] + np.random.uniform(low=-lm_pixel_half,
#                                                                                           high=lm_pixel_half)
#         true_pos[i, 1] = ms[grid_pos[i, 0]] + np.random.uniform(low=0,
#                                                                                             high=lm_pixel_half)
#
# elif opts.l_offset or opts.m_offset:
#     if opts.l_offset:
#         l_off = int(mid_l*opts.l_offset)
#     else:
#         l_off = 0
#
#     if opts.m_offset:
#         m_off = int(mid_m*opts.m_offset)
#     else:
#         m_off = 0
#
#     for i in range(grid_pos.shape[0]):
#         grid_pos[i, 0] = mid_m + m_off
#         grid_pos[i, 1] = mid_l + l_off
#         true_pos[i, 0] = ls[grid_pos[i, 1]] + np.random.uniform(low=-lm_pixel_half,
#                                                                                           high=lm_pixel_half)
#         true_pos[i, 1] = ms[grid_pos[i, 0]] + np.random.uniform(low=-lm_pixel_half,
#                                                                                             high=lm_pixel_half)
#
# elif opts.uniform_sky:
#     nsources = nlm
#     grid_pos = np.array([[0, 0]])
#     for i in range(0, opts.npix_side):
#         for j in range(0, opts.npix_side):
#             grid_pos = np.append(grid_pos, np.array([[i, j]]), axis=0)
#     grid_pos = grid_pos[1:]
#     true_pos = np.stack([ms_vec, ls_vec]).T
#
# else:
#     for i in range(grid_pos.shape[0]):
#         grid_pos[i, 0] = np.random.randint(0, ms.shape[0])
#         grid_pos[i, 1] = np.random.randint(0, ls.shape[0])
#         true_pos[i, 0] = ls[grid_pos[i, 1]] + np.random.uniform(low=-lm_pixel_half,
#                                                                                           high=lm_pixel_half)
#         true_pos[i, 1] = ms[grid_pos[i, 0]] + np.random.uniform(low=-lm_pixel_half,
#                                                                                             high=lm_pixel_half)
#
# if opts.beam:
#     if nsources == 1:
#         pos_vec = np.where(np.logical_and(ls_vec == ls[grid_pos[:, 1]],
#                                                             ms_vec == ms[grid_pos[:, 0]]))[0][0]
#     else:
#         pos_vec = np.zeros(0, dtype=int)
#         for (l, m) in np.stack([ls[grid_pos[:, 1]], ms[grid_pos[:, 0]]], axis=1):
#             pos_vec = np.append(pos_vec, np.where(np.logical_and(ls_vec == l,
#                                                                                               ms_vec == m))[0][0])
#
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
    beam_E = fits.getdata(opts.beam, extname='BEAM_E')
    beam_freqs = fits.getdata(opts.beam, extname='FREQS')
    beam_grid = np.zeros((nfreqs, thetas.size))

    if opts.fit_beam:
        L, M = np.meshgrid(ls, ms)
        fit_beam_grid = np.zeros_like(beam_grid)

    # Get beam on grid
    ref_freq_ind = np.where(beam_freqs == opts.ref_freq)[0][0]
    ref_beam = hp.get_interp_val(beam_E[:, ref_freq_ind], thetas, phis)
    ref_beam /= ref_beam.max()
    for i, freq in enumerate(freqs):
        freq_ind = np.where(beam_freqs == freq)[0][0]
        beam_grid[i] = hp.get_interp_val(beam_E[:, freq_ind], thetas, phis)
        beam_grid[i] /= beam_grid[i].max()

        if opts.fractional_fit:
            beam_grid[i] /= ref_beam

        if opts.fit_beam:
            initial_guess = (0., 0., 1., 1.)
            # popt, pcov = curve_fit(lambda (l, m), l0, m0, sl, sm: twoD_Gaussian((l, m), l0, m0, sl, sm), (L, M), beam_grid[i], p0=initial_guess)
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
# if opts.uniform_sky:
print 'Constructing visibilities using FFT...'
for i in range(nfreqs):
    # Construct visibilities faster as FT of beam*sky
    if opts.beam:
        Vs[i] = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift((beam_grid[i]*Sky_vec[i]).reshape([npix_side]*2)))).flatten()
    else:
        Vs[i] = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift((Sky_vec[i]).reshape([npix_side]*2)))).flatten()
# else:
#     for i in range(nfreqs):
#         for j in range(us_vec.size):
#             if opts.uvdata:
#                 inv_wavelength = freqs[i]*1.e6/C
#             else:
#                 inv_wavelength = 1.
#             if opts.grid_pos:
#                 if opts.beam:
#                     Vs[i, j] = np.sum(beam_grid[i, pos_vec]*Vs_func(us_vec[j]*inv_wavelength,
#                                                            ls[grid_pos[:, 1]],
#                                                            vs_vec[j]*inv_wavelength,
#                                                            ms[grid_pos[:, 0]]))
#                 else:
#                     Vs[i, j] = np.sum(Vs_func(us_vec[j]*inv_wavelength,
#                                                            ls[grid_pos[:, 1]],
#                                                            vs_vec[j]*inv_wavelength,
#                                                            ms[grid_pos[:, 0]]))
#             else:
#                 if opts.beam:
#                     Vs[i, j] = np.sum(beam_grid[i, pos_vec]*Vs_func(us_vec[j]*inv_wavelength,
#                                                                                   true_pos[:, 0],
#                                                                                   vs_vec[j]*inv_wavelength,
#                                                                                   true_pos[:, 1]))
#                 else:
#                     Vs[i, j] = np.sum(Vs_func(us_vec[j]*inv_wavelength,
#                                                            true_pos[:, 0],
#                                                            vs_vec[j]*inv_wavelength,
#                                                            true_pos[:, 1]))


## ---------------------------------- Construct MaxL Sky ---------------------------------- ##
print 'Constructing maximum likelihood sky...'
# Construct solution using analytic solution for maximum likelihood
# Assumes a Gaussian log likelihood function
# Requires noise injection into data (visibilities above)
a = np.zeros_like(Vs)
# Vs_maxL = np.zeros_like(a)
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
        complex_noise = np.random.normal(0, opts.rms, 1) + 1j*np.random.normal(0, opts.rms, 1)
        d[freq_ind, j] += complex_noise
        d[freq_ind, neg_ind] += complex_noise.conjugate()

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

    # Generate visibilities from maximum liklihood solution
    # if opts.beam:
    #     Vs_maxL[freq_ind] = np.dot(DftP,  a[freq_ind])
    #
    #     # del(DFT, inv_part, right_part, P, DftP)
    # else:
    #     Vs_maxL[freq_ind] = np.dot(DFT, a[freq_ind])
    #
    #     # del(DFT, inv_part, right_part)

print 'For loop finished...'

if opts.write:
    # Write fitted RMS data
    if os.path.exists('./sim_vis/'):
        if nfreqs > 1:
            filename = 'sim_vis/maxL_visdata_%sMHz_%sMHz_%.0fdfov' %(opts.freq,
                                                                                                             opts.freq_res,
                                                                                                             np.rad2deg(FOV))
        else:
            filename = 'sim_vis/maxL_visdata_%sMHz_%.0fdfov' %(opts.freq,
                                                                                                  np.rad2deg(FOV))
    else:
        if nfreqs > 1:
            filename = 'maxL_visdata_%sMHz_%sMHz_%.0fdfov' %(opts.freq,
                                                                                                  opts.freq_res,
                                                                                                  np.rad2deg(FOV))
        else:
            filename = 'maxL_visdata_%sMHz_%.0fdfov' %(opts.freq,
                                                                                      np.rad2deg(FOV))

    filename += '_%dnpix-side' %opts.npix_side
    filename += '_rms%.0e' %opts.rms
    if opts.zenith_source:
        filename += '_zenith-source'
    elif opts.horizon_source:
        filename += '_horizon-source'
    elif opts.uniform_sky:
        filename += '_uniform-sky'
    elif opts.noise_sky:
        filename += '_noise-sky'
    else:
        if opts.l_offset:
            if opts.m_offset:
                filename += '_loff-%.3f_moff-%.3f' %(opts.l_offset, opts.m_offset)
            else:
                filename += '_loff-%.3f' %opts.l_offset
        elif opts.m_offset:
            filename += '_moff-%.3f' %opts.m_offset
        else:
            filename += '_%dsources' %opts.nsources

    if opts.fractional_fit:
        filename += '_fractional-fit'
    else:
        filename += '_absolute-fit'



    out_dic = {}
    out_dic['sky'] = Sky
    out_dic['vis'] = Vs
    out_dic['maxL_sky'] = a
    # out_dic['maxL_vis'] = Vs_maxL
    out_dic['input_rms'] = opts.rms
    out_dic['freqs'] = freqs
    if opts.beam:
        out_dic['beam_file'] = opts.beam
    if opts.fit_beam:
        out_dic['fitted_beam'] = True
    else:
        out_dic['fitted_beam'] = False

    # try:
    #     os.makedirs(filename)
    #     print 'Made directory ' + filename
    #     print 'Writing ' + filename + '/data_dic.npy ...\n'
    #     np.save(filename + '/data_dic.npy', out_dic)
    # except:
    #     print 'Writing ' + filename + '.npy ...\n'
    #     np.save(filename + '.npy', out_dic)
    print 'Writing ' + filename + '.npy ...\n'
    np.save(filename + '.npy', out_dic)


    sys.exit()


# Fit for RMS of residuals
diff_data = np.abs(Vs) - np.abs(Vs_maxL)
counts, bins = np.histogram(diff_data[0], bins=50)
bin_width = np.mean(np.diff(bins))
fit_xs = bins[:-1] + bin_width/2
guess_params = [np.max(counts), 0.0, opts.rms]
fit_params, fit_cov = curve_fit(Gaussian, fit_xs, counts, p0=guess_params)
# 0: amplitude, 1: mean, 2:std dev
gauss_fit = Gaussian(fit_xs, *fit_params)

## ---------------------------------- Plotting ---------------------------------- ##
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib.pyplot import *

freq_ind = 0
inv_wavelength = 1.e6*freqs[freq_ind]/C
fontsize = 14
extent_uv = np.array([us_grid.min() - uv_pixel_half, us_grid.max() + uv_pixel_half,
                    vs_grid.min() - uv_pixel_half, vs_grid.max() + uv_pixel_half])
extent_uv *= inv_wavelength
extent_lm = [ls.min() - lm_pixel_half, ls.max() + lm_pixel_half,
                    ms.min() - lm_pixel_half, ms.max() + lm_pixel_half]

fig = figure(figsize=(16,8))
gs = gridspec.GridSpec(2, 3)


# Plot sky
skyax = fig.add_subplot(gs[0,0])
skyax.scatter(true_pos[:, 0], true_pos[:, 1], marker='.', c='r', alpha=0.5, s=25)
if opts.log_scale:
    skyim = skyax.imshow(np.log10(Sky[freq_ind]),
                                      extent=extent_lm,
                                      origin='lower')
    skyax.set_title('Log Sky, %.1fMHz' %freqs[freq_ind], fontsize=fontsize)
else:
    skyim = skyax.imshow(Sky[freq_ind],
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
    visim = visax.imshow(np.log10(np.abs(Vs[freq_ind])).reshape([npix_side]*2),
                                    extent=extent_uv,
                                    vmin=vmin, vmax=vmax,
                                    origin='lower')
    visax.set_title('Log |Vs|, %.1fMHz ' %freqs[freq_ind], fontsize=fontsize)
else:
    visim = visax.imshow(np.abs(Vs[freq_ind]).reshape([npix_side]*2),
                                    extent=extent_uv,
                                    vmin=vmin, vmax=vmax,
                                    origin='lower')
    visax.set_title('|Vs|, %.1fMHz ' %freqs[freq_ind], fontsize=fontsize)
visax.set_xlabel('u [wavelengths]')
visax.set_ylabel('v [wavelengths]')


# Plot sky solution for maximum likelihood
anskyax = fig.add_subplot(gs[1,0])
if opts.log_scale:
    anskyim = anskyax.imshow(np.log10(np.real(a[freq_ind])).reshape([npix_side]*2),
                                             extent=extent_lm,
                                             origin='lower')
    anskyax.set_title('Log MaxL Sky, %.1fMHz' %freqs[freq_ind], fontsize=fontsize)
else:
    anskyim = anskyax.imshow(np.real(a[freq_ind]).reshape([npix_side]*2),
                                             extent=extent_lm,
                                             origin='lower')
    anskyax.set_title('MaxL Sky, %.1fMHz' %freqs[freq_ind], fontsize=fontsize)
anskyax.set_xlabel('l')
anskyax.set_ylabel('m')


# Plot visibilities for maximum likelihood solution
anvisax = fig.add_subplot(gs[1,1])
if opts.log_scale:
    anvisim = anvisax.imshow(np.log10(np.abs(Vs_maxL[freq_ind])).reshape([npix_side]*2),
                                           extent=extent_uv,
                                           origin='lower')
    anvisax.set_title('Log|MaxL Visibilities|, %.1fMHz' %freqs[freq_ind], fontsize=fontsize)
else:
    anvisim = anvisax.imshow(np.abs(Vs_maxL[freq_ind]).reshape([npix_side]*2),
                                           extent=extent_uv,
                                           origin='lower')
    anvisax.set_title('|MaxL Visibilities|, %.1fMHz' %freqs[freq_ind], fontsize=fontsize)
anvisax.set_xlabel('u [wavelengths]')
anvisax.set_ylabel('v [wavelengths]')


# Plot |Vs| - |Vs_maxL|
diffax = fig.add_subplot(gs[:, -1])
if opts.log_scale:
    diffim = diffax.imshow(np.log10(diff_data[freq_ind]).reshape([npix_side]*2),
                                      extent=extent_uv,
                                      origin='lower')
    diffax.set_title('Log(|Vs| - |MaxL Vs|), %.1fMHz\nFitted RMS: %.2e\nMean: %.2e'
                          %(freqs[freq_ind], fit_params[2], fit_params[1]), fontsize=fontsize)
else:
    diffim = diffax.imshow(diff_data[freq_ind].reshape([npix_side]*2),
                                      extent=extent_uv,
                                      origin='lower')
    diffax.set_title('|Vs| - |MaxL Vs|, %.1fMHz\nFitted RMS: %.2e\nMean: %.2e'
                          %(freqs[freq_ind], fit_params[2], fit_params[1]), fontsize=fontsize)
diffax.set_xlabel('u [wavelengths]')
diffax.set_ylabel('v [wavelengths]')


# imgs = [skyim, visim, ftvisim, anskyim, anvisim, anftvisim]
imgs = [skyim, visim, anskyim, anvisim, diffim]

for i,ax in enumerate(fig.axes):
    ax_divider = make_axes_locatable(ax)
    cax = ax_divider.append_axes("right", size="5%", pad="2%")
    cb = fig.colorbar(imgs[i], cax=cax)

gs.tight_layout(fig)
show()
