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

if opts.data is '':
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
    L, M = np.meshgrid(ls, ms)
    ls_vec, ms_vec = L.flatten(), M.flatten()
    ls_vec = ls_vec.reshape((1, -1))
    ms_vec = ms_vec.reshape((1, -1))

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

        print '\tAdding source at l=%.1f, m=%.1f...' %(np.rad2deg(ls[mid_l + l_off]),
                                                                             np.rad2deg(ms[mid_m + m_off]))
        Sky[:, mid_m + m_off, mid_l + l_off] = 1.0

    elif opts.uniform_sky or opts.noise_sky:
        nsources = nlm
        if opts.noise_sky:
            Sky = np.ones_like(Sky)*np.random.normal(0, opts.rms, Sky.shape)
        else:
            Sky = np.ones_like(Sky)

    elif opts.gaussian_source:
        stddev = 0.03
        if opts.l_offset:
            l_off = int(mid_l*opts.l_offset)
        else:
            l_off = 0

        if opts.m_offset:
            m_off = int(mid_m*opts.m_offset)
        else:
            m_off = 0
        print '\tAdding Gaussian source at l=%.1f, m=%.1f...' %(np.rad2deg(ls[mid_l + l_off]),
                                                                                           np.rad2deg(ms[mid_m + m_off]))
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
        print 'Reading in beam file: %s' %opts.beam
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
            fit_beam_grid = np.zeros_like(beam_grid)

        # Get beam on grid
        ref_beam = None
        for i, freq in enumerate(freqs):
            freq_ind = np.where(beam_freqs == freq)[0][0]
            beam_grid[i] = hp.get_interp_val(beam_E[:, freq_ind], thetas, phis)
            beam_grid[i] = beam_grid[i]**2
            beam_grid[i] /= beam_grid[i].max()

            if opts.fractional_fit:
                if ref_beam is None:
                    print '\tSetting reference beam at %.3fMHz' %freq
                    ref_beam = np.copy(beam_grid[i])
                beam_grid[i] /= ref_beam

            if opts.fit_beam:
                initial_guess = (0., 0., 1., 1.)
                popt, pcov = curve_fit(twoD_Gaussian, (L, M), beam_grid[i], p0=initial_guess)
                print '\tFitting beam...'
                fit_beam_grid[i] = twoD_Gaussian((L, M), *popt)

        inset_fov = np.deg2rad(opts.inset_fov)
        inset_fov_l_inds = np.logical_and(ls_vec >= -inset_fov/2, ls_vec <= inset_fov/2)
        inset_fov_m_inds = np.logical_and(ms_vec >= -inset_fov/2, ms_vec <= inset_fov/2)
        inset_fov_inds = (inset_fov_l_inds*inset_fov_m_inds).flatten()
        ls_inset_vec = ls_vec[:, inset_fov_inds]
        ms_inset_vec = ms_vec[:, inset_fov_inds]
        ls_inset_vec = ls_inset_vec.reshape((1, -1))
        ms_inset_vec = ms_inset_vec.reshape((1, -1))
        ls_inset_grid = np.unique(ls_inset_vec)
        ms_inset_grid = np.copy(ls_inset_grid)

    ## ----------------------------------- Construct Visibilities ----------------------------------- ##
    print 'Constructing visibilities...'

    # Construct u,v grid
    us_inset_grid = np.fft.fftshift(np.fft.fftfreq(ls_inset_grid.shape[0], d=np.mean(np.diff(ls_inset_grid))))
    vs_inset_grid = np.fft.fftshift(np.fft.fftfreq(ms_inset_grid.shape[0], d=np.mean(np.diff(ms_inset_grid))))
    U_inset, V_inset = np.meshgrid(us_inset_grid, vs_inset_grid)
    us_inset_vec, vs_inset_vec = U_inset.flatten(), V_inset.flatten()
    us_inset_vec = us_inset_vec.reshape((-1, 1))
    vs_inset_vec = vs_inset_vec.reshape((-1, 1))
    uv_pixel_half = np.mean(np.diff(us_inset_grid))/2.
    nvis = us_inset_vec.shape[0]

    # Use analytical solution to get visibilities using true positions
    Vs = np.zeros((nfreqs, nvis), dtype=complex)
    # Use (u, v) corresponding to the inset fov to get "data"
    DFT = np.exp(-1j*2*np.pi*(np.dot(us_inset_vec, ls_vec)
                        +
                        np.dot(vs_inset_vec, ms_vec)))
    # if opts.uniform_sky:
    for i in range(nfreqs):
        if opts.beam:
            Vs[i] = np.dot(DFT, beam_grid[i]*Sky_vec[i])
        else:
            Vs[i] = np.dot(DFT, Sky_vec[i])



    ## ---------------------------------- Construct MaxL Sky ---------------------------------- ##
    print 'Constructing maximum likelihood sky...'
    # Construct solution using analytic solution for maximum likelihood
    # Assumes a Gaussian log likelihood function
    # Requires noise injection into data (visibilities above)
    a = np.zeros((nfreqs, ls_inset_vec.size), dtype=complex)
    N_inv = np.eye(us_inset_vec.size)/opts.rms**2

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
        half_ind = int(us_inset_vec.size/2.) + 1
        for j, [u,v] in enumerate(np.stack((us_inset_vec, vs_inset_vec), axis=1)[:half_ind]):
            neg_ind = np.where(np.logical_and(us_inset_vec == -u, vs_inset_vec == -v))[0][0]
            complex_noise = np.random.normal(0, opts.rms, 1) + 1j*np.random.normal(0, opts.rms, 1)
            d[freq_ind, j] += complex_noise
            d[freq_ind, neg_ind] += complex_noise.conjugate()

            # Construct DFT matrix
            DFT = np.exp(-1j*2*np.pi*(np.dot(us_inset_vec, ls_inset_vec)
                                +
                                np.dot(vs_inset_vec, ms_inset_vec)))

        if opts.beam:
            if opts.fit_beam:
                P = np.diag(fit_beam_grid[freq_ind, inset_fov_inds])
            else:
                P = np.diag(beam_grid[freq_ind, inset_fov_inds])
            DftP = np.dot(DFT, P)
            inv_part = np.linalg.inv(np.dot(np.dot(DftP.conj().T, N_inv), DftP))
            right_part = np.dot(np.dot(DftP.conj().T, N_inv), d[freq_ind])
        else:
            inv_part = np.linalg.inv(np.dot(np.dot(DFT.conj().T, N_inv), DFT))
            right_part = np.dot(np.dot(DFT.conj().T, N_inv), d[freq_ind])

        # Maximum likelihood solution for the sky
        a[freq_ind] = np.dot(inv_part, right_part)

    if opts.write:
        # Write fitted RMS data
        if nfreqs > 1:
            filename = 'ML_sky_%sMHz_%sMHz_%.0fdfov' %(opts.freq,
                                                                                      opts.freq_res,
                                                                                      np.rad2deg(FOV))
        else:
            filename = 'ML_sky_%sMHz_%.0fdfov' %(opts.freq,
                                                                         np.rad2deg(FOV))

        filename += '_%dnpix-side' %opts.npix_side

        if opts.zenith_source:
            filename += '_zenith-source'
        elif opts.horizon_source:
            filename += '_horizon-source'
        elif opts.uniform_sky:
            filename += '_uniform-sky'
        else:
            if opts.gaussian_source:
                filename += '_gaussian-source'
            if opts.l_offset:
                if opts.m_offset:
                    filename += '_loff-%.3f_moff-%.3f' %(opts.l_offset, opts.m_offset)
                else:
                    filename += '_loff-%.3f' %opts.l_offset
            elif opts.m_offset:
                filename += '_moff-%.3f' %opts.m_offset
            elif not (opts.l_offset or opts.m_offset or opts.gaussian_source):
                filename += '_%dsources' %opts.nsources

        if opts.fractional_fit:
            filename += '_fractional-fit'
        else:
            filename += '_absolute-fit'

        filename += '_%.0fdinset-fov' %opts.inset_fov

        print 'Writing ' + filename + '.npy ...\n'
        out_dic = {}
        out_dic['sky'] = Sky
        out_dic['inset_fov'] = opts.inset_fov
        out_dic['vis'] = Vs
        out_dic['maxL_sky'] = a
        out_dic['input_rms'] = opts.rms
        out_dic['freqs'] = freqs
        if opts.beam:
            out_dic['beam_file'] = beam_grid
        if opts.fit_beam:
            out_dic['fitted_beam'] = fit_beam_grid
        np.save(filename + '.npy', out_dic)

        sys.exit()

if not opts.data is '':
        print 'Reading in %s...' %opts.data

        # Read in rms data from opts.data
        data_dic = np.load(opts.data).item()
        a = data_dic['maxL_sky']
        Sky = data_dic['sky']
        beam_grid = data_dic['beam_file']
        fit_beam_grid = data_dic['fitted_beam']
        npix_ind = opts.data.find('npix')
        back_ind = opts.data.find('_', npix_ind - 4)
        npix_side = int(opts.data[back_ind+1:npix_ind])
        dfov_ind = opts.data.find('dfov')
        back_ind = opts.data.find('_', dfov_ind - 4)
        FOV = np.deg2rad(float(opts.data[back_ind + 1:dfov_ind]))
        ls = np.linspace(-FOV/2, FOV/2, npix_side)
        ms = np.copy(ls)
        nlm = ls.size*ms.size
        lm_pixel_half = np.diff(ls)[0]/2.
        ls_vec, ms_vec = np.zeros(0), np.zeros(0)
        for m in ms:
            for l in ls:
                ls_vec = np.append(ls_vec, l)
                ms_vec = np.append(ms_vec, m)
        if 'dinset' in opts.data:
            dinset_ind = opts.data.find('dinset')
            back_ind = opts.data.find('_', dinset_ind - 4)
            inset_fov = np.round(np.deg2rad(float(opts.data[back_ind + 1:dinset_ind])), decimals=14)
        inset_fov_l_inds = np.logical_and(ls_vec >= -inset_fov/2, ls_vec <= inset_fov/2)
        inset_fov_m_inds = np.logical_and(ms_vec >= -inset_fov/2, ms_vec <= inset_fov/2)
        inset_fov_inds = inset_fov_l_inds*inset_fov_m_inds
        ls_inset_grid = np.unique(ls_vec[inset_fov_inds])
        ms_inset_grid = np.copy(ls_inset_grid)
        ls_inset_grid = ls[np.logical_and(ls >= -inset_fov/2, ls <= inset_fov/2)]
        ms_inset_grid = np.copy(ls_inset_grid)
        if not any(x in opts.data for x in ['uniform', 'zenith', 'horizon']):
            nsources = np.sum(data_dic['sky'])
        nfreqs = 1
        Sky_vec = Sky.reshape((nfreqs, npix_side**2))
        us_grid = np.copy(ls_inset_grid)
        if inset_fov < FOV:
            diff_data = a[0].real - (Sky_vec[0]*beam_grid[0]/fit_beam_grid[0])[inset_fov_inds]
        else:
            diff_data = a[0].real - Sky_vec[0]*beam_grid[0]/fit_beam_grid[0]


## ---------------------------------- Plotting ---------------------------------- ##
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib.pyplot import *

print 'Plotting...'

nplots = 4
plot_size = 5
nrows = 1
ncols = nplots
fig = figure(figsize=(plot_size*(nplots + 1), plot_size))

fontsize = 14
extent_lm = [ls.min() - lm_pixel_half, ls.max() + lm_pixel_half,
                    ms.min() - lm_pixel_half, ms.max() + lm_pixel_half]
extent_lm_inset_fov = [ls_inset_grid.min() - lm_pixel_half, ls_inset_grid.max() + lm_pixel_half,
                    ms_inset_grid.min() - lm_pixel_half, ms_inset_grid.max() + lm_pixel_half]
extent_lm = np.rad2deg(extent_lm)
extent_lm_inset_fov = np.rad2deg(extent_lm_inset_fov)

gs = gridspec.GridSpec(nrows, nplots)
imgs = []

for plot_ind in range(nplots):
    rms_ax = fig.add_subplot(gs[plot_ind])

    if plot_ind == 0:
        if opts.fit_beam:
            if opts.log_scale:
                im = rms_ax.imshow(np.log10((Sky_vec[0]*beam_grid[0]/fit_beam_grid[0])[inset_fov_inds]).reshape([ls_inset_grid.size]*2),
                                                       origin='lower',
                                                       extent=extent_lm_inset_fov)
            else:
                im = rms_ax.imshow((Sky_vec[0]*beam_grid[0]/fit_beam_grid[0])[inset_fov_inds].reshape([ls_inset_grid.size]*2),
                                                origin = 'lower',
                                                extent = extent_lm_inset_fov)
                                                # vmin = a[0].real.min(),
                                                # vmax = a[0].real.max())
            rms_ax.set_title('Sky*Beam/Fitted_Beam')
        else:
            if opts.log_scale:
                im = rms_ax.imshow(np.log10(Sky_vec[0])[inset_fov_inds].reshape([ls_inset_grid.size]*2),
                                                       origin='lower',
                                                       extent=extent_lm_inset_fov)
            else:
                im = rms_ax.imshow((Sky_vec[0])[inset_fov_inds].reshape([ls_inset_grid.size]*2),
                                                origin = 'lower',
                                                extent = extent_lm_inset_fov)
                                                # vmin = a[0].real.min(),
                                                # vmax = a[0].real.max())
            rms_ax.set_title('Sky')
        # if inset_fov < FOV:
        #     vlines(np.rad2deg([ls_inset_grid.min(), ls_inset_grid.max()]),
        #                 np.rad2deg(ms_inset_grid.min()),
        #                 np.rad2deg(ms_inset_grid.max()))
        #     hlines(np.rad2deg([ms_inset_grid.min(), ms_inset_grid.max()]),
        #                 np.rad2deg(ls_inset_grid.min()),
        #                 np.rad2deg(ls_inset_grid.max()))
    else:
        if inset_fov < FOV:
            extent = extent_lm_inset_fov
        else:
            extent = extent_lm
        if plot_ind ==1:
            if opts.log_scale:
                im = rms_ax.imshow(np.log10(a[0].real).reshape([ls_inset_grid.size]*2),
                                                origin = 'lower',
                                                extent = extent)
            else:
                im = rms_ax.imshow((a[0].real).reshape([ls_inset_grid.size]*2),
                                                origin = 'lower',
                                                extent = extent)
            rms_ax.set_title('ML Sky Solution')
        elif plot_ind == 2:
            if not opts.fit_beam:
                rms_ax.set_title('ML Sky Solution - Sky')
                if opts.inset_fov < opts.fov:
                    diff_data = a[0].real - Sky_vec[0, inset_fov_inds]
                else:
                    diff_data = a[0].real - Sky_vec[0]
            else:
                rms_ax.set_title('ML Sky Solution - Sky*Beam/Fitted_Beam', size=8)
                if opts.inset_fov < opts.fov:
                    diff_data = a[0].real - (Sky_vec[0]*beam_grid[0]/fit_beam_grid[0])[inset_fov_inds]
            if opts.log_scale:
                im = rms_ax.imshow(np.log10(diff_data).reshape([ls_inset_grid.size]*2),
                                                origin = 'lower',
                                                extent = extent)
            else:
                im = rms_ax.imshow((diff_data).reshape([ls_inset_grid.size]*2),
                                                origin = 'lower',
                                                extent = extent)
        elif plot_ind == 3:
            if not opts.fit_beam:
                rms_ax.set_title('ML Sky Solution / Sky')
                if opts.inset_fov < opts.fov:
                    ratio_data = a[0].real/Sky_vec[0, inset_fov_inds]
                else:
                    ratio_data = a[0].real/Sky_vec[0]
            else:
                rms_ax.set_title('ML Sky Solution / Sky*Beam/Fitted_Beam', size=8)
                if opts.inset_fov < opts.fov:
                    ratio_data = a[0].real/((Sky_vec[0]*beam_grid[0]/fit_beam_grid[0])[inset_fov_inds])
            if opts.log_scale:
                im = rms_ax.imshow(np.log10(ratio_data).reshape([ls_inset_grid.size]*2),
                                                origin = 'lower',
                                                extent = extent)
            else:
                im = rms_ax.imshow((ratio_data).reshape([ls_inset_grid.size]*2),
                                                origin = 'lower',
                                                extent = extent)


    rms_ax.set_xlabel('l [deg]')
    rms_ax.set_ylabel('m [deg]')
    rms_ax.set_aspect('equal')

    imgs.append(im)

# Set title based on params
title = r'FOV: %d$\,^o$' %np.round(np.rad2deg(FOV))
if inset_fov < FOV:
    title += r', Fit FOV: %d$\,^o$' %np.round(np.rad2deg(inset_fov))
if opts.data is '':
    if opts.uniform_sky:
        title = r'Uniform Sky, ' + title
    elif opts.zenith_source:
        title = r'Zenith Source, ' + title
    elif opts.horizon_source:
        title = r'Horizon Source, ' + title
    else:
        title = r'Nsources: %d, ' %nsources + title
else:
    if 'uniform' in opts.data:
        title = r'Uniform Sky, ' + title
    elif 'zenith' in opts.data:
        title = r'Zenith Source, ' + title
    elif 'horizon' in opts.data:
        title = r'Horizon Source, ' + title
    else:
        title = r'Nsources: %d, ' %nsources + title

fig.suptitle(title)
gs.tight_layout(fig)
# gs.update(top=0.8, wspace=0.25)
gs.update(left=0.05, right=0.95, wspace=0.4, top=0.9)

if opts.force_lim:
    # Append master colorbar
    gs.update(right=0.95)
    ax_divider = make_axes_locatable(rms_ax)
    cbar_ax = ax_divider.append_axes("right", size="5%", pad="2%")
    cb = fig.colorbar(imgs[0], cax=cbar_ax)
    cb.set_clim(vmin, vmax)

    if opts.log_scale:
        cb.set_label('Log Visibilities', size=fontsize, labelpad=10)
    else:
        cb.set_label('Visibilities', size=fontsize, labelpad=10)

    cb.ax.tick_params(labelsize=fontsize)
else:
    for i,ax in enumerate(fig.axes):
        ax_divider = make_axes_locatable(ax)
        cax = ax_divider.append_axes("right", size="5%", pad="2%")
        cb = fig.colorbar(imgs[i], cax=cax)
        if i == (nplots - 1):
            if opts.log_scale:
                cb.set_label('log10', size=fontsize)
            # else:
            #     cb.set_label('Sky Amplitude', size=fontsize)

show()
