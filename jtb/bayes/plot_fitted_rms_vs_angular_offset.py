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
    help = 'Channel width in MHz if --freq passed with \'-\'.')

o.add_option('--rms',
    type = float,
    default = 1.e-5,
    help = 'RMS for noise injection.')

o.add_option('--l_walk',
    action = 'store_true',
    default = True,
    help = 'If passed, walk along pixels with constant m.')

o.add_option('--m_walk',
    action = 'store_true',
    help = 'If passed, walk along pixels with consant l.')

o.add_option('--l_offset',
    type = float,
    help = 'Moves source in l-direction by l_offset*ls.max().  Must be between 0 and 1.')

o.add_option('--m_offset',
    type = float,
    help = 'Moves source in m-direction by m_offset*ms.max().  Must be between 0 and 1.')

o.add_option('--write',
    action = 'store_true',
    help = 'If passed, write fitted RMS and associated errors to a .npy file.')

o.add_option('--rms_data',
    type = str,
    help = 'Path to file with generated fitted RMS data.')

o.add_option('--beam',
    type = str,
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

o.add_option('--npix_centers',
    type = str,
    default = '10',
    help = ('Number of pixel centers to use as source is moved from zenith to horizon. '
              +
              'Cannot exceed opts.npix_side/2 + 1 (i.e. <= 16 for opts.npix_side = 31).'))

o.add_option('--noplot',
    action = 'store_true',
    help = 'If passed, suppress plotting.')

opts,args = o.parse_args(sys.argv[1:])

# If no data passed, create fitted RMS data
if not opts.rms_data:

    if opts.m_walk:
        opts.l_walk = False


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
    npix_side = opts.npix_side
    npix = npix_side**2

    # Noise injection params
    N_inv = np.eye(npix)/opts.rms**2


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
    if opts.l_offset or opts.m_offset:
        if len(ls) % 2 == 0:
            mid_l = len(ls)/2
        else:
            mid_l = int(len(ls)/2. - 0.5)
        if len(ms) % 2 == 0:
            mid_m = len(ms)/2
        else:
            mid_m = int(len(ms)/2. - 0.5)

        if opts.l_offset:
            l_off = int(mid_l*opts.l_offset)
        else:
            l_off = 0

        if opts.m_offset:
            m_off = int(mid_m*opts.m_offset)
        else:
            m_off = 0

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

    # Construct u,v grid
    us_grid = np.fft.fftshift(np.fft.fftfreq(ls.shape[0], d=np.mean(np.diff(ls))))
    vs_grid = np.fft.fftshift(np.fft.fftfreq(ms.shape[0], d=np.mean(np.diff(ms))))
    us_vec, vs_vec = np.zeros(0), np.zeros(0)
    for v in vs_grid:
        for u in us_grid:
            us_vec = np.append(us_vec, u)
            vs_vec = np.append(vs_vec, v)

    # Set up angular offsets, arrays to store fitted RMS values
    ls_pos = np.copy(ls[ls >= 0])
    if opts.npix_centers.isdigit() and int(opts.npix_centers) > len(ls_pos):
        print 'npix_centers cannot exceed npix_side/2 + 1'
        sys.exit()

    if len(ls_pos) < 10:
        npix_centers = len(ls_pos)
    elif opts.npix_centers.isdigit():
        npix_centers = int(opts.npix_centers)
    else:
        npix_centers = len(ls_pos)

    angular_offsets = np.zeros(npix_centers)
    angular_offsets[0] = ls_pos[0]
    angular_offsets[-1] = ls_pos[-1]
    angular_offsets[1:-1] = np.sort(np.random.choice(ls_pos[1:-1], npix_centers - 2, replace=False))

    sys.exit()

    fitted_RMS = np.zeros((nfreqs, angular_offsets.size))
    fitted_RMS_err = np.zeros_like(fitted_RMS)

    # Store source positions for reference
    true_pos = np.zeros((angular_offsets.size, 2))

    if opts.write:
        # Write fitted RMS data
        if os.path.exists('./sim_vis/'):
            if nfreqs > 1:
                filename = 'sim_vis/maxL_fitted_RMS_%sMHz_%sMHz_%.0fdfov' %(opts.freq,
                                                                                                                    opts.freq_res,
                                                                                                                    np.rad2deg(FOV))
            else:
                filename = 'sim_vis/maxL_fitted_RMS_%sMHz_%.0fdfov' %(opts.freq,
                                                                                                        np.rad2deg(FOV))
        else:
            if nfreqs > 1:
                filename = 'maxL_fitted_RMS_%sMHz_%sMHz_%.0fdfov' %(opts.freq,
                                                                                                        opts.freq_res,
                                                                                                        np.rad2deg(FOV))
            else:
                filename = 'maxL_fitted_RMS_%sMHz_%.0fdfov' %(opts.freq,
                                                                                            np.rad2deg(FOV))

        print 'Writing ' + filename + '.npy ...\n'
        out_dic = {}
        out_dic['input_rms'] = opts.rms
        out_dic['angular_offsets'] = angular_offsets
        if opts.beam:
            out_dic['beam'] = opts.beam
        if opts.fit_beam:
            filename += '_fitbeam'

    # Precompute DFT matrix
    DFT = np.exp(-1j*2*np.pi*(np.outer(us_vec, ls_vec) + np.outer(vs_vec, ms_vec)))

    # Create arrays to store visibilities and maxL solutions
    Vs = np.zeros((nfreqs, npix), dtype=complex)
    d = np.copy(Vs) # copy of visibilties for noise injection
    a = np.zeros_like(Vs)
    Vs_maxL = np.zeros_like(a)

    print 'Frequency [MHz]: ',
    for freq_ind in range(nfreqs):
        if freq_ind == nfreqs - 1:
            print '%.0f' %freqs[freq_ind]
        else:
            print '%.0f, ' %freqs[freq_ind] ,

        # Compute components of maxL solution
        if opts.beam:
            if opts.fit_beam:
                P = np.diag(fit_beam_grid[freq_ind])
            else:
                P = np.diag(beam_grid[freq_ind])
            DFTP = np.dot(DFT, P)
            inv_part = np.linalg.inv(np.dot(np.dot(DFTP.conj().T, N_inv), DFTP))
        else:
            inv_part = np.linalg.inv(np.dot(np.dot(DFT.conj().T, N_inv), DFT))

        # Loop over sky positions
        for offset_ind, angular_offset in enumerate(angular_offsets):
            if opts.l_walk:
                if opts.m_offset:
                    true_pos[offset_ind, 1] = ms[mid_m + m_off]
                true_pos[offset_ind, 0] = angular_offset
            elif opts.m_walk:
                if opts.l_offset:
                    true_pos[offset_ind, 0] = ls[mid_l + l_off]
                true_pos[offset_ind, 1] = angular_offset

            # Get location of source in flattened image array
            pos_vec = np.where(np.logical_and(ls_vec == true_pos[offset_ind, 0],
                                                                ms_vec == true_pos[offset_ind, 1]))[0][0]

            for j in range(us_vec.size):
                if opts.beam:
                    Vs[freq_ind, j] = np.sum(beam_grid[freq_ind, pos_vec]*Vs_func(us_vec[j],
                                                                                                        true_pos[offset_ind, 0],
                                                                                                        vs_vec[j],
                                                                                                        true_pos[offset_ind, 1]))
                else:
                    Vs[freq_ind, j] = np.sum(Vs_func(us_vec[j],
                                                                      true_pos[offset_ind, 0],
                                                                      vs_vec[j],
                                                                      true_pos[offset_ind, 1]))

            # Vs = Vs.reshape((nfreqs, npix_side, npix_side))
            d[freq_ind] = np.copy(Vs[freq_ind])

            half_ind = int(us_vec.size/2.) + 1
            for j,[u,v] in enumerate(np.stack((us_vec, vs_vec), axis=1)[:half_ind]):
                neg_ind = np.where(np.logical_and(us_vec == -u, vs_vec == -v))[0][0]
                d[freq_ind, [j, neg_ind]] += (np.random.normal(0, opts.rms, 1)
                                                          +
                                                          1j*np.random.normal(0, opts.rms, 1))

            if opts.beam:
                right_part = np.dot(np.dot(DFTP.conj().T, N_inv), d[freq_ind])
            else:
                right_part = np.dot(np.dot(DFT.conj().T, N_inv), d[freq_ind])

            # Maximum likelihood solution for the sky
            a[freq_ind] = np.dot(inv_part, right_part)

            # Generate visibilities from maximum liklihood solution
            if opts.beam:
                Vs_maxL[freq_ind] = np.dot(DFTP,  a[freq_ind])
            else:
                Vs_maxL[freq_ind] = np.dot(DFT, a[freq_ind])

            # Compute fitted RMS
            diff_data = np.abs(Vs[freq_ind]) - np.abs(Vs_maxL[freq_ind])
            counts, bins = np.histogram(diff_data, bins=50)
            bin_width = np.mean(np.diff(bins))
            fit_xs = bins[:-1] + bin_width/2
            guess_params = [np.max(counts), 0.0, opts.rms]
            # fit_params: 0, amplitude; 1, mean; 2, std dev
            fit_params, fit_cov = curve_fit(Gaussian, fit_xs, counts, p0=guess_params)

            # Store fitted RMS and fit error
            fitted_RMS[freq_ind, offset_ind] = fit_params[2]
            fitted_RMS_err[freq_ind, offset_ind] = fit_cov[-1, -1]

            if opts.write:
                out_dic['fit_rms'] = fitted_RMS
                out_dic['fit_rms_err'] = fitted_RMS_err
                np.save(filename + '.npy', out_dic)


# If data already exists, load it in to plot it
else:
    data_dic = np.load(opts.rms_data).item()
    fitted_RMS = data_dic['fit_rms']
    fitted_RMS_err = data_dic['fit_rms_err']
    angular_offsets = data_dic['angular_offsets']
    opts.rms = data_dic['input_rms']
    nfreqs = fitted_RMS.shape[0]

    freq_range = map(float, opts.rms_data.split('_')[4].strip('MHz').split('-'))
    if len(freq_range) > 1:
        freq_res = float(opts.rms_data.split('_')[5].strip('MHz'))
        freqs = np.arange(freq_range[0], freq_range[1] + freq_res, freq_res)
    else:
        freqs = freq_range[0]

if not opts.noplot:
    # Plotting
    from matplotlib.pyplot import *
    from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

    fig = figure(figsize=(16,8))
    gs = gridspec.GridSpec(1, 1)
    # gs = gridspec.GridSpec(1, 2)

    # path_ax = fig.add_subplot(gs[0, 0])
    # path_ax.scatter(true_pos[:, 0], true_pos[:, 1])
    # path_ax.set_xlim([ls.min(), ls.max()])
    # path_ax.set_ylim([ms.min(), ms.max()])
    # path_ax.set_title('Walk Pattern', size=16)

    # rms_ax = fig.add_subplot(gs[0, 1])
    rms_ax = fig.add_subplot(gs[0])
    plot_xs = np.rad2deg(np.copy(angular_offsets))
    colors = np.copy(cm.jet(np.linspace(0, 1, nfreqs)))
    for freq_ind in range(nfreqs):
        im = rms_ax.errorbar(plot_xs, fitted_RMS[freq_ind]/opts.rms,
                                        yerr=fitted_RMS_err[freq_ind]/opts.rms,
                                        label='%.0f MHz' %freqs[freq_ind],
                                        color=colors[freq_ind])
    rms_ax.set_xlabel('Source Location [deg]', size=16)
    rms_ax.set_ylabel('Fitted RMS/$10^{%.0f}$' %np.log10(opts.rms), size=16)
    rms_ax.set_ylim([0.5, 1.5])
    rms_ax.set_xlim([plot_xs.min(), plot_xs.max()])
    # rms_ax.set_title(' '.join(sys.argv))
    # rms_ax.set_title(sys.argv)

    for ax in fig.axes:
        ax.tick_params(which='both', labelsize=16)
        ax.legend(loc='best')

    gs.tight_layout(fig)
    show()
