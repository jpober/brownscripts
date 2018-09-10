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

o.add_option('--write',
    action = 'store_true',
    help = 'If passed, write fitted RMS and associated errors to a .npy file.')

o.add_option('--uvdata',
    type = str,
    help = 'Filename for input uvw binary .npy file.')

o.add_option('--poly_order',
    type = str,
    default = '2',
    help = 'If n, fit beam variation as a function of frequency with an n-th order polynomial.')

o.add_option('--rms_data',
    type = str,
    help = 'File path for previously generated rms data via this script.')

opts,args = o.parse_args(sys.argv[1:])


if not opts.rms_data:
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

    # Get polynomal order(s)
    if ',' in opts.poly_order:
        poly_orders = np.sort(map(int, opts.poly_order.split(',')))
    elif '-' in opts.poly_order:
        poly_array = np.sort(map(int, opts.poly_order.split('-')))
        poly_orders = np.arange(poly_array[0], poly_array[1] + 1, 1)
        poly_orders = poly_orders[np.where(poly_orders <= poly_array[1])]
    elif not (',' in opts.poly_order and '-' in opts.poly_order):
        poly_orders = [int(opts.poly_order)]

    npolys = len(poly_orders)



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

    # Build sky from opts params
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
    for i in range(nfreqs):
        # Construct visibilities faster as FT of beam*sky
        if opts.beam:
            Vs[i] = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift((beam_grid[i]*Sky_vec[i]).reshape([npix_side]*2)))).flatten()
        else:
            Vs[i] = np.fft.fftshift(np.fft.fft2(np.fft.ifftshift((Sky_vec[i]).reshape([npix_side]*2)))).flatten()



    ## ---------------------------------- Construct MaxL Sky ---------------------------------- ##
    print 'Constructing maximum likelihood sky...'
    # Construct solution using analytic solution for maximum likelihood
    # Assumes a Gaussian log likelihood function
    # Requires noise injection into data (visibilities above)
    a = np.zeros_like(Vs)
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

    # Create data containers for storing polynomial fit data
    rms_data = np.zeros([npolys, nlm])
    poly_fits = np.zeros([npolys, nfreqs, nlm])
    # rms_data = np.zeros_like(a[0])
    # poly_fits = np.zeros_like(a)

    print 'Performing polynomial fit...'
    print 'Polynomial order: ',
    for poly_ind, poly_order in enumerate(poly_orders):
        if not poly_order == poly_orders[-1]:
            print '%d, ' %poly_order ,
        else:
            print '%d' %poly_order
        # Fit quadratics in frequency along frequency axis if a
        fit_coeffs = np.polyfit(freqs, np.abs(a), poly_order)
        # shape (poly_order + 1, npix), i.e. (poly_order + 1) coefficients for each pixel

        for pix_ind in range(npix):
            for j in range(poly_order + 1):
                poly_fits[poly_ind, :, pix_ind] += fit_coeffs[j, pix_ind]*freqs**(poly_order - j)

        residuals = np.abs(a) - np.abs(poly_fits[poly_ind])
        rms_data[poly_ind] = np.std(residuals, axis=0)

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

        filename += '_%dnpix-side' %opts.npix_side
        filename += '_%dpoly' %opts.poly_order

        if opts.zenith_source:
            filename += '_zenith-source'
        elif opts.horizon_source:
            filename += '_horizon-source'
        elif opts.uniform_sky:
            filename += '_uniform-sky'
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

        print 'Writing ' + filename + '.npy ...\n'
        out_dic = {}
        out_dic['sky'] = Sky
        out_dic['vis'] = Vs
        out_dic['maxL_sky'] = a
        out_dic['input_rms'] = opts.rms
        out_dic['freqs'] = freqs
        if opts.beam:
            out_dic['beam_file'] = opts.beam
        if opts.fit_beam:
            out_dic['fitted_beam'] = True
        else:
            out_dic['fitted_beam'] = False
        out_dic['rms_data'] = rms_data
        out_dic['poly_order'] = opts.poly_order
        np.save(filename + '.npy', out_dic)

        sys.exit()

else:
    print 'Reading in %s...' %opts.rms_data
    # Read in rms data from opts.rms_data
    data_dic = np.load(opts.rms_data).item()
    rms_data = data_dic['rms_data']
    poly_order_str = data_dic['poly_order']
    if '-' in poly_order_str:
        poly_order_arr = map(int, poly_order_str.split('-'))
    elif ',' in poly_order_str:
        poly_order_arr = map(int, poly_order_str.split(','))
    else:
        poly_order_arr = int(poly_order_str)
    if len(poly_order_arr) > 1:
        poly_orders = np.arange(poly_order_arr[0], poly_order_arr[-1] + 1, 1)
    else:
        poly_orders = poly_order_arr
    npix_ind = opts.rms_data.find('npix')
    back_ind = opts.rms_data.find('_', npix_ind - 4)
    npix_side = int(opts.rms_data[back_ind+1:npix_ind])
    dfov_ind = opts.rms_data.find('dfov')
    back_ind = opts.rms_data.find('_', dfov_ind - 4)
    FOV = np.deg2rad(float(opts.rms_data[back_ind + 1:dfov_ind]))
    ls = np.linspace(-FOV/2, FOV/2, npix_side)
    ms = np.copy(ls)
    nlm = ls.size*ms.size
    lm_pixel_half = np.diff(ls)[0]/2.
    if not any(x in opts.rms_data for x in ['uniform', 'zenith', 'horizon']):
        nsources = np.sum(data_dic['sky'])


## ---------------------------------- Plotting ---------------------------------- ##
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib.pyplot import *

print 'Polynomial orders: ',
print map(str, poly_orders)
print 'Plotting...'

nplots = len(poly_orders)
plot_size = 5
nrows = 1
ncols = nplots
fig = figure(figsize=(plot_size*nplots, plot_size))

fontsize = 14
extent_lm = [ls.min() - lm_pixel_half, ls.max() + lm_pixel_half,
                    ms.min() - lm_pixel_half, ms.max() + lm_pixel_half]
extent_lm = np.rad2deg(extent_lm)

gs = gridspec.GridSpec(nrows, nplots)
imgs = []
if nplots == 1:
    rms_ax = fig.add_subplot(gs[0])
    if opts.log_scale:
        rms_im = rms_ax.imshow(10*np.log10(rms_data).reshape([npix_side]*2))
    else:
        rms_im = rms_ax.imshow(rms_data.reshape([npix_side]*2))
    rms_ax.set_xlabel('l')
    rms_ax.set_ylabel('m')
    rms_ax.set_title('Polynomial Order: %d' %poly_orders[0])

    imgs .append(rms_im)

else:
    if opts.force_lim:
        # Add master labels for all x and y axes
        master_ax = fig.add_subplot(gs[:,:])
        master_ax.set_xticks([])
        master_ax.set_yticks([])
        master_ax.spines['top'].set_visible(False)
        master_ax.spines['right'].set_visible(False)
        master_ax.spines['bottom'].set_visible(False)
        master_ax.spines['left'].set_visible(False)

    for plot_ind in range(nplots):
        rms_ax = fig.add_subplot(gs[plot_ind])

        if plot_ind == 0:
            if opts.force_lim:
                if opts.log_scale:
                    vmin = np.log10(rms_data).min()
                    vmax = np.log10(rms_data).max()
                else:
                    vmin = rms_data.min()
                    vmax = rms_data.max()

        if opts.force_lim:
            if opts.log_scale:
                rms_im = rms_ax.imshow(np.log10(rms_data)[plot_ind].reshape([npix_side]*2),
                                                       origin='lower',
                                                       extent=extent_lm,
                                                       vmin=vmin,
                                                       vmax=vmax)
            else:
                rms_im = rms_ax.imshow(rms_data[plot_ind].reshape([npix_side]*2),
                                                       origin='lower',
                                                       extent=extent_lm,
                                                       vmin=vmin,
                                                       vmax=vmax)

            # remove axis labels from interior plots
            if not plot_ind == 0:
                rms_ax.set_yticks([])
            else:
                rms_ax.set_ylabel('m')

        else:
            if opts.log_scale:
                rms_im = rms_ax.imshow(np.log10(rms_data)[plot_ind].reshape([npix_side]*2),
                                                       origin='lower',
                                                       extent=extent_lm)
            else:
                rms_im = rms_ax.imshow(rms_data[plot_ind].reshape([npix_side]*2),
                                                       origin='lower',
                                                       extent=extent_lm)

        rms_ax.set_title('Polynomial Order: %d' %poly_orders[plot_ind])
        rms_ax.set_xlabel('l')

        imgs.append(rms_im)

# Set title based on params
title = r'FOV: %d$\,^o$' %np.round(np.rad2deg(FOV))
if opts.rms_data is None:
    if opts.uniform_sky:
        title = r'Uniform Sky, ' + title
    elif opts.zenith_source:
        title = r'Zenith Source, ' + title
    elif opts.horizon_source:
        title = r'Horizon Source, ' + title
    else:
        title = r'Nsources: %d, ' %nsources + title
else:
    if 'uniform' in opts.rms_data:
        title = r'Uniform Sky, ' + title
    elif 'zenith' in opts.rms_data:
        title = r'Zenith Source, ' + title
    elif 'horizon' in opts.rms_data:
        title = r'Horizon Source, ' + title
    else:
        title = r'Nsources: %d, ' %nsources + title

if opts.fractional_fit:
    title = r'Fractional fit, ' + title
else:
    title = r'Absolute fit, ' + title

fig.suptitle(title)
gs.tight_layout(fig)
gs.update(top=0.8)

if opts.force_lim:
    # Append master colorbar
    gs.update(right=0.85)
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
                cb.set_label('log10(Fitted RMS)', size=fontsize)
            else:
                cb.set_label('Fitted RMS', size=fontsize)

show()
