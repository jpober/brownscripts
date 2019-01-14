import numpy as np
import optparse, sys, os, copy

from pyuvdata import UVData

## ----------------- Option Parser ----------------- ##
o = optparse.OptionParser()

o.add_option('--nfreqs',
    type = int,
    default = 38,
    help='Frequency(ies) of observation in MHz.')

o.add_option('--ntimes',
    type = int,
    default = 60,
    help='Frequency(ies) of observation in MHz.')

o.add_option('--data',
    type = str,
    default = '/Users/jburba/hera_things/data/golden_day/zen.2458042.12552.xx.HH.uvOR/',
    help = 'Filename for input HERA data file.')

o.add_option('--spoof_int_time',
    action = 'store_true',
    help = 'Spoof integration time to 30 seconds for phasing of data. '
           'Defualt for HERA is 10 seconds.')

o.add_option('--phase',
    action = 'store_true',
    help = 'If passed, phase data to first visibility lst and telescope lat.')

o.add_option('--rms',
    type = 'float',
    default=1.e-5,
    help = 'RMS to use for Gaussian noise in dataset.  Actual rms used is opts.rms/np.sqrt(2).')

o.add_option('--filepath',
    type = 'str',
    help = 'Filepath for saving output file.')

o.add_option('--dic',
    action = 'store_true',
    help = 'If passed, output data to numpy dictionary format.')

o.add_option('--random',
    action = 'store_true',
    help = 'If passed, use a random UV sampling with bounds set by --data.')

o.add_option('--grid',
    action = 'store_true',
    help = 'If passed, use a gridded data set based on BayesEoR/params parameters.')

o.add_option('--beam',
    type = 'str',
    help = 'Filepath for beam file in .npy format or \'gaussian\' to produce data '
           +
           'with a gaussian beam and a fwhm given by --fwhm.')

o.add_option('--fwhm',
    type = float,
    help = 'Full width half maximum for the gaussian beam in degrees.')

opts,args = o.parse_args(sys.argv[1:])
print o.values

if not os.path.exists(opts.filepath):
    print opts.filepath + ' is not a valid path'
    sys.exit()

rms = opts.rms/np.sqrt(2)
# rms = opts.rms
# Power too low, so try multiplying by npix in eor cube to get right normalizaiton
# rms *= 512

if opts.grid:
    nu, nv = [9]*2
    nuv = nu*nv - 1
    ntimes = 1

    if opts.beam:
        if not opts.beam == 'gaussian':
            # Read in interpolated beam and compute 2d DFT using FOV and npix_side from BayesEoR.Params
            print 'Reading in beam from ' + opts.beam
            interp_beam = np.load(opts.beam)
            nf = interp_beam.shape[1]

            # DFT beam
            interp_beam = interp_beam.reshape((nu, nv, nf))
            # THIS ORDERING IS CONFUSING NEED TO SWITCH TO (nf, nx, ny)
        else:
            print 'Making Gaussian beam with FWHM: %.1f [deg]' %opts.fwhm
            # Generate sky coordinates
            nf = opts.nfreqs
            simulation_FoV_deg = 12.0
            nx, ny = [9]*2
            FOV = np.round(np.deg2rad(simulation_FoV_deg), decimals=14)
            ls = np.linspace(-FOV/2, FOV/2, nx)
            ms = np.copy(ls)
            nlm = ls.size*ms.size
            L, M = np.meshgrid(ls, ms)
            ls_vec, ms_vec = L.flatten(), M.flatten()

            # Generate Gaussian beam normalized to 1 in each frequency channel
            fwhm_deg = opts.fwhm
            fwhm = np.deg2rad(fwhm_deg)
            stddev = fwhm/(2*np.sqrt(2*np.log(2)))
            gauss_beam = np.exp(-(ls_vec**2 + ms_vec**2)/(2*stddev**2))
            gauss_beam = gauss_beam.reshape((nx, ny))
            gauss_beam = np.tile(gauss_beam, (nf, 1, 1))
            interp_beam = gauss_beam.copy()

        # Set up white noise image cube
        noise_cube = np.random.normal(0, rms, interp_beam.shape)

        # FT beam*noise cube to get visibilities
        axes_tuple = (1, 2)
        s_before_ZM = np.fft.ifftshift(noise_cube*interp_beam, axes=axes_tuple)
        s_before_ZM = np.fft.fftn(s_before_ZM, axes=axes_tuple)
        s_before_ZM = np.fft.fftshift(s_before_ZM, axes=axes_tuple)
        s_before_ZM = s_before_ZM.reshape((nf, nx*ny))

        # Remove (u, v) = (0, 0) pixel
        half_ind = nuv/2
        s_before_ZM /= s_before_ZM[:, 0].size**0.5
        data_array = np.delete(s_before_ZM, half_ind, axis=1)
        # Convert to Fortran ordering to match pyuvdata formatting of files with (nblts, nfreqs)
        data_array = data_array.flatten().reshape((nx*ny-1, nf), order='F')
    else:
        # Generate white noise visibilties
        nf = opts.nfreqs
        data_array = np.zeros((nuv*ntimes, nf), dtype='complex')
        half_ind = nuv/2
        for ntime in range(ntimes):
            for i in range(half_ind):
                noise = np.random.normal(0, rms, nf) + 1j*np.random.normal(0, rms, nf)
                data_array[ntime*nuv + i] = noise
                data_array[ntime*nuv + (nuv - i - 1)] = noise.conjugate()

    # Write output file
    if not opts.filepath[-1] == '/':
        opts.filepath += '/'
    filename = opts.filepath + 'noise_gridded_%.2erms' %rms
    filename += '_%dnfreqs' %nf
    filename += '_nu_%d' %nu
    filename += '_nv_%d' %nv
    if opts.beam:
        if not opts.beam == 'gaussian':
            filename += '_w-cst-beam'
        else:
            filename += '_w-gauss-beam_%.1fdfwhm' %fwhm_deg

    file_counter = 0
    file_path_exists = os.path.exists(filename + '.npy')
    while file_path_exists:
        file_counter += 1
        testfile = copy.copy(filename)
        testfile += '_%d' %file_counter
        file_path_exists = os.path.exists(testfile + '.npy')
    if not file_counter == 0:
        filename += '_%d' %file_counter
    filename += '.npy'

    out_dic = {}
    out_dic['data_array'] = data_array

    print 'Writing %s' %filename
    np.save(filename, out_dic)

else:
    # Read in data
    uvd = UVData()
    uvd.read_miriad(opts.data)

    # Filter data by antenna positions
    antpos, ants = uvd.get_ENU_antpos()
    pos_moduli = np.sqrt(antpos[:, 0]**2 + antpos[:, 1]**2)
    hex7_inds = np.where(pos_moduli <= 15)[0]
    antenna_nums = ants[hex7_inds]

    # Get frequencies
    freq_inds = np.where(np.logical_and(uvd.freq_array[0] <= 161e6, uvd.freq_array[0] >= 150e6))[0]
    frequencies = uvd.freq_array[0, freq_inds][:opts.nfreqs]

    # Get times
    times = np.unique(uvd.time_array)[:opts.ntimes]

    # Filter data by ants, freqs, and times
    print 'Applying select...'
    uvd.select(ant_str = 'cross')
    uvd.select(antenna_nums=antenna_nums, frequencies=frequencies, times=times)

    # Phasing
    if opts.phase:
        print '\nPhasing data...'
        from astropy.time import Time
        if opts.spoof_int_time:
            unique_times = np.unique(uvd.time_array)
            ntimes = len(unique_times)
            jd_30s_interval = unique_times[3] - unique_times[0]
            times = np.arange(unique_times[0], ntimes*jd_30s_interval + unique_times[0], jd_30s_interval)
            for ind, jd in enumerate(unique_times):
                inds = np.where(uvd.time_array == jd)[0]
                uvd.time_array[inds] = times[ind]
        uvd.phase_to_time(Time(uvd.time_array[0], format='jd'))
        # uvd.phase(uvd.lst_array[0], uvd.telescope_location_lat_lon_alt[1], use_ant_pos=True)

    # Generate noise data
    # uvd.data_array = (np.random.normal(0, rms, uvd.data_array.shape)
    #                            +
    #                            np.random.normal(0, rms, uvd.data_array.shape)*1j)
    if opts.random:
        nuv = uvd.uvw_array.shape[0]
        max_u = uvd.uvw_array[:, 0].max()
        max_v = uvd.uvw_array[:, 1].max()
        u_array = np.random.uniform(-max_u, max_u, nuv)
        v_array = np.random.uniform(-max_v, max_v, nuv)
        uvw_array = np.stack((u_array, v_array)).T
        data_array = (np.random.normal(0, rms, (nuv, opts.nfreqs))
                      +
                      1j*np.random.normal(0, rms, (nuv, opts.nfreqs)))

    else:
        uvw_array = np.copy(uvd.uvw_array)
        uvw_array = np.vstack((uvw_array, -uvw_array))
        data_array = np.zeros((uvw_array.shape[0], opts.nfreqs), dtype=complex)
        half_ind = uvw_array.shape[0]/2
        print uvw_array.shape
        for i in range(half_ind):
            print uvw_array[i, :2], uvw_array[half_ind + i, :2]
            complex_noise = np.random.normal(0, rms, (1, opts.nfreqs)) + 1j*np.random.normal(0, rms, (1, opts.nfreqs))
            data_array[i] = complex_noise
            data_array[half_ind + i] = complex_noise.conjugate()

    # lexsort_inds = np.lexsort((uvw_array[:, 0], uvw_array[:, 1]))
    # uvw_array = uvw_array[lexsort_inds]
    # data_array = data_array[lexsort_inds]

    # Save data
    filename = opts.filepath + 'noise_hera7-hex_%.1erms' %opts.rms
    if opts.phase:
        filename += '_phased'
    if opts.random:
        filename += '_random-uv'
    if opts.spoof_int_time:
        filename += '_30s-int-time'
    if opts.nfreqs:
        filename += '_%dnfreqs' %opts.nfreqs
    if opts.ntimes:
        filename += '_%dntimes' %opts.ntimes
    if not opts.dic:
        filename += '.uvfits'
    else:
        filename += '.npy'
    print 'Writing %s' %filename
    if not opts.dic:
        # This is not working rn btw
        uvd.write_uvfits(filename, spoof_nonessential=True)
    else:
        out_dic = {}
        out_dic['data_array'] = data_array.copy()
        out_dic['uvw_array'] = uvw_array.copy()
        out_dic['freq_array'] = uvd.freq_array.squeeze().copy()
        np.save(filename, out_dic)
