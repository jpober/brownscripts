import numpy as np
import optparse, sys, os, copy

# Need to reinstall pyuvdata
# from pyuvdata import UVData

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
    default = 1.e-5,
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
    help = 'Filepath for beam file in .npy format, \'gaussian\', or \'uniform\' '
           +
           'to produce data with a gaussian beam and a fwhm given by --fwhm.')

o.add_option('--fwhm',
    type = float,
    help = 'Full width half maximum for the gaussian beam in degrees.')

o.add_option('--amp',
    type = float,
    default = 1.0,
    help = 'Amplitude of beam maximum if making a gaussian or uniform beam.')

o.add_option('--seed',
    action = 'store_true',
    help = 'If passed, seed np.random so datasets have the same noise cubes.')

o.add_option('--no_write',
    action = 'store_true',
    help = 'If passed, do not write the data to a file.')

o.add_option('--slope',
    type = 'float',
    default = 0.0,
    help = 'If --slope=n, apply a slope of k_z^(-n) to the delay transformed visibilities'
           +
           ' (default = 0.0).')

o.add_option('--cube_side_Mpc',
    type = 'float',
    default = 512.0,
    help = 'If passing --slope, need to specify a cube side length in Mpc to set scaling for k_z'
           +
           ' (default = 512).')

o.add_option('--cube_side_pix',
    type = 'int',
    default = 128,
    help = 'If passing --slope, need to specify a cube side length in pixels to set scaling for k_z'
           +
           ' (default = 128).')

opts,args = o.parse_args(sys.argv[1:])
# print o.values

if not (opts.no_write or os.path.exists(opts.filepath)):
    print opts.filepath + ' is not a valid path'
    sys.exit()

rms = opts.rms/np.sqrt(2)

if opts.grid:
    nu, nv = [9]*2
    nuv = nu*nv - 1
    ntimes = 1
    nf = opts.nfreqs

    if opts.beam:
        if not opts.beam == 'gaussian' and not opts.beam == 'uniform':
            # Read in interpolated beam and compute 2d DFT using FOV and npix_side from BayesEoR.Params
            print 'Reading in beam from ' + opts.beam
            interp_beam = np.load(opts.beam)
            nf = interp_beam.shape[1]

            # DFT beam
            beam_cube = interp_beam.reshape((nu, nv, nf))
            # THIS ORDERING IS CONFUSING NEED TO SWITCH TO (nf, nx, ny)
        elif opts.beam == 'gaussian':
            # print 'Making Gaussian beam with FWHM: %.1f [deg]' %opts.fwhm
            # Generate sky coordinates
            simulation_FoV_deg = 12.0
            nx, ny = (nu, nv)
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
            # gauss_beam = np.exp(-(ls_vec**2 + ms_vec**2)/(2*stddev**2))
            # gauss_beam = np.tile(gauss_beam, nf)
            # gauss_beam = gauss_beam.reshape((nf, nx, ny))
			
            # Try making beam slice and tiling to make sure the beam
            # cube is what I think it is
            gauss_beam_slice = np.exp(-(L**2 + M**2)/(2*stddev**2))
            gauss_beam = np.zeros((nf, nu, nv))
            for i_f in range(nf):
				gauss_beam[i_f] = gauss_beam_slice
			
            beam_cube = gauss_beam.copy()*opts.amp
        else:
            uniform_beam = np.ones((nf, nu, nv))*opts.amp
            beam_cube = uniform_beam.copy()


        # Set up white noise image cube
        if opts.seed:
            # print 'Setting seed for np.random'
            np.random.seed(12345)
        noise_cube = np.random.normal(0, rms, beam_cube.shape)

        # FT beam*noise cube to get visibilities
        axes_tuple = (1, 2)
        s_before_ZM = np.fft.ifftshift(noise_cube*beam_cube, axes=axes_tuple)
        s_before_ZM = np.fft.fftn(s_before_ZM, axes=axes_tuple)
        s_before_ZM = np.fft.fftshift(s_before_ZM, axes=axes_tuple)
        s_before_ZM /= s_before_ZM[0].size**0.5
		
        if not opts.slope == 0.0:
			# Need to make this symmetric around kz=0 at kz[nf/2]
			z, y, x = np.mgrid[-(nf/2):(nf/2),-(nu/2):(nu/2)+1,-(nv/2):(nv/2)+1]
			z = z.astype('float64')
			data_side_Mpc_scaled = opts.cube_side_Mpc*nf/opts.cube_side_pix
			deltakpara = 2*np.pi/data_side_Mpc_scaled
			kz = z*deltakpara

			kz_weights = 1./np.abs(kz)**opts.slope
			kz_weights[[0, 1, nf/2, -1]] = 0.0 # no power in modes that are removed
			kz_weights /= kz_weights.max()

			# Apply weighting in 3D k space before flattening
			# FFT from uvf to uvnu
			axes_tuple = (0,)
			vfft_s = np.fft.ifftshift(s_before_ZM.copy(), axes=axes_tuple)
			vfft_s = np.fft.fftn(vfft_s, axes=axes_tuple)
			vfft_s = np.fft.fftshift(vfft_s, axes=axes_tuple)
			# Aplly weights
			vfft_s *= kz_weights
			# Inverse FFT back to uvf cube
			s_before_ZM = np.fft.ifftshift(vfft_s.copy(), axes=axes_tuple)
			s_before_ZM = np.fft.ifftn(s_before_ZM, axes=axes_tuple)
			s_before_ZM = np.fft.fftshift(s_before_ZM, axes=axes_tuple)
		
		# Current method
        s_before_ZM_flat = s_before_ZM.flatten()
        ZM_mask = np.ones_like(s_before_ZM_flat).astype('bool')
        ZM_mask[(nu*nv/2)::(nu*nv)] = False
        data_array = s_before_ZM_flat[ZM_mask]
		
        # Old method
        # s_before_ZM = s_before_ZM.reshape((nf, nx*ny))

        # Remove (u, v) = (0, 0) pixel
        # half_ind = nuv/2
        # data_array = np.delete(s_before_ZM, half_ind, axis=1)
        # Convert to Fortran ordering to match pyuvdata formatting of files with (nblts, nfreqs)
        # data_array = data_array.flatten().reshape((nx*ny-1, nf), order='F')
        # data_array = data_array.flatten()

    else:
        # Set up white noise image cube
        if opts.seed:
            print 'Setting seed for np.random'
            np.random.seed(12345)
        noise_cube = np.random.normal(0, rms, (nf, nu, nv)) + 0j
        # noise_cube = np.ones((nf, nu, nv)) + 0j
        # noise_cube_face = np.arange(1, nu*nv + 1)
        # noise_cube = np.tile(noise_cube_face, (nf, 1))
        # noise_cube = noise_cube.reshape((nf, nu, nv))

        # FT noise cube to get visibilities (uvf cube)
        axes_tuple = (1, 2)
        s_before_ZM = np.fft.ifftshift(noise_cube, axes=axes_tuple)
        s_before_ZM = np.fft.fftn(s_before_ZM, axes=axes_tuple)
        s_before_ZM = np.fft.fftshift(s_before_ZM, axes=axes_tuple)
        s_before_ZM /= s_before_ZM[0].size**0.5

        if not opts.slope == 0.0:
			# Need to make this symmetric around kz=0 at kz[nf/2]
			z, y, x = np.mgrid[-(nf/2):(nf/2),-(nu/2):(nu/2)+1,-(nv/2):(nv/2)+1]
			z = z.astype('float64')
			data_side_Mpc_scaled = opts.cube_side_Mpc*nf/opts.cube_side_pix
			deltakpara = 2*np.pi/data_side_Mpc_scaled
			kz = z*deltakpara

			kz_weights = 1./np.abs(kz)**opts.slope
			kz_weights[[0, 1, nf/2, -1]] = 0.0 # no power in modes that are removed
			kz_weights /= kz_weights.max()

			# Apply weighting in 3D k space before flattening
			# FFT from uvf to uvnu
			axes_tuple = (0,)
			vfft_s = np.fft.ifftshift(s_before_ZM.copy(), axes=axes_tuple)
			vfft_s = np.fft.fftn(vfft_s, axes=axes_tuple)
			vfft_s = np.fft.fftshift(vfft_s, axes=axes_tuple)
			# Aplly weights
			vfft_s *= kz_weights
			# Inverse FFT back to uvf cube
			s_before_ZM = np.fft.ifftshift(vfft_s.copy(), axes=axes_tuple)
			s_before_ZM = np.fft.ifftn(s_before_ZM, axes=axes_tuple)
			s_before_ZM = np.fft.fftshift(s_before_ZM, axes=axes_tuple)
        
        s_before_ZM_flat = s_before_ZM.flatten()
        ZM_mask = np.ones_like(s_before_ZM_flat).astype('bool')
        ZM_mask[(nu*nv/2)::(nu*nv)] = False
        data_array = s_before_ZM_flat[ZM_mask]
    
    
    if not opts.no_write:
        # Write output file
        if not opts.filepath[-1] == '/':
            opts.filepath += '/'
        filename = opts.filepath + 'noise_gridded_%.2erms' %rms
        filename += '_%dnfreqs' %nf
        filename += '_nu_%d' %nu
        filename += '_nv_%d' %nv
        if opts.beam:
            if not opts.beam == 'gaussian' and not opts.beam == 'uniform':
                filename += '_w-cst-beam'
            elif opts.beam == 'gaussian':
                filename += '_w-gauss-beam_%.1fdfwhm' %fwhm_deg
            elif opts.beam == 'uniform':
                filename += '_w-uniform-beam'
            # append beam amplitude
            filename += '_%.1famp' %opts.amp
        if not opts.slope == 0.0:
            filename += ('_%.1elogk-slope' %opts.slope).replace('.', 'd')
        if len(data_array.shape) < 2:
            filename += '_vectorized'

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

		
		
		
		
###
# OLD CODE
###

#             kz_weights = np.ones_like(kz)
#             # Force a bunch of delays to zero
#             kz_weights[:11] = 0.0
#             kz_weights[-10:] = 0.0
#             kz_weights[nf/2] = 0.0
#             kz_weights = 1 - np.abs(kz)/np.abs(kz).max()
#             kz_weights[:2] = 0.0
#             kz_weights[-1] = 0.0
#             kz_weights = kz.copy()
#             kz_weights = 1.e3*(np.abs(kz).max() - np.abs(kz))/np.abs(kz).max()
#             # Don't want power in the highest modes (2 for kz<0, 1 for kz>0)
#             axes_tuple = (0, 1, 2)
#             noise_cube_ft = np.fft.ifftshift(noise_cube.copy(), axes=axes_tuple)
#             noise_cube_ft = np.fft.fftn(noise_cube_ft, axes=axes_tuple)
#             noise_cube_ft = np.fft.fftshift(noise_cube_ft, axes=axes_tuple)
#             # noise_cube_ft *= np.abs(opts.slope*kz)
# #             noise_cube_ft *= 1 - np.abs(kz)/np.abs(kz).max()
#             noise_cube_ft += kz_weights # addition or multiplication?
#             # Don't need to FT back just to FT again
#             noise_cube = np.fft.ifftshift(noise_cube_ft.copy(), axes=axes_tuple)
#             noise_cube = np.fft.ifftn(noise_cube, axes=axes_tuple)
#             noise_cube = np.fft.fftshift(noise_cube, axes=axes_tuple)
#         # sys.exit()
#         """

# Reshaping code
# Reshape array to mimic pyuvdata formatting
# s_before_ZM = s_before_ZM.reshape((nf, nu*nv))

# Remove (u, v) = (0, 0) pixel
# half_ind = nuv/2
# data_array = np.delete(s_before_ZM, half_ind, axis=1)
# Convert to Fortran ordering to match pyuvdata formatting of files with (nblts, nfreqs)
# data_array = data_array.flatten().reshape((nu*nv-1, nf), order='F')