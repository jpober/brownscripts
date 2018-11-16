import numpy as np
import optparse, sys, os

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

opts,args = o.parse_args(sys.argv[1:])
print o.values

if not os.path.exists(opts.filepath):
    print opts.filepath + ' is not a valid path'
    sys.exit()

rms = opts.rms/np.sqrt(2)

if opts.grid:
    nu, nv = 9, 9
    nf = 38
    nuv = nu*nv - 1
    ntimes = 1
    data_array = np.zeros((nuv*ntimes, nf), dtype='complex')
    half_ind = nuv/2
    for ntime in range(ntimes):
        for i in range(half_ind):
            noise = np.random.normal(0, rms, nf) + 1j*np.random.normal(0, rms, nf)
            data_array[ntime*nuv + i] = noise
            data_array[ntime*nuv + (nuv - i - 1)] = noise.conjugate()

    filename = opts.filepath + 'noise_gridded_%.1erms' %opts.rms
    filename += '_%dnfreqs' %nf
    filename += '_%dntimes' %ntimes
    filename += '.npy'

    out_dic = {}
    out_dic['data_array'] = data_array.copy()

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
