import numpy as np
import optparse, sys, os

from pyuvdata import UVData

## ----------------- Option Parser ----------------- ##
o = optparse.OptionParser()

o.add_option('--nfreqs',
    type = int,
    default = 36,
    help='Frequency(ies) of observation in MHz.')

o.add_option('--ntimes',
    type = int,
    default = 10,
    help='Frequency(ies) of observation in MHz.')

o.add_option('--data',
    type = str,
    default = '/Users/jburba/hera_things/data/golden_day/zen.2458042.12552.xx.HH.uvOR/',
    help = 'Filename for input HERA data file.')

o.add_option('--phase',
    action = 'store_true',
    help = 'If passed, phase data to first visibility lst and telescope lat.')

o.add_option('--rms',
    type = 'float',
    default=1.e-5,
    help = 'RMS to use for Gaussian noise in dataset.  Actual rms used is opts.rms/np.sqrt(2).')

o.add_option('--filepath',
    type = 'str',
    help = 'Filepath for saving output .uvfits file.')

o.add_option('--dic',
    action = 'store_true',
    help = 'If passed, output data to numpy dictionary format.')

opts,args = o.parse_args(sys.argv[1:])
print o.values

if not os.path.exists(opts.filepath):
    print opts.filepath + ' is not a valid path'
    sys.exit()

# Read in data
uvd = UVData()
uvd.read_miriad(opts.data)

if opts.phase:
    print '\nPhasing data...'
    from astropy.time import Time
    uvd.phase_to_time(Time(uvd.time_array[0], format='jd'))
    # uvd.phase(uvd.lst_array[0], uvd.telescope_location_lat_lon_alt[1], use_ant_pos=True)

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

# Generate noise data
rms = opts.rms/np.sqrt(2)
# uvd.data_array = (np.random.normal(0, rms, uvd.data_array.shape)
#                            +
#                            np.random.normal(0, rms, uvd.data_array.shape)*1j)

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
if opts.nfreqs:
    filename += '_%dnfreqs' %opts.nfreqs
if opts.ntimes:
    filename += '_%dntimes' %opts.ntimes
if opts.phase:
    filename += '_phased'
if not opts.dic:
    filename += '.uvfits'
else:
    filename += '.npy'
print 'Writing %s' %filename
if not opts.dic:
    uvd.write_uvfits(filename, spoof_nonessential=True)
else:
    out_dic = {}
    out_dic['data_array'] = uvd.data_array.copy()
    out_dic['uvw_array'] = uvd.uvw_array.copy()
    out_dic['freq_array'] = uvd.freq_array.copy()
    np.save(filename, out_dic)
