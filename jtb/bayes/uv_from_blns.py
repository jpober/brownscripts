import numpy as np
import optparse, sys, os
import matplotlib.gridspec as gridspec
import pyuvdata.utils as uvu

from pyuvdata import UVData
from pyuvdata.data import DATA_PATH
from matplotlib.pyplot import *

# Option parser
o = optparse.OptionParser()
o.add_option('--data',
    type=str,
    help='Filename for input data.')
o.add_option('--positions',
    type=str,
    help='Filename with antenna positions.')
o.add_option('--freq',
    type=str,
    help='Frequency(ies) of observation in MHz.')
o.add_option('--freq_res',
    type = float,
    default = 1.0,
    help = 'Channel width in MHz if --freq passed with \'-\'.')
o.add_option('--width',
    type=float,
    help='Width in degrees of observational period.')
o.add_option('--snapshot',
    action='store_true',
    help='If passed, only a snapshot of (u,v,w) at zenith is generated.')
o.add_option('--write',
    action='store_true',
    help='If passed, write uvws array as .npy file.')
o.add_option('--plot',
    action='store_true',
    help = 'If passed, a plot of the generated (u,v) coordinates is generated.')
opts,args = o.parse_args(sys.argv[1:])

if opts.data is None and opts.positions is None:
    opts.data = os.path.join(DATA_PATH,'hera_testfile/')

# Setup
c = 3.e8

# Get frequencies to generate uv sampling for
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

if opts.snapshot:
    hours = np.array([np.pi/2])
else:
    low_lim = np.pi/2-0.5*opts.width*np.pi/180
    high_lim = np.pi/2+0.5*opts.width*np.pi/180
    hours = np.linspace(low_lim, high_lim, 10)
decs = np.ones_like(hours)*np.pi/2 #zenith transitting point source

uv = UVData()
if opts.data is not None:
    # Read in telescope data
    if opts.data.endswith('.uvfits'):
        uv.read_uvfits(opts.data)
    else:
        uv.read_miriad(opts.data)
    ant_pos = uv.antenna_positions + uv.telescope_location
    lat,lon,alt = uv.telescope_location_lat_lon_alt
    ant_pos = uvu.ENU_from_ECEF(ant_pos.T, lat, lon, alt).T
    blns = uv.get_baseline_nums()

elif opts.positions is not None:
    # Read in antenna positions from file
    ant_pos = np.load(opts.positions)
    blns = np.zeros(0, dtype=int)
    # ant_1_array = np.zeros(0, dtype=int)
    # ant_2_array = np.zeros(0, dtype=int)
    num_ants = ant_pos.shape[0]
    ant_nums = range(1, num_ants+1)
    # for i in range(len(ant_nums)):
    for ant1 in ant_nums:
        # ant1 = ant_nums[i]
        # for j in range(len(ant_nums)):
        for ant2 in ant_nums:
            # ant2 = ant_nums[j]
            if ant1 == ant2:
                continue
            # else:
            #     ant_1_array = np.append(ant_1_array, ant1)
            #     ant_2_array = np.append(ant_2_array, ant2)
            if (uv.antnums_to_baseline(ant1, ant2) not in blns
                and uv.antnums_to_baseline(ant2, ant1) not in blns):
                blns = np.append(blns, uv.antnums_to_baseline(ant1, ant2))

# New line, better output formatting
print ''
if len(freqs) == 1:
    print 'Frequency (MHz): %f' %freqs[0]
else:
    freqs_str = ','.join(map(str, np.round(freqs, decimals=2)))
    print 'Frequency range (MHz): %s' %','.join(map(str, [freqs.min(), freqs.max()]))

# uvws = np.zeros((3, hours.shape[0], blns.shape[0], nfreqs))
uvws = np.zeros((nfreqs, hours.shape[0], blns.shape[0], 3))

for j,freq in enumerate(freqs):
    # Get baseline
    l = c/(freq*1e6) #m, wavelength
    for i, bln in enumerate(blns):
        antpair = uv.baseline_to_antnums(bln)
        # make all baselines point in same direction
        bln_vec = (ant_pos[antpair[1]-1,:] - ant_pos[antpair[0]-1,:])/l
        uvws[j,:,i,0] = np.sin(hours)*bln_vec[0] + np.cos(hours)*bln_vec[1]
        uvws[j,:,i,1] = (-np.sin(decs)*(np.cos(hours)*bln_vec[0] - np.sin(hours)*bln_vec[1])
                                + np.cos(decs)*bln_vec[2])
        uvws[j,:,i,2] = (np.cos(decs)*(np.cos(hours)*bln_vec[0] - np.sin(hours)*bln_vec[2])
                                + np.sin(decs)*bln_vec[2])


# Keep only the unique (u,v,w)
uvws *= -1.
uvws = np.unique(np.round(uvws, decimals=10), axis=2)

if opts.snapshot:
    # Only take UV points with angles between -pi/4 and 3*pi/4 (symmetry line for visibilities)
    uvws = np.append(-uvws, uvws, axis=2)
    uvws_meters = uvws[0, 0, :]*3.e8/(freqs[0]*1.e6)
    angles = np.angle(uvws_meters[:, 0] + 1j*uvws_meters[:, 1])
    sym_inds = np.where(np.logical_and(angles >= -np.pi/4, angles <= 3*np.pi/4))[0]
    # sym_inds = np.where(np.logical_and(angles >= -np.pi/2, angles <= np.pi/2))[0]
    # sym_inds = np.where(np.logical_and(angles >= 0, angles <= np.pi))[0]
    uvws = uvws[:, :, sym_inds]
    uvws_meters = uvws_meters[sym_inds]

if opts.write:
    if os.path.exists('./uvw_data/'):
        filename = 'uvw_data/uvws_%sMHz_%sMHz.npy' %(opts.freq, opts.freq_res)
    else:
        filename = 'uvws_%sMHz_%sMHz.npy' %(opts.freq, opts.freq_res)
    print 'Writing ' + filename + ' ...\n'

    # Only keep unique (u,v), i.e. assume perfect degeneracy
    np.save(filename, uvws)

if opts.plot:
    # Plotting
    fig = figure(figsize=(10,4.5))
    gs = gridspec.GridSpec(1,2)

    # Plot baseline(s) in position space
    ax1 = subplot(gs[0], aspect='equal')
    ax1.scatter(ant_pos[:,0], ant_pos[:,1], marker='o')
    ax1.set_xlabel('x [m]', size=16)
    ax1.set_ylabel('y [m]', size=16)

    ax2 = subplot(gs[1], aspect='equal')
    ax2.scatter(-uvws_meters[:, 0], -uvws_meters[:, 1], marker='.', c='r', label='Discarded')
    ax2.scatter(uvws_meters[:, 0], uvws_meters[:, 1], marker='.', c='g', label='Kept')
    ax2.set_xlabel(r'u [m]', size=16)
    ax2.set_ylabel(r'v [m]', size=16)
    ax2.legend(loc='upper center', frameon=False, ncol=2)

    xs = np.linspace(-uvws_meters[:, 0].max(), uvws_meters[:, 0].max(), 10)
    ax2.plot(xs, -xs, 'k--')

    for ax in fig.axes:
        ax.tick_params(axis='both',labelsize=16)

    gs.tight_layout(fig, rect=[None, None, None, 0.95])

    show()
