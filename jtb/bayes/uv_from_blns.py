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
    freqs = map(float, opts.freq.split(','))
elif '-' in opts.freq:
    freqs = map(float, opts.freq.split('-'))
    # freqs = np.linspace(freqs[0], freqs[1], int((freqs[1] - freqs[0])/opts.freq_res) + 1)
    freqs = np.arange(freqs[0], freqs[1] + opts.freq_res, opts.freq_res)
elif not (',' in opts.freq and '-' in opts.freq):
    freqs = [float(opts.freq)]

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
    blns = np.zeros(0)
    num_ants = ant_pos.shape[0]
    ant_nums = range(1, num_ants+1)
    for i in range(len(ant_nums)):
        ant1 = ant_nums[i]
        for j in range(len(ant_nums)):
            ant2 = ant_nums[j]
            if ant1 == ant2:
                continue
            if (uv.antnums_to_baseline(ant1, ant2) not in blns
                and uv.antnums_to_baseline(ant2, ant1) not in blns):
                blns = np.append(blns, uv.antnums_to_baseline(ant1, ant2))

# New line, better output formatting
print ''
for freq in freqs:
    # Get baseline
    l = c/(freq*1e6) #m, wavelength
    uvws = np.zeros((hours.shape[0], 3, blns.shape[0]))
    for i, bln in enumerate(blns):
        antpair = uv.baseline_to_antnums(bln)
        bln_vec = (ant_pos[antpair[1]-1,:] - ant_pos[antpair[0]-1,:])/l
        uvws[:,0,i] = np.sin(hours)*bln_vec[0] + np.cos(hours)*bln_vec[1]
        uvws[:,1,i] = -np.sin(decs)*(np.cos(hours)*bln_vec[0] - np.sin(hours)*bln_vec[1]) + np.cos(decs)*bln_vec[2]
        uvws[:,2,i] = np.cos(decs)*(np.cos(hours)*bln_vec[0] - np.sin(hours)*bln_vec[2]) + np.sin(decs)*bln_vec[2]

    if opts.write:
        filename = 'uvws_%sMHz.npy' %freq
        print 'Writing ' + filename + ' ...\n'
        np.save(filename, uvws)

    if opts.plot:
        # Plotting
        fig = figure(figsize=(10,4.5))
        gs = gridspec.GridSpec(1,2)

        # Plot baseline(s) in position space
        ax1 = subplot(gs[0], aspect='equal')
        ax1.plot(ant_pos[:,0],ant_pos[:,1],'o')
        ax1.set_xlabel('x [m]', size=18)
        ax1.set_ylabel('y [m]', size=18)

        ax2 = subplot(gs[1], aspect='equal')
        for i in range(blns.shape[0]):
            ax2.plot(uvws[:,0,i], uvws[:,1,i], 'b.')
            ax2.plot(-uvws[:,0,i], -uvws[:,1,i], 'b.')
        ax2.set_xlabel('u', size=18)
        ax2.set_ylabel('v', size=18)

        for ax in fig.axes:
            ax.tick_params(axis='both',labelsize=16)

        if opts.data:
            suptitle(opts.data + ', ' + str(freq) + ' MHz', size=16)
        else:
            suptitle(opts.positions + ', ' + str(freq) + ' MHz', size=16)

        gs.tight_layout(fig, rect=[None, None, None, 0.95])
        # tight_layout()

        show()
