import numpy as np
import healpy as hp
import matplotlib.gridspec as gridspec
import time, optparse, sys, os

from astropy.io import fits
from matplotlib.colors import Normalize
from matplotlib.pyplot import *

## ----------------- Option Parser ----------------- ##
o = optparse.OptionParser()

o.add_option('--freq',
    type=str,
    default='50-250',
    help='Frequency(ies) of observation in MHz.  Nicolas\' beams run from 50-250MHz.' )

o.add_option('--freq_res',
    type = float,
    default = 1.,
    help = 'Channel width in MHz if --freq passed with \'-\'.')

o.add_option('--beam',
    type = 'str',
    default = ('/Users/jburba/hera_things/hera-team/HERA-Beams/'
                    +
                    'NicolasFagnoniBeams/healpix_beam.fits'),
    help = 'Healpix compatible fits file containing the primary beam.')

opts,args = o.parse_args(sys.argv[1:])


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


# --------------------------------- Beam Stuff --------------------------------- #
beam_freqs = fits.getdata(opts.beam, extname='FREQS')
start_ind = np.where(beam_freqs == freqs[0])[0][0]

beam_E = fits.getdata(opts.beam, extname='BEAM_E')[:, start_ind:start_ind + nfreqs]
beam_freqs = beam_freqs[start_ind:start_ind + nfreqs]
nside = hp.npix2nside(beam_E.shape[0])

max_inds = np.argmax(beam_E, axis=0)
max_locs = np.zeros((nfreqs, 2))
max_locs[:, 0], max_locs[:, 1] = hp.pix2ang(nside, max_inds)

max_locs_lm = np.zeros_like(max_locs)
# max_locs_lm[:, 0] = max_locs[:, 0]*np.cos(max_locs[:, 1])
# max_locs_lm[:, 1] = max_locs[:, 0]*np.sin(max_locs[:, 1])
max_locs_lm[:, 0] = np.cos(max_locs[:, 1])*np.sin(max_locs[:, 0])
max_locs_lm[:, 1] = np.sin(max_locs[:, 1])*np.sin(max_locs[:, 0])
max_locs_lm = np.rad2deg(max_locs_lm)
# for freq_ind in range(nfreqs):
#     max_ind = np.argmax(beam_E[:, freq_ind])
#     max_loc= np.rad2deg(hp.pixind2ang(nside, max_ind))

# ---------------------------------- Plotting ---------------------------------- #

# Color scale stuff
norm = Normalize()
norm.autoscale(freqs)
s_m = cm.ScalarMappable(cmap=cm.jet_r)
s_m.set_array([])
colors = s_m.to_rgba(np.arange(freqs[0], freqs[-1] + 1, 1))

figure()
scatter(max_locs_lm[:, 0], max_locs_lm[:, 1], c=colors, marker='o', linestyle='-')
xlabel('l [deg]', size=16)
ylabel('m [deg]', size=16)
cb = colorbar(s_m)
cb.set_label('Frequency [MHz]', size=16)

grid(which='both', axis='both')
scatter(max_locs_lm[:, 0], max_locs_lm[:, 1], c=colors, marker='o', linestyle='-')

show()
