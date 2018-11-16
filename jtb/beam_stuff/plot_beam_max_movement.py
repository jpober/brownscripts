import numpy as np
import healpy as hp
import time, optparse, sys, os, re

from astropy.io import fits
from scipy import interpolate


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
    # default = ('/Users/jburba/hera_things/hera-team/HERA-Beams/'
    #                 +
    #                 'NicolasFagnoniBeams/healpix_beam.fits'),
    help = 'Healpix compatible fits file containing the primary beam.')

o.add_option('--nside',
    type = 'int',
    default = 64,
    help = 'Nside parameter for healpix map.')

o.add_option('--linear',
    action = 'store_true',
    help = 'If passed, use a linear nearest neighbor interpolation.')

o.add_option('--use_cst',
    action = 'store_true',
    help = 'If passed, use CST output, not interpolated output.')

opts,args = o.parse_args(sys.argv[1:])
nside = opts.nside

if opts.beam:
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
else:
    filenames = args
    print 'Reading in %d files...' %len(filenames)
    freqs = np.array([float(re.findall(r'\d+', f.split('_')[-1])[0]) for f in filenames])
    order = np.argsort(freqs)
    freqs = freqs[order]
    nfreqs = len(freqs)
    filenames = np.array(filenames)[order]

    nfiles = len(filenames)

    if opts.use_cst:
        data = np.loadtxt(filenames[0], skiprows=2, usecols=[0, 1])
        thetas, phis = np.deg2rad(data[:, 0]), np.deg2rad(data[:, 1])
        npix = thetas.size
        beam_E = np.zeros((npix, nfiles))
        for fi, f in enumerate(filenames):
            beam_E[:, fi] = np.loadtxt(f, skiprows=2, usecols=3)

        max_inds = np.argmax(beam_E, axis=0)
        max_locs = np.zeros((nfreqs, 2))
        max_locs[:, 0] = thetas[max_inds]
        max_locs[:, 1] = phis[max_inds]

    else:
        beam_E = np.zeros((hp.nside2npix(nside), len(freqs)))
        thetai, phii = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))

        for fi, f in enumerate(filenames):
            data = np.loadtxt(f, skiprows=2)
            if opts.linear:
                lat = data[:, 0] * np.pi / 180.0
                nlat = len(lat)
                lon = data[:, 1] * np.pi / 180.0
                nlon = len(lon)
                coords = np.stack((lat, lon)).T
                gain = data[:, 2]
                lut = interpolate.LinearNDInterpolator(coords, gain)
            else:
                lat = np.unique(data[:, 0]) * np.pi / 180.0
                nlat = len(lat)
                lon = np.unique(data[:, 1]) * np.pi / 180.0
                nlon = len(lon)
                gain = data[:, 2].reshape(nlon, nlat).transpose()
                lut = interpolate.RectBivariateSpline(lat, lon, gain)
            for i in np.arange(hp.nside2npix(nside)):
                beam_E[i, fi] = lut(thetai[i], phii[i])

        max_inds = np.argmax(beam_E, axis=0)
        max_locs = np.zeros((nfreqs, 2))
        max_locs[:, 0], max_locs[:, 1] = hp.pix2ang(nside, max_inds)

        # Write to .fits file
        if opts.linear:
            filename = 'healpix_beam_%.0f-%.0fMHz_nside-%d_linear-interp.fits' %(freqs[0], freqs[-1], nside)
        else:
            filename = 'healpix_beam_%.0f-%.0fMHz_nside-%d.fits' %(freqs[0], freqs[-1], nside)
        new_hdul = fits.HDUList()
        new_hdul.append(fits.ImageHDU(data=beam_E, name='BEAM_E'))
        new_hdul.append(fits.ImageHDU(data=freqs, name='FREQS'))
        print 'Writing ' + filename
        new_hdul.writeto(filename)

max_locs_lm = np.zeros_like(max_locs)
max_locs_lm[:, 0] = np.cos(max_locs[:, 1])*np.sin(max_locs[:, 0])
max_locs_lm[:, 1] = np.sin(max_locs[:, 1])*np.sin(max_locs[:, 0])
max_locs_lm = np.rad2deg(max_locs_lm)
max_locs = np.rad2deg(max_locs)

# ---------------------------------- Plotting ---------------------------------- #
import matplotlib.gridspec as gridspec

from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib.pyplot import *

fig = figure(figsize=(9, 4))
gs = GridSpec(2, 2)

scat_ax = fig.add_subplot(gs[:, 0])
im = scat_ax.scatter(max_locs_lm[:, 0], max_locs_lm[:, 1], c=freqs, marker='o', linestyle='-')
scat_ax.set_xlabel('l [deg]', size=16)
scat_ax.set_ylabel('m [deg]', size=16)
scat_ax.grid(which='both')

lat_ax = fig.add_subplot(gs[0, 1])
lat_ax.scatter(freqs, max_locs[:, 0], linestyle='-', c='k', marker='.')
lat_ax.set_ylabel('Theta [deg]', size=16)
lat_ax.set_xticks([])

lon_ax = fig.add_subplot(gs[1, 1])
lon_ax.scatter(freqs, max_locs[:, 1], linestyle='-', c='k', marker='.')
lon_ax.set_ylabel('Phi [deg]', size=16)
lon_ax.set_xlabel('Frequency [MHz]', size=16)

ax_divider = make_axes_locatable(scat_ax)
cbar_ax = ax_divider.append_axes("right", size="5%", pad="2%")
cb = fig.colorbar(im, cax=cbar_ax)
cb.set_label('Frequency [MHz]', size=16)

for ax in fig.axes:
    ax.set_aspect('auto')

gs.tight_layout(fig)

gs.update(top=0.9, hspace=0., wspace=0.6)

if opts.use_cst:
    fig.suptitle('Raw CST Output')
else:
    fig.suptitle('Nside: %d' %nside)


show()
