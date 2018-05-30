import numpy as np
import healpy as hp
import matplotlib.gridspec as gridspec
import optparse, sys

from astropy.io import fits
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib.pyplot import *


## ----------------- Option Parser ----------------- ##
o = optparse.OptionParser()

o.add_option('--npix_side',
    type = int,
    default = 31,
    help = 'Number of pixels on a side for the lm and uv grids.')

o.add_option('--fov',
    type = float,
    default = 10,
    help = 'Field of view in degrees.')

opts,args = o.parse_args(sys.argv[1:])

# Construct l,m grid
NPIX_SIDE = opts.npix_side
FOV = np.deg2rad(opts.fov)
ls = np.linspace(-FOV/2, FOV/2, NPIX_SIDE)
ms = np.copy(ls)
extent_lm = [ls.min(), ls.max(), ms.min(), ms.max()]
ls_vec, ms_vec = np.zeros(0), np.zeros(0)
for m in ms:
    for l in ls:
        ls_vec = np.append(ls_vec, l)
        ms_vec = np.append(ms_vec, m)

thetas = np.sqrt(ls_vec**2 + ms_vec**2)
thetas *= FOV/(2.*thetas.max())
phis = np.arctan2(np.copy(ms_vec), np.copy(ls_vec))
phis[phis < 0.0] += 2*np.pi

# Get beam on my (l, m) grid
fits_file = ('/Users/jburba/hera_things/hera-team/HERA-Beams/'
                +
                'NicolasFagnoniBeams/healpix_beam.fits')
beam_E = fits.getdata(fits_file, extname='BEAM_E')
freqs = fits.getdata(fits_file, extname='FREQS')
freq_ind = np.where(freqs == 150.)[0][0]

beam_grid = hp.get_interp_val(beam_E[:, freq_ind], thetas, phis)

# Plotting
fig = figure(figsize=(16, 4))
gs = gridspec.GridSpec(1, 3)

thax = fig.add_subplot(gs[0])
thim = imshow(np.rad2deg(thetas).reshape((NPIX_SIDE, NPIX_SIDE)),
                        origin='lower',
                        extent=extent_lm)

phax = fig.add_subplot(gs[1])
phim = imshow(np.rad2deg(phis).reshape((NPIX_SIDE, NPIX_SIDE)),
                        origin='lower',
                        extent=extent_lm)

beamax = fig.add_subplot(gs[2])
beamim = imshow(10*np.log10(beam_grid).reshape((NPIX_SIDE, NPIX_SIDE)),
                            origin='lower',
                            extent=extent_lm)

imgs = [thim, phim, beamim]

for i, ax in enumerate(fig.axes):
    ax_divider = make_axes_locatable(ax)
    cax = ax_divider.append_axes("right", size="5%", pad="2%")
    cb = fig.colorbar(imgs[i], cax=cax)
    if imgs[i] == beamim:
        cb.set_label('dB', size=16)

    ax.set_xlabel('l')
    ax.set_ylabel('m')

gs.tight_layout(fig)
show()
