import numpy as np
import healpy as hp
import matplotlib.gridspec as gridspec

from astropy.io import fits
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib.pyplot import *


fits_file = ('/Users/jburba/hera_things/hera-team/HERA-Beams/'
                +
                'NicolasFagnoniBeams/healpix_beam.fits')
beam_E = fits.getdata(fits_file, extname='BEAM_E')
beam_freqs = fits.getdata(fits_file, extname='FREQS')

nside = hp.npix2nside(beam_E.shape[0])
thetas, phis = hp.pix2ang(nside, np.arange(beam_E.shape[0]))
plot_phi_inds = np.where(np.logical_and(phis <= np.deg2rad(1), thetas <= np.pi/2.))[0]

nfreqs = 11
start_ind = np.where(beam_freqs == 150.0)[0][0]
colors = cm.jet_r(np.linspace(0, 1, nfreqs))
colors_phis = cm.jet_r(np.linspace(0, 1, plot_phi_inds.shape[0]))

# thetas_inds, freq_inds = np.meshgrid(plot_phi_inds, np.arange(start_ind, start_ind+nfreqs+1))
freq_inds, thetas_inds = np.meshgrid(np.arange(start_ind, start_ind+nfreqs), plot_phi_inds)
grid_data = beam_E[thetas_inds, freq_inds]

plot_freqs = beam_freqs[start_ind : start_ind + nfreqs+1]
extent = [plot_freqs.min(), plot_freqs.max(),
               np.rad2deg(thetas[plot_phi_inds]).min(),
               np.rad2deg(thetas[plot_phi_inds]).max()]

# fig = figure(figsize=(8,3))
# fig = figure()
# im = imshow(10*np.log10(grid_data), extent=extent, origin='upper', aspect='auto')
# xlabel('Frequency [MHz]', size=16)
# ylabel('Polar Angle [deg]', size=16)
# title('Cut of beam through phi=0', size=16)
# ax_divider = make_axes_locatable(gca())
# cax = ax_divider.append_axes("right", size="5%", pad="2%")
# cb = fig.colorbar(im, cax=cax)
# cb.set_label('dB', size=16)


# for i in range(start_ind, start_ind + nfreqs):
#     plot(np.rad2deg(thetas[plot_phi_inds]),
#           10*np.log10(beam_E[plot_phi_inds, i]),
#           label='%.0fMHz' %beam_freqs[i],
#           color=colors[i - 100])
# xlabel('Polar angle [deg]', size=16)
for i in range(plot_phi_inds.shape[0]):
    plot(beam_freqs[start_ind:start_ind + nfreqs],
           10*np.log10(beam_E[plot_phi_inds[i], start_ind:start_ind + nfreqs]),
           label='%.1f deg' %np.rad2deg(thetas[plot_phi_inds[i]]),
           color=colors_phis[i])
xlabel('Frequency [MHz]', size=16)
ylabel('dB', size=16)
tick_params(which='x', labelsize=16)
# legend(loc='best', ncol=2)
tight_layout()
show()
