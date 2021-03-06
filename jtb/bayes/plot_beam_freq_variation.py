import numpy as np
import healpy as hp
import matplotlib.gridspec as gridspec
import optparse, copy

from astropy.io import fits
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from matplotlib.colors import Normalize
from matplotlib.pyplot import *

## ----------------- Option Parser ----------------- ##
o = optparse.OptionParser()

o.add_option('--freq',
    type=str,
    help='Frequency(ies) of observation in MHz.')

o.add_option('--freq_res',
    type = float,
    default = 1.,
    help = 'Channel width in MHz if --freq passed with \'-\'.')

o.add_option('--beam',
    type = str,
    # default = ('/Users/jburba/hera_things/hera-team/HERA-Beams/'
    #                 +
    #                 'NicolasFagnoniBeams/healpix_beam.fits'),
    help = 'Healpix compatible fits file containing the primary beam.')

o.add_option('--poly_order',
    type = int,
    help = 'If n, fit beam variation as a function of frequency with an n-th order polynomial.')

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
    freqs = [float(opts.freq)]

nfreqs = len(freqs)


## ----------------- Main ----------------- ##

# fits_file = ('/Users/jburba/hera_things/hera-team/HERA-Beams/'
#                 +
#                 'NicolasFagnoniBeams/healpix_beam.fits')
beam_E = fits.getdata(opts.beam, extname='BEAM_E')
beam_freqs = fits.getdata(opts.beam, extname='FREQS')

nside = hp.npix2nside(beam_E.shape[0])
thetas, phis = hp.pix2ang(nside, np.arange(beam_E.shape[0]))
mean_diff = np.mean(np.diff(thetas)[np.diff(thetas) != 0.0])
fit_thetas = np.arange(0, np.pi + mean_diff, mean_diff)
fit_phis = np.zeros_like(fit_thetas)

start_ind = np.where(beam_freqs == freqs[0])[0][0]
# Normalize beam and get interpolated value along phi=0
beam_E /= beam_E[:, start_ind:start_ind + nfreqs + 1].max()
fit_beam = np.zeros((fit_thetas.size, freqs.size))
for i in range(freqs.size):
    fit_beam[:, i] = hp.get_interp_val(beam_E[:, start_ind + i], fit_thetas, fit_phis)

# Color scale stuff
norm = Normalize()
norm.autoscale(np.rad2deg(fit_thetas))
s_m = cm.ScalarMappable(cmap=cm.jet_r)
s_m.set_array([])
colors = s_m.to_rgba(np.linspace(0, 1, nfreqs))
colors_phis = s_m.to_rgba(np.linspace(0, 1, fit_thetas.size))

# Plotting stuff
plot_freqs = beam_freqs[start_ind : start_ind + nfreqs+1]
extent = [plot_freqs.min(),
               plot_freqs.max(),
               np.rad2deg(fit_thetas).min(),
               np.rad2deg(fit_thetas).max()]

fig = figure(figsize=(11,7))
for i in range(fit_thetas.size):
    if opts.poly_order:
        fit_xs = np.copy(beam_freqs[start_ind:start_ind + nfreqs])
        fit_ys = np.copy(fit_beam[i])
        coeffs = np.polyfit(fit_xs, fit_ys, opts.poly_order)
        quad_fit = np.zeros_like(fit_xs)
        for j in range(opts.poly_order + 1):
            quad_fit += coeffs[j]*fit_xs**(opts.poly_order - j)

        residuals = fit_ys - quad_fit
        plot(fit_xs, residuals,
               label='%.1f deg' %np.rad2deg(fit_thetas[i]),
               color=colors_phis[i])
    else:
        # plot(beam_freqs[start_ind:start_ind + nfreqs],
        #        beam_E[plot_phi_inds[i], start_ind:start_ind + nfreqs],
        #        label='%.1f deg' %np.rad2deg(thetas[plot_phi_inds[i]]),
        #        color=colors_phis[i])
        plot(beam_freqs[start_ind:start_ind + nfreqs],
               10*np.log10(fit_beam[i]),
               label='%.1f deg' %np.rad2deg(fit_thetas[i]),
               color=colors_phis[i])
        ylabel('dB', size=16)
xlabel('Frequency [MHz]', size=16)
tick_params(which='x', labelsize=16)
cb = colorbar(s_m)
cb.set_label(r'Polar Angle / $\pi$', size=16)

tight_layout()
show()







## ----------------- Extra Code ----------------- ##
# thetas_inds, freq_inds = np.meshgrid(plot_phi_inds, np.arange(start_ind, start_ind+nfreqs+1))
# freq_inds, thetas_inds = np.meshgrid(np.arange(start_ind, start_ind+nfreqs), plot_phi_inds)
# grid_data = beam_E[thetas_inds, freq_inds]

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

# colorbar(s_m)
# # Shrink current axis by 20%
# ax = gca()
# box = ax.get_position()
# ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
#
# # Put a legend to the right of the current axis
# ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=2)

# legend(loc='best', ncol=2)

# for i in range(start_ind, start_ind + nfreqs):
#     plot(np.rad2deg(thetas[plot_phi_inds]),
#           10*np.log10(beam_E[plot_phi_inds, i]),
#           label='%.0fMHz' %beam_freqs[i],
#           color=colors[i - 100])
# xlabel('Polar angle [deg]', size=16)
