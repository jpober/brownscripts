import numpy as np
import matplotlib.gridspec as gridspec

from matplotlib.pyplot import *
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from scipy.interpolate import griddata, LinearNDInterpolator, interp2d, Rbf

# Create noise "image" cube
noise_rms = 1.
nl, nm, nf = 9, 9, 10
ls = np.arange(-nl/2 + 1, nl/2 + 1)
ms = ls.copy()
lm_pixel_half = np.diff(ls)[0]/2.
L, M = np.meshgrid(ls, ms)
ls_vec, ms_vec = L.flatten(), M.flatten()
ls_vec = ls_vec.reshape((1, -1))
ms_vec = ms_vec.reshape((1, -1))
# 81 pixels, 10 frequency channels
noise_cube = np.random.normal(0, noise_rms, (nf, nl, nl))
noise_cube_vec = noise_cube.reshape((nf, nl*nm))
npix = nl*nm*nf

# Create u,v grid for dft
us_dft = np.fft.fftshift(np.fft.fftfreq(nl))
vs_dft = us_dft.copy()
U_dft, V_dft = np.meshgrid(us_dft, vs_dft)
us_vec_dft, vs_vec_dft = U_dft.flatten(), V_dft.flatten()
us_vec_dft = us_vec_dft.reshape((-1, 1))
vs_vec_dft = vs_vec_dft.reshape((-1, 1))
uv_pixel_half = np.mean(np.diff(us_dft))/2.
nuv_dft = us_vec_dft.size
interp_coords = (us_vec_dft.squeeze(), vs_vec_dft.squeeze())
delta_u_dft = np.diff(us_dft).mean()

# Construct dft matrix
DFT = np.exp(-1j*2*np.pi*(np.dot(us_vec_dft, ls_vec)
                          +
                          np.dot(vs_vec_dft, ms_vec)))
DFTinv = np.exp(1j*2*np.pi*(np.outer(ls_vec, us_vec_dft)
                            +
                            np.outer(ms_vec, vs_vec_dft)))

#Perform 1000 iterations for random integer multiples of nuv_dft points for nuv_nudft
Nit = 1
hist_info = np.zeros((Nit, 3))
for iter_ind in range(Nit):
    print iter_ind + 1,
    # Create set of u,v floating point locations for nudft
    nratio_points = np.random.randint(1, 10)
    nuv_nudft = nratio_points*nuv_dft
    max_uv = us_dft.max()
    us_vec_nudft = np.random.uniform(-max_uv, max_uv, nuv_nudft)
    vs_vec_nudft = np.random.uniform(-max_uv, max_uv, nuv_nudft)
    us_vec_nudft = us_vec_nudft.reshape((-1, 1))
    vs_vec_nudft = vs_vec_nudft.reshape((-1, 1))
    points = (us_vec_nudft.squeeze(), vs_vec_nudft.squeeze())
    npad = 100
    pad_us = np.append(np.append(np.random.uniform(-max_uv - delta_u_dft, -max_uv, npad/4),
                       np.random.uniform(max_uv, max_uv + delta_u_dft, npad/4)),
                       np.random.uniform(-max_uv - delta_u_dft, max_uv + delta_u_dft, npad/2))
    pad_vs = np.append(np.random.uniform(-max_uv - delta_u_dft, max_uv + delta_u_dft, npad/2),
                       np.append(np.random.uniform(-max_uv - delta_u_dft, -max_uv, npad/4),
                       np.random.uniform(max_uv, max_uv + delta_u_dft, npad/4)))
    padded_us = np.append(us_vec_nudft.squeeze(), pad_us)
    padded_vs = np.append(vs_vec_nudft.squeeze(), pad_vs)
    padded_points = (padded_us, padded_vs)

    # Create nudft matrix
    NUDFT = np.exp(-1j*2*np.pi*(np.dot(us_vec_nudft, ls_vec)
                                +
                                np.dot(vs_vec_nudft, ms_vec)))
    NUDFTinv = np.exp(1j*2*np.pi*(np.outer(ls_vec, us_vec_nudft)
                                  +
                                  np.outer(ms_vec, vs_vec_nudft)))

    us_vec_nudft = us_vec_nudft.squeeze()
    vs_vec_nudft = vs_vec_nudft.squeeze()
    us_dft = us_dft.squeeze()
    vs_dft = vs_dft.squeeze()

    # Compute visibilities in two cases
    Vs_dft = np.zeros((nf, nuv_dft), dtype=complex)
    Vs_nudft = np.zeros((nf, nuv_nudft), dtype=complex)
    Vs_nudft_interp = np.zeros_like(Vs_dft)
    Vs_nudft_padded = np.zeros((nf, nuv_nudft + npad), dtype=complex)
    for i in range(nf):
        Vs_dft[i] = np.dot(DFT, noise_cube_vec[i])
        Vs_nudft[i] = np.dot(NUDFT, noise_cube_vec[i])
        # Vs_nudft_padded[i, :-npad] = Vs_nudft[i]
        # Vs_nudft_padded[i, -npad:] = (np.random.uniform(Vs_nudft[i].real.mean(), np.std(Vs_nudft[i]), npad)
        #                               +
        #                               np.random.uniform(0, np.std(Vs_nudft[i]), npad)*1j)
        # Interpolate Vs_nudft to us_dft, vs_dft grid before performing the inverse transform
        Vs_nudft_interp[i] = griddata(points,
                                      Vs_nudft[i],
                                      interp_coords,
                                      method = 'cubic',
                                      rescale = True)
        Vs_nudft_interp[np.isnan(Vs_nudft_interp)] = 0.0 + 1j*0.0
        # Vs_nudft_interp[i] = griddata(padded_points,
        #                               Vs_nudft_padded[i],
        #                               interp_coords,
        #                               method = 'linear')
        #
        #
        # # Interp2d
        # interp_func_real = interp2d(us_vec_nudft,
        #                             vs_vec_nudft,
        #                             Vs_nudft[i].real,
        #                             kind = 'linear',
        #                             fill_value = None)
        # interp_func_imag = interp2d(us_vec_nudft,
        #                             vs_vec_nudft,
        #                             Vs_nudft[i].imag,
        #                             kind = 'linear',
        #                             fill_value = None)
        # Vs_interp_real = interp_func_real(us_dft, vs_dft)
        # Vs_interp_imag = interp_func_imag(us_dft, vs_dft)
        # Vs_nudft_interp[i] = (Vs_interp_real + 1j*Vs_interp_imag).flatten()
        # rbfi = Rbf(us_vec_nudft, vs_vec_nudft, Vs_nudft[i], function='cubic')
        # Vs_nudft_interp[i] = rbfi(us_vec_dft, vs_vec_dft).squeeze()

    # Compare standard deviations of two cases
    # print 'Stddev Vs_dft: %.3e' %np.std(Vs_dft)
    # print 'Stddev Vs_nudft: %.3e' %np.std(Vs_nudft)

    # bad_inds = np.where(Vs_nudft_interp[0].imag == 0j)[0]
    # Vs_nudft_interp[:, bad_inds] = 0.0 + 0.0*1j

    # Get noise cube back from Visibilities
    noise_cube_dft = np.zeros((nf, nl*nm))
    noise_cube_nudft = np.zeros_like(noise_cube_dft)
    noise_cube_nudft_interp = np.zeros_like(noise_cube_dft)
    for i in range(nf):
        noise_cube_dft[i] = np.dot(DFTinv, Vs_dft[i]).real/nuv_dft
        noise_cube_nudft[i] = np.dot(NUDFTinv, Vs_nudft[i]).real/(nuv_nudft)
        noise_cube_nudft_interp[i] = np.dot(DFTinv, Vs_nudft_interp[i]).real/nuv_dft

    dft_stddev = np.std(noise_cube_dft)
    nudft_stddev = np.std(noise_cube_nudft)
    nudft_interp_stddev = np.std(noise_cube_nudft_interp)
    hist_info[iter_ind, 2] = nudft_interp_stddev/dft_stddev
    hist_info[iter_ind, 1] = nudft_stddev/dft_stddev
    hist_info[iter_ind, 0] = nratio_points
    # print 'Stddev noise_cube_dft: %.3e' %dft_stddev
    # print 'Stddev noise_cube_nudft: %.3e' %nudft_stddev
    # print 'Ratio nudft/dft: %.3f' %(nudft_stddev/dft_stddev)

    # Free up memory of N and Ninv
    # del(NUDFT, NUDFTinv)

# # write hist_info to numpy array
# outfile = 'stddev_ratios_%diterations' %Nit
# np.save(outfile, hist_info)
#
sys.exit()

# PLOTTING
fig = figure(figsize=(13,8))
gs = gridspec.GridSpec(2,4)

imgs = []
aspect = 'equal'

extent_lm = [ls.min() - lm_pixel_half, ls.max() + lm_pixel_half,
             ms.min() - lm_pixel_half, ms.max() + lm_pixel_half]
extent_uv = [us_dft.min() - uv_pixel_half, us_dft.max() + uv_pixel_half,
             vs_dft.min() - uv_pixel_half, vs_dft.max() + uv_pixel_half]

im_ax = subplot(gs[:, :2])
im = im_ax.imshow(noise_cube[0],
                  origin = 'lower',
                  extent = extent_lm,
                  aspect = aspect)
im_ax.set_title('Input')
im_ax.set_xlabel('l', size=16)
im_ax.set_ylabel('m', size=16, rotation=0)
imgs.append(im)

dft_ax = subplot(gs[0, 2])
im = dft_ax.imshow(np.abs(Vs_dft[0]).reshape([us_dft.size]*2),
                   origin = 'lower',
                   extent = extent_uv,
                   aspect = aspect)
dft_ax.set_title('Vs_dft')
imgs.append(im)

nudft_ax = subplot(gs[0, 3])
im = nudft_ax.scatter(us_vec_nudft[:, 0], vs_vec_nudft[:, 0],
                      c = np.abs(Vs_nudft[0]),
                      s = 2)
# x0,x1 = nudft_ax.get_xlim()
# y0,y1 = nudft_ax.get_ylim()
# nudft_ax.set_aspect(abs(x1-x0)/abs(y1-y0))
nudft_ax.set_aspect(1)
nudft_ax.set_title('Vs_nudft')
imgs.append(im)

im_dft_ax = subplot(gs[1, 2])
im = im_dft_ax.imshow(noise_cube_dft[0].reshape((nl, nm)),
                      origin = 'lower',
                      extent = extent_lm,
                      aspect = aspect)
im_dft_ax.set_title('Recovered Input (dft)')
imgs.append(im)

im_nudft_ax = subplot(gs[1, 3])
im = im_nudft_ax.imshow(noise_cube_nudft[0].reshape((nl, nm)),
                        origin = 'lower',
                        extent = extent_lm,
                        aspect = aspect)
im_nudft_ax.set_title('Recovered Input (nudft)')
imgs.append(im)

for i in np.arange(1, 5):
    a = fig.axes[i]
    if i < 3:
        a.set_xlabel('u', size=16)
        a.set_ylabel('v', size=16, rotation=0)
    else:
        a.set_xlabel('l', size=16)
        a.set_ylabel('m', size=16, rotation=0)

for i,axe in enumerate(fig.axes):
    divider = make_axes_locatable(axe)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    colorbar(imgs[i], cax=cax)

gs.tight_layout(fig)

show()
