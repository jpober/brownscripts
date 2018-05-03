import numpy as np
import matplotlib.gridspec as gridspec
import time, optparse, sys, os

from matplotlib.pyplot import *



## ----------------- Option Parser ----------------- ##
o = optparse.OptionParser()
o.add_option('--visdata',
    type = str,
    help = 'Filename for input visibility binary .npy file.')
o.add_option('--plot_type',
    type = str,
    default = 'abs',
    help = 'Plot type, can be abs, phase, or real.')
o.add_option('--log_scale',
    action = 'store_true',
    default = False,
    help = 'If passed, plot in log10 for visibilities.')
o.add_option('--grid_pos',
    action='store_true',
    help = 'If passed, compare with grid centers. '
               +
               'Otherwise compare with true source positions.')
o.add_option('--multi-freq',
    action = 'store_true',
    help = 'If passed, perform analysis over multiple frequencies.')
opts,args = o.parse_args(sys.argv[1:])



# parse frequency information from filename
filename = opts.visdata
if '/' in filename:
    filename = opts.visdata.split('/')[-1]
freqs = np.array(map(float, filename.split('_')[1].strip('MHz').split('-')))
if not "zenith" in filename:
    freq_res = float(filename.split('_')[2].split('.npy')[0].strip('MHz'))
else:
    freq_res = float(filename.split('_')[2].strip('MHz'))
freqs = np.round(np.arange(freqs.min(), freqs.max() + freq_res, freq_res), decimals=3)



# Visibility function
def Vs_func(u, l, v, m):
    return np.exp(-2*np.pi*1j*(u*l+ v*m))
# Vs_func = lambda u,l,v,m: np.exp(-2*np.pi*1j*(u*l+ v*m))


# Constants
C = 3.e8 # meters per second


## ----------------- Read Data ----------------- ##
data_dic = np.load(opts.visdata).item()
nfreqs = 1
I = data_dic['sky']
true_pos = data_dic['catalog_true']
grid_pos = data_dic['catalog_grid']


# Need to consider gridded to gridded case
npix_side = 31
npix = npix_side**2
freq_ind = 0
wavelength = C/(freqs[freq_ind]*1.e6)

#Use 10 degrees on a side for the sky
# LOOK UP DEFINITION OF L AND M IN TMS

# Construct l,m grid
FOV = np.deg2rad(10) # 10 degree FOV
ls = np.linspace(-FOV/2, FOV/2, npix_side)
ms = np.copy(ls)
ls_vec, ms_vec = np.zeros(0), np.zeros(0)
for l in ls:
    for m in ms:
        ls_vec = np.append(ls_vec, l)
        ms_vec = np.append(ms_vec, m)

# Construct u,v grid
us_grid = np.fft.fftshift(np.fft.fftfreq(ls.shape[0], d=np.mean(np.diff(ls))))
vs_grid = np.fft.fftshift(np.fft.fftfreq(ms.shape[0], d=np.mean(np.diff(ms))))
us_vec, vs_vec = np.zeros(0), np.zeros(0)
for u in us_grid:
    for v in vs_grid:
        us_vec = np.append(us_vec, u)
        vs_vec = np.append(vs_vec, v)


# Use analytical solution to get visibilities using true positions
Vs = np.zeros((nfreqs, npix), dtype=complex)
for i in range(nfreqs):
    wavelength = C/(freqs[i]*1.e6)
    for j in range(us_vec.size):
        if opts.grid_pos:
            Vs[i, j] = np.sum(Vs_func(us_vec[j], ls[grid_pos[:, 0]],
                                                   vs_vec[j], ms[grid_pos[:, 1]]))
        else:
            Vs[i, j] = np.sum(Vs_func(us_vec[j], true_pos[:, 0],
                                                   vs_vec[j], true_pos[:, 1]))
Vs = Vs.reshape((nfreqs, npix_side, npix_side))

# Construct solution using analytic solution for maximum likelihood
# Assumes a Gaussian log likelihood function
# Requires noise injection into data (visibilities above)
a = np.zeros_like(Vs)
Vs_analytic = np.zeros_like(a)
RMS = 1.e-5
N_inv = np.eye(npix)/RMS**2

# Need to add the SAME noise to V(u,v) and V(-u, -v)
angles = np.angle(us_vec + 1j*vs_vec)
center_ind = np.where(np.logical_and(us_vec == 0.0, vs_vec == 0.0))[0]
top_right_inds = np.where(np.logical_and(angles > -np.pi/4, angles < 3*np.pi/4))[0]
top_right_inds = top_right_inds[top_right_inds != center_ind]
bottom_left_inds = np.where(np.logical_or(angles < -np.pi/4, angles > 3*np.pi/4))[0]
diag_inds = np.where(np.logical_or(angles == -np.pi/4, angles == 3*np.pi/4))[0]

# Create data from visibilities with injected Gaussian noise
d = np.copy(Vs)

for i in range(nfreqs):
    # Noise must be added so that d is still Hermitian

    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    # np.dot(DFT.conj().T, d.flatten()) not yielding a strictly real answer
    # Currently (4/17/18) noise is NOT being added correctly
    # Should test on a small subset of this data maybe 10pix on a side
    # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    d_flat = d[i].flatten()
    Vs_flat = Vs[i].flatten()
    # # Add noise to off diagonals
    # n = np.random.normal(0, RMS, top_right_inds.size)
    # d_flat[top_right_inds] += n
    # d_flat[bottom_left_inds] += n
    #
    # # Add noise to diagonal
    # n = np.random.normal(0, RMS, diag_inds.size)
    # d_flat[diag_inds] += n
    #
    # # Add noise to central pixel
    # n = np.random.normal(0, RMS, center_ind.size)
    # d_flat[center_ind] += n
    #
    # d[i] = d_flat.reshape((npix_side, npix_side))

    half_ind = int(us_vec.size/2.) + 1
    for j,[u,v] in enumerate(np.stack((us_vec, vs_vec), axis=1)[:half_ind]):
        neg_ind = np.where(np.logical_and(us_vec == -u, vs_vec == -v))[0][0]
        d_flat[[j, neg_ind]] += np.random.normal(0, RMS, 1)
        # print j, neg_ind, (u, v),
        # print (Vs_flat[j], Vs_flat[neg_ind]),
        # print (d_flat[j], d_flat[neg_ind])

    d[i] = d_flat.reshape([npix_side]*2)


# Iterate over frequencies
print 'Entering for loop...'
for i in range(nfreqs):
    DFT = np.exp(-1j*2*np.pi*(np.outer(us_vec, ls_vec)
                        +
                        np.outer(vs_vec, ms_vec)))
    inv_part = np.linalg.inv(np.dot(np.dot(DFT.conj().T, N_inv), DFT))
    right_part = np.dot(np.dot(DFT.conj().T, N_inv), d[i].flatten())
    a[i] = np.dot(inv_part, right_part).reshape((npix_side, npix_side))
    Vs_analytic[i] = np.dot(DFT, a[i].flatten()).reshape((npix_side, npix_side))
print 'For loop finished...'


## ----------------- Compute Power Spectrum ----------------- ##

# P(k_u,k_v) = sum |Visibilities(u,v)**2| for a given baseline length
U, V = np.meshgrid(us_grid, vs_grid)
rs = np.round(np.sqrt(U**2 + V**2), decimals=8)
rs_unique = np.unique(rs)


## ----------------- Perform FT Along Frequency Axis----------------- ##
Vs_ft = np.zeros_like(Vs)
for i in range(Vs.shape[1]):
    for j in range(Vs.shape[2]):
        Vs_ft[:, i, j] = np.fft.fftshift(np.fft.fft(Vs[:, i, j]))
ft_freqs = np.fft.fftshift(np.fft.fftfreq(nfreqs, d=1.e6*np.mean(np.diff(freqs))))
ft_freqs /= 1.e-9 # ns

# blns_unique = rs_unique*3.e8/(freqs[0]*1e6)
pspec = np.zeros((nfreqs, rs_unique.size))
for i in range(nfreqs):
    # rs = np.round(np.sqrt(uvs[i, 0, :, 0]**2 + uvs[i, 0, :, 1]**2), decimals=8)
    # rs_unique = np.unique(rs)
    for j,r in enumerate(rs_unique):
        r_inds = np.where(rs == r)
        pspec[i, j] = np.mean(np.abs(Vs_ft[i][r_inds])**2)

# FT analytic visibilities
Vs_aft = np.zeros_like(Vs_analytic)
for i in range(Vs_analytic.shape[1]):
    for j in range(Vs_analytic.shape[2]):
        Vs_aft[:, i, j] = np.fft.fftshift(np.fft.fft(Vs_analytic[:, i, j]))

# Compute delay power spectrum of analytic visibilities
pspec_analytic = np.zeros_like(pspec)
for i in range(nfreqs):
    # rs = np.round(np.sqrt(uvs[i, 0, :, 0]**2 + uvs[i, 0, :, 1]**2), decimals=8)
    # rs_unique = np.unique(rs)
    for j,r in enumerate(rs_unique):
        r_inds = np.where(rs == r)[0]
        pspec_analytic[i, j] = np.mean(np.abs(Vs_aft[i][r_inds])**2)





sys.exit()

# Old plotting with power spectrum

## ----------------- Plotting ----------------- ##
fig = figure(figsize = (15,6))
w = 10
gs = gridspec.GridSpec(2,5,
                                   height_ratios = [1, 10])#,
                                   # width_ratios = [w, w, w, w])
aspect = 'auto'

extent = [rs_unique.min(), rs_unique.max(), ft_freqs.min(), ft_freqs.max()]
# extent = [blns_unique.min(), blns_unique.max(), ft_freqs.min(), ft_freqs.max()]

master_ax = fig.add_subplot(gs[0,:])
master_ax.set_xticks([])
master_ax.set_yticks([])
master_ax.spines['top'].set_visible(False)
master_ax.spines['right'].set_visible(False)
master_ax.spines['bottom'].set_visible(False)
master_ax.spines['left'].set_visible(False)
ttl = master_ax.set_title(opts.visdata)
ttl.set_position([0.5, 0.5])

# Plot sky
skyax = subplot(gs[1, 0])
skyax.scatter(ls[grid_pos[:, 0]], ms[grid_pos[:, 1]], c='b', marker='o', label='Grid Position')
skyim = skyax.scatter(true_pos[:, 0], true_pos[:, 1], c='r', marker='.', label='Floating Position')
# skyax.legend(loc='best')
skyax.set_xlabel('l', size=16)
skyax.set_ylabel('m', size=16)
skyax.set_xlim([ls.min(), ls.max()])
skyax.set_ylim([ms.min(), ms.max()])

# # Shrink current axis by 20%
# box = skyax.get_position()
# skyax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
skyax.legend(loc='best')#, bbox_to_anchor=(1, 0.5))


# Plot numerical delay power spectrum
myax = subplot(gs[1, 1])
if opts.log_scale:
    vmin = np.log10(pspec).min()
    vmax = np.log10(pspec).max()
    if vmin < -10:
        vmin = -10
    myim = myax.imshow(np.log10(pspec),
                         interpolation='nearest',
                         origin='lower',
                         extent=extent,
                         aspect=aspect,
                         vmin=vmin, vmax=vmax)
    myax.set_title('Log Numerical (grid centers)', size=16)
else:
    print 'Not plotting in log scale'
    vmin = pspec.min()
    vmax = pspec.max()
    myim = myax.imshow(np.abs(pspec),
                         interpolation='nearest',
                         origin='lower',
                         extent=extent,
                         aspect=aspect,
                         vmin=vmin, vmax=vmax)
    myax.set_title('Numerical (grid centers)', size=16)


# Plot analytical delay power spectrum
anax = subplot(gs[1, 2])
if opts.log_scale:
    vmin = np.log10(pspec_analytic).min()
    vmax = np.log10(pspec_analytic).max()
    anim = anax.imshow(np.log10(pspec_analytic),
                        interpolation='nearest',
                        origin='lower',
                        extent=extent,
                        aspect=aspect,
                        vmin=vmin, vmax=vmax)
    if opts.grid_pos:
        anax.set_title('Log Analytic (grid centers)', size=16)
    else:
        anax.set_title('Log Analytic (true pos)', size=16)
else:
    vmin = pspec_analytic.min()
    vmax = pspec_analytic.max()
    anim = anax.imshow(pspec_analytic,
                        interpolation='nearest',
                        origin='lower',
                        extent=extent,
                        aspect=aspect,
                        vmin=vmin, vmax=vmax)
    if opts.grid_pos:
        anax.set_title('Analytic (grid centers)', size=16)
    else:
        anax.set_title('Analytic (true pos)', size=16)


# Plot difference of numerical and analytical delay power spectra
diffax = subplot(gs[1, 3])
diff_data = pspec - pspec_analytic
# diff_data = (pspec_analytic - pspec)/pspec_analytic
if opts.log_scale:
    vmin = np.log10(diff_data).min()
    vmax = np.log10(diff_data).max()
    diffim = diffax.imshow(np.log10(diff_data),
                         interpolation='nearest',
                         origin='lower',
                         extent=extent,
                         aspect=aspect,
                         vmin=vmin, vmax=vmax)
    diffax.set_title('Log(Grid - True)', size=16)
else:
    vmin = diff_data.min()
    vmax = diff_data.max()
    diffim = diffax.imshow(diff_data,
                         interpolation='nearest',
                         origin='lower',
                         extent=extent,
                         aspect=aspect,
                         vmin=vmin, vmax=vmax)
diffax.set_title('Grid - True', size=16)
# else:
#     vmin = np.log10(diff_data).min()
#     vmax = np.log10(diff_data).max()
#     diffim = diffax.imshow(diff_data,
#                          interpolation='nearest',
#                          origin='lower',
#                          extent=extent,
#                          aspect=aspect,
#                          vmin=vmin, vmax=vmax)
#     diffax.set_title('Log(Numerical - Analytic)', size=16)

fdiffax = subplot(gs[1, 4])
fdiff_data = (pspec - pspec_analytic)/pspec_analytic
# if opts.log_scale:
vmin = np.nanmin(np.log10(fdiff_data))
vmax = np.nanmax(np.log10(fdiff_data))
fdiffim = fdiffax.imshow(np.log10(fdiff_data),
                     interpolation='nearest',
                     origin='lower',
                     extent=extent,
                     aspect=aspect,
                     vmin=vmin, vmax=vmax)
fdiffax.set_title('Log[(Grid - True)/True]', size=16)
# else:
#     vmin = fdiff_data.min()
#     vmax = fdiff_data.max()
#     fdiffim = fdiffax.imshow(fdiff_data,
#                          interpolation='nearest',
#                          origin='lower',
#                          extent=extent,
#                          aspect=aspect,
#                          vmin=vmin, vmax=vmax)
fdiffax.set_title('(Grid - True)/True', size=16)

imgs = [myim, anim, diffim, fdiffim]
tau_max = rs_unique/3.e8*1e9
# tau_max = blns_unique/3.e8*1e9

for i,ax in enumerate(fig.axes[2:]):
    cb = fig.colorbar(imgs[i], ax=ax)#, format='%.1f')
    cb.ax.tick_params(labelsize=16)
    ax.set_xlabel(r'$\left| \vec{b} \right|$ [m]', size=16)
    ax.set_ylabel(r'$\tau$ [ns]', size=16)
    ax.plot(rs_unique, tau_max, 'w--', lw=2)
    ax.plot(rs_unique, -tau_max, 'w--', lw=2)
    # ax.plot(blns_unique, tau_max, 'w--', lw=2)
    # ax.plot(blns_unique, -tau_max, 'w--', lw=2)


gs.tight_layout(fig)

show()
#
# filename = 'delay_pspec_'+ '_'.join(opts.visdata.strip('.npy').split('_')[1:])
# savefig(filename+'.png')
