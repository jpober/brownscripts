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
opts,args = o.parse_args(sys.argv[1:])



## ----------------- Read Data ----------------- ##
data_dic = np.load(opts.visdata).item()
Vs = data_dic['vis']
uvs = data_dic['uv']
us = uvs[:, 0, :, 0]
vs = uvs[:, 0, :, 1]
I = data_dic['sky']
lms = data_dic['lm']
ls = lms[:, 0]
ms = lms[:, 1]
true_pos = data_dic['catalog_true']
grid_pos = data_dic['catalog_grid']
nfreqs = Vs.shape[0]
print Vs.shape

# parse frequency information from filename
if '/' in opts.visdata:
    opts.visdata = opts.visdata.split('/')[-1]
freqs = np.array(map(float, opts.visdata.split('_')[1].strip('MHz').split('-')))
if not "zenith" in opts.visdata:
    freq_res = float(opts.visdata.split('_')[2].split('.npy')[0].strip('MHz'))
else:
    freq_res = float(opts.visdata.split('_')[2].strip('MHz'))
freqs = np.round(np.arange(freqs.min(), freqs.max() + freq_res, freq_res), decimals=3)

# # Average over redundant baselines
# uvs_unique = np.unique(uvs, axis = 0)
# Vs_unique = np.zeros((nfreqs, uvs_unique.shape[0]))
#
# for u,v in uvs_unique:
#     avg_inds = np.where(np.logical_and(uvs[:, 0] == u, uvs[:, 1] == v))[0]
#     Vs_unique = np.append(Vs_unique, np.mean(Vs[avg_inds]))


## ----------------- Perform FT Along Frequency Axis----------------- ##
Vs_ft = np.zeros_like(Vs)
for i in range(Vs.shape[1]):
    Vs_ft[:, i] = np.fft.fftshift(np.fft.fft(Vs[:, i]))
ft_freqs = np.fft.fftshift(np.fft.fftfreq(nfreqs, d=1.e6*np.mean(np.diff(freqs))))
ft_freqs /= 1.e-9 # ns


## ----------------- Compute Power Spectrum ----------------- ##

# P(k_u,k_v) = sum |Visibilities(u,v)**2| for a given baseline length

rs = np.round(np.sqrt(uvs[0, 0, :, 0]**2 + uvs[0, 0, :, 1]**2), decimals=8)
rs_unique = np.unique(rs)
blns_unique = rs_unique*3.e8/(freqs[0]*1e6)
pspec = np.zeros((nfreqs, rs_unique.size))
for i in range(nfreqs):
    rs = np.round(np.sqrt(uvs[i, 0, :, 0]**2 + uvs[i, 0, :, 1]**2), decimals=8)
    rs_unique = np.unique(rs)
    for j,r in enumerate(rs_unique):
        r_inds = np.where(rs == r)[0]
        pspec[i, j] = np.mean(np.abs(Vs_ft[i, r_inds])**2)


## ----------------- Analytic solution comparison ----------------- ##
# Construct sky parameters
Vs_func = lambda u,l,v,m: np.exp(-2*np.pi*1j*(u*l+ v*m))

# Construct analytical solution
Vs_analytic = np.zeros((nfreqs, us.shape[1]), dtype=complex)
for i in range(nfreqs):
    for j in range(uvs.shape[2]):
        if opts.grid_pos:
            # Use nearest grid centers
            Vs_analytic[i, j] = np.sum(Vs_func(uvs[i, 0, j, 0], ls[grid_pos[:, 0]],
                                                               uvs[i, 0, j, 1], ms[grid_pos[:, 1]]))
        else:
            # Use true positions
            Vs_analytic[i, j] = np.sum(Vs_func(uvs[i, 0, j, 0], true_pos[:, 0],
                                                               uvs[i, 0, j, 1], true_pos[:, 1]))

# FT analytic visibilities
Vs_aft = np.zeros_like(Vs_analytic)
for i in range(Vs_analytic.shape[1]):
    Vs_aft[:, i] = np.fft.fftshift(np.fft.fft(Vs_analytic[:, i]))

# Compute delay power spectrum of analytic visibilities
pspec_analytic = np.zeros_like(pspec)
for i in range(nfreqs):
    rs = np.round(np.sqrt(uvs[i, 0, :, 0]**2 + uvs[i, 0, :, 1]**2), decimals=8)
    rs_unique = np.unique(rs)
    for j,r in enumerate(rs_unique):
        r_inds = np.where(rs == r)[0]
        pspec_analytic[i, j] = np.mean(np.abs(Vs_aft[i, r_inds])**2)


# Use analytical solution to get delay power spectrum using grid centers
#
# Vs = np.zeros((nfreqs, us.shape[1]), dtype=complex)
# for i in range(nfreqs):
#     for j in range(uvs.shape[2]):
#         if opts.grid_pos:
#             # Use nearest grid centers
#             Vs[i, j] = np.sum(Vs_func(uvs[i, 0, j, 0], ls[grid_pos[:, 0]],
#                                                    uvs[i, 0, j, 1], ms[grid_pos[:, 1]]))
#
# Vs_ft = np.zeros_like(Vs)
# for i in range(Vs.shape[1]):
#     Vs_ft[:, i] = np.fft.fftshift(np.fft.fft(Vs[:, i]))
# ft_freqs = np.fft.fftshift(np.fft.fftfreq(nfreqs, d=1.e6*np.mean(np.diff(freqs))))
# ft_freqs /= 1.e-9 # ns
#
# rs = np.round(np.sqrt(uvs[0, 0, :, 0]**2 + uvs[0, 0, :, 1]**2), decimals=8)
# rs_unique = np.unique(rs)
# blns_unique = rs_unique*3.e8/(freqs[0]*1e6)
# pspec = np.zeros((nfreqs, rs_unique.size))
# for i in range(nfreqs):
#     rs = np.round(np.sqrt(uvs[i, 0, :, 0]**2 + uvs[i, 0, :, 1]**2), decimals=8)
#     rs_unique = np.unique(rs)
#     for j,r in enumerate(rs_unique):
#         r_inds = np.where(rs == r)[0]
#         pspec[i, j] = np.mean(np.abs(Vs_ft[i, r_inds])**2)

# ------------------------------------------------------------------------------



## ----------------- Plotting ----------------- ##
fig = figure(figsize = (15,6))
w = 10
gs = gridspec.GridSpec(2,5,
                                   height_ratios = [1, 10])#,
                                   # width_ratios = [w, w, w, w])
aspect = 'auto'

# extent = [rs_unique.min(), rs_unique.max(), ft_freqs.min(), ft_freqs.max()]
extent = [blns_unique.min(), blns_unique.max(), ft_freqs.min(), ft_freqs.max()]

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
    # vmin = np.log10(pspec_analytic).min()
    # vmax = np.log10(pspec_analytic).max()
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
    # vmin = pspec_analytic.min()
    # vmax = pspec_analytic.max()
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
tau_max = blns_unique/3.e8*1e9

for i,ax in enumerate(fig.axes[2:]):
    cb = fig.colorbar(imgs[i], ax=ax)#, format='%.1f')
    cb.ax.tick_params(labelsize=16)
    ax.set_xlabel(r'$\left| \vec{b} \right|$ [m]', size=16)
    ax.set_ylabel(r'$\tau$ [ns]', size=16)
    ax.plot(blns_unique, tau_max, 'w--', lw=2)
    ax.plot(blns_unique, -tau_max, 'w--', lw=2)


gs.tight_layout(fig)

show()
#
# filename = 'delay_pspec_'+ '_'.join(opts.visdata.strip('.npy').split('_')[1:])
# savefig(filename+'.png')
