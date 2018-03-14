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
I = data_dic['sky']
true_pos = data_dic['catalog_true']
grid_pos = data_dic['catalog_grid']
nfreqs = Vs.shape[0]
extent_uv = [uvs[:, 0].min(), uvs[:, 0].max(), uvs[:, 1].min(), uvs[:, 1].max()]
U, V = np.meshgrid(uvs[:,0], uvs[:,1])
print Vs.shape

freqs = np.array(map(float, opts.visdata.split('_')[1].strip('MHz').split('-')))
freq_res = float(opts.visdata.split('_')[2].split('.')[0].strip('MHz'))
freqs = np.arange(freqs.min(), freqs.max() + freq_res, freq_res)

# # Average over redundant baselines
# uvs_unique = np.unique(uvs, axis = 0)
# Vs_unique = np.zeros((nfreqs, uvs_unique.shape[0]))
#
# for u,v in uvs_unique:
#     avg_inds = np.where(np.logical_and(uvs[:, 0] == u, uvs[:, 1] == v))[0]
#     Vs_unique = np.append(Vs_unique, np.mean(Vs[avg_inds]))


## ----------------- Perform FT Along Frequency Axis----------------- ##
Vs_ft = np.fft.fftshift(np.fft.fft(Vs, axis=0))
ft_freqs = np.fft.fftshift(np.fft.fftfreq(nfreqs, d=1.e6*np.mean(np.diff(freqs))))
ft_freqs /= 1.e-9 # ns


## ----------------- Compute Power Spectrum ----------------- ##
# P(k_u,k_v) = sum_r(Visibilities(u,v)**2) s.t. sqrt(u**2+v**2) = r
# Averaged over circular annuli

rs = np.round(np.sqrt(U**2 + V**2), decimals=14)
rs_unique = np.unique(rs)
pspec = np.zeros((nfreqs, rs_unique.shape[0]))*1.0j
# pspec[:,0] gives all frequency fft amplitudes at r=rs_unique[0]

for i,r in enumerate(rs_unique):
    r_inds = np.where(rs == r)
    pspec[:, i] = np.mean(np.abs(Vs_ft[:, r_inds[0], r_inds[1]])**2, axis=1)


## ----------------- Analytic solution comparison ----------------- ##
# Construct sky parameters
ls = np.linspace(-1, 1, I.shape[1])
ms = np.linspace(-1, 1, I.shape[2])

# Point source, Flat beam
Vs_func = lambda u,l,v,m: np.exp(-2*np.pi*1j*(u*l+ v*m))
Vs_analytic = np.zeros_like(Vs, dtype=complex)
for i in range(Vs.shape[0]):
    for j, u in enumerate(uvs[:, 0]):
        for k, v in enumerate(uvs[:, 1]):
            if not opts.grid_pos:
                Vs_analytic[i, j, k] = np.sum(Vs_func(u, true_pos[:,0],
                                                                       v, true_pos[:,1]))
            else:
                Vs_analytic[i, j, k] = np.sum(Vs_func(u, ls[grid_pos[:,0]],
                                                                    v, ms[grid_pos[:,1]]))

Vs_aft = np.fft.fftshift(np.fft.fft(Vs_analytic, axis=0))
pspec_analytic = np.zeros_like(pspec)*1.0j
for i,r in enumerate(rs_unique):
    inds = np.where(rs == r)
    pspec_analytic[:, i] = np.mean(np.abs(Vs_aft[:, r_inds[0], r_inds[1]])**2, axis=1)


## ----------------- Plotting ----------------- ##
fig = figure(figsize = (15,6))
gs = gridspec.GridSpec(2,3, height_ratios=[1,10])
aspect = 'auto'

extent = [rs_unique.min(), rs_unique.max(), ft_freqs.min(), ft_freqs.max()]

master_ax = fig.add_subplot(gs[0,:])
master_ax.set_xticks([])
master_ax.set_yticks([])
master_ax.spines['top'].set_visible(False)
master_ax.spines['right'].set_visible(False)
master_ax.spines['bottom'].set_visible(False)
master_ax.spines['left'].set_visible(False)
ttl = master_ax.set_title(opts.visdata)
ttl.set_position([0.5, 0.5])

myax = subplot(gs[1,0])
if opts.log_scale:
    vmin = np.log10(np.abs(pspec)).min()
    vmax = np.log10(np.abs(pspec)).max()
    myim = myax.imshow(np.log10(np.abs(pspec)),
                         interpolation='nearest',
                         origin='lower',
                         extent=extent,
                         aspect=aspect,
                         vmin=vmin, vmax=vmax)
    myax.set_title('Log Numerical (grid centers)', size=16)
else:
    vmin = np.abs(pspec).min()
    vmax = np.abs(pspec).max()
    myim = myax.imshow(np.abs(pspec),
                         interpolation='nearest',
                         origin='lower',
                         extent=extent,
                         aspect=aspect,
                         vmin=vmin, vmax=vmax)
    myax.set_title('Numerical (grid centers)', size=16)


anax = subplot(gs[1,1])
if opts.log_scale:
    vmin = np.log10(np.abs(pspec_analytic)).min()
    vmax = np.log10(np.abs(pspec_analytic)).max()
    anim = anax.imshow(np.abs(pspec_analytic),
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
    vmin = np.abs(pspec_analytic).min()
    vmax = np.abs(pspec_analytic).max()
    anim = anax.imshow(np.abs(pspec_analytic),
                        interpolation='nearest',
                        origin='lower',
                        extent=extent,
                        aspect=aspect,
                        vmin=vmin, vmax=vmax)
    if opts.grid_pos:
        anax.set_title('Analytic (grid centers)', size=16)
    else:
        anax.set_title('Analytic (true pos)', size=16)


diffax = subplot(gs[1,2])
diff_data = np.abs(pspec) - np.abs(pspec_analytic)
if opts.log_scale:
    vmin = np.log10(diff_data).min()
    vmax = np.log10(diff_data).max()
    diffim = diffax.imshow(np.log10(diff_data),
                         interpolation='nearest',
                         origin='lower',
                         extent=extent,
                         aspect=aspect,
                         vmin=vmin, vmax=vmax)
    diffax.set_title('Log(Numerical - Analytic)', size=16)
else:
    vmin = diff_data.min()
    vmax = diff_data.max()
    diffim = diffax.imshow(diff_data,
                         interpolation='nearest',
                         origin='lower',
                         extent=extent,
                         aspect=aspect,
                         vmin=vmin, vmax=vmax)
    diffax.set_title('Numerical - Analytic', size=16)

imgs = [myim, anim, diffim]

for i,ax in enumerate(fig.axes[1:]):
    cb = fig.colorbar(imgs[i], ax=ax)#, format='%.1f')
    cb.ax.tick_params(labelsize=16)
    ax.set_xlabel(u'\"$k_\perp$\"', size=16)
    ax.set_ylabel(u'$t$ [ns]', size=16)


gs.tight_layout(fig)

show()
