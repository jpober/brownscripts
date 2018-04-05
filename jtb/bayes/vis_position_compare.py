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
opts,args = o.parse_args(sys.argv[1:])



## ----------------- Read Data ----------------- ##
data_dic = np.load(opts.visdata).item()
Vs = data_dic['vis']
uvs = data_dic['uv']
I = data_dic['sky']
lms = data_dic['lm']
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


## ----------------- Analytic solution comparison ----------------- ##
# Construct sky parameters
ls = lms[:, 0]
ms = lms[:, 1]

# Point source, Flat beam
Vs_func = lambda u,l,v,m: np.exp(-2*np.pi*1j*(u*l+ v*m))
Vs_analytic_grid = np.zeros_like(Vs, dtype=complex)
Vs_analytic_float = np.zeros_like(Vs, dtype=complex)
for i in range(nfreqs):
    for j in range(uvs.shape[2]):
        Vs_analytic_float[i, j] = np.sum(Vs_func(uvs[i, 0, j, 0], true_pos[:, 0],
                                                           uvs[i, 0, j, 1], true_pos[:, 1]))
        Vs_analytic_grid[i, j] = np.sum(Vs_func(uvs[i, 0, j, 0], ls[grid_pos[:, 0]],
                                                           uvs[i, 0, j, 1], ms[grid_pos[:, 1]]))


# Compute difference of visibilities for first frequency
# Vs_diff_grid = Vs - Vs_analytic_grid
# Vs_diff_true = Vs - Vs_analytic_float

abs_Vs_diff_grid = np.abs(Vs) - np.abs(Vs_analytic_grid)
abs_Vs_diff_true = np.abs(Vs) - np.abs(Vs_analytic_float)
phase_Vs_diff_grid = np.angle(Vs) - np.angle(Vs_analytic_grid)
phase_Vs_diff_true = np.angle(Vs) - np.angle(Vs_analytic_float)

fig = figure(figsize=(15,5))
gs = gridspec.GridSpec(2, 3, height_ratios=[1,30])

master_ax = fig.add_subplot(gs[0,:])
master_ax.set_xticks([])
master_ax.set_yticks([])
master_ax.spines['top'].set_visible(False)
master_ax.spines['right'].set_visible(False)
master_ax.spines['bottom'].set_visible(False)
master_ax.spines['left'].set_visible(False)
ttl = master_ax.set_title(opts.visdata)
ttl.set_position([0.5, 0.5])


skyax = fig.add_subplot(gs[1, 0])
skyax.scatter(ls[grid_pos[:, 0]], ms[grid_pos[:, 1]], c='b', marker='o', label='Grid Position')
skyim = skyax.scatter(true_pos[:, 0], true_pos[:, 1], c='r', marker='.', label='Floating Position')
# skyax.legend(loc='best')
skyax.set_xlabel('l', size=16)
skyax.set_ylabel('m', size=16)
skyax.set_xlim([ls.min(), ls.max()])
skyax.set_ylim([ms.min(), ms.max()])


absax = fig.add_subplot(gs[1, 1])
for i in range(nfreqs):
    plot(abs_Vs_diff_grid[i], 'o')
    plot(abs_Vs_diff_true[i], '*')
    if i == range(nfreqs)[-1]:
        plot(abs_Vs_diff_grid[i], 'o', label='Grid Centers')
        plot(abs_Vs_diff_true[i], '*', label='True Pos')
absax.set_ylabel('| Numerical | - | Analytic |', fontsize=16, labelpad=10)
absax.legend(loc='lower right', fontsize=16, frameon=False)


phaseax = fig.add_subplot(gs[1, 2])
for i in range(nfreqs):
    plot(phase_Vs_diff_grid[i], 'o')
    plot(phase_Vs_diff_true[i], '*')
    if i == range(nfreqs)[-1]:
        plot(phase_Vs_diff_grid[i], 'o', label='Grid Centers')
        plot(phase_Vs_diff_true[i], '*', label='True Pos')
phaseax.set_ylabel('Phase(Numerical) - Phase(Analytic)', fontsize=16, labelpad=10)
# phaseax.legend(loc='lower right', fontsize=16, frameon=False)

# title('Single point source at zenith')
gs.tight_layout(fig)
show()

sys.exit()



## ----------------- Plotting ----------------- ##
fig = figure(figsize = (15, 6))
gs = gridspec.GridSpec(2, 2, height_ratios=[1, 10])
aspect = 'auto'

extent = [uvs[0, 0, :, 0].min(), uvs[0, 0, :, 0].max(),
               uvs[0, 0, :, 1].min(), uvs[0, 0, :, 1].max()]

master_ax = fig.add_subplot(gs[0, :])
master_ax.set_xticks([])
master_ax.set_yticks([])
master_ax.spines['top'].set_visible(False)
master_ax.spines['right'].set_visible(False)
master_ax.spines['bottom'].set_visible(False)
master_ax.spines['left'].set_visible(False)
ttl = master_ax.set_title(opts.visdata)
ttl.set_position([0.5, 0.5])

diff_grid_ax = subplot(gs[1, 0])
if opts.log_scale:
    vmin = np.log10(np.abs(Vs_diff_grid)).min()
    vmax = np.log10(np.abs(Vs_diff_grid)).max()
    diff_grid_im = diff_grid_ax.scatter(uvs[0, 0, :, 0], uvs[0, 0, :, 1],
                         c=np.log10(np.abs(Vs_diff_grid)),
                         vmin=vmin, vmax=vmax)
    diff_grid_ax.set_title('Log( Numerical - Analytic (grid centers) )', size=16)
else:
    vmin = np.abs(Vs_diff_grid).min()
    vmax = np.abs(Vs_diff_grid).max()
    diff_grid_im = diff_grid_ax.scatter(uvs[0, 0, :, 0], uvs[0, 0, :, 1],
                         c=np.abs(Vs_diff_grid),
                         vmin=vmin, vmax=vmax)
    diff_grid_ax.set_title('Numerical - Analytic (grid centers)', size=16)


diff_true_ax = subplot(gs[1,1])
if opts.log_scale:
    diff_true_im = diff_true_ax.scatter(uvs[0, 0, :, 0], uvs[0, 0, :, 1],
                        c=np.log10(np.abs(Vs_diff_true)),
                        vmin=vmin, vmax=vmax)
    diff_true_ax.set_title('Log( Numerical - Analytic (true pos) )', size=16)
else:
    diff_true_im = diff_true_ax.scatter(uvs[0, 0, :, 0], uvs[0, 0, :, 1],
                        c=np.abs(Vs_diff_true),
                        vmin=vmin, vmax=vmax)
    diff_true_ax.set_title('Numerical - Analytic (true pos)', size=16)


imgs = [diff_grid_im, diff_true_im]

for i,ax in enumerate(fig.axes[1:]):
    cb = fig.colorbar(imgs[i], ax=ax)#, format='%.1f')
    cb.ax.tick_params(labelsize=16)
    ax.set_xlabel(r'$u$', size=16)
    ax.set_ylabel(r'$v$', size=16)


gs.tight_layout(fig)

show()
