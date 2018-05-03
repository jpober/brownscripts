import numpy as np
import matplotlib.gridspec as gridspec
import time, optparse, sys, os

from matplotlib.pyplot import *

## ----------------- Option Parser ----------------- ##
o = optparse.OptionParser()
o.add_option('--uvdata',
    type = str,
    help = 'Filename for input uvw binary .npy file.')
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
o.add_option('--nsources',
    type = int,
    default = 1,
    help= 'Number of sources to add to image')
o.add_option('--zenith_source',
    action = 'store_true',
    help = 'If passed, only place one source at zenith.')
o.add_option('--horizon_source',
    action = 'store_true',
    help = 'If passed, only place one source at the horizon.')
o.add_option('--spec_index',
    type = float,
    default = '0.0',
    help = 'Spectral index for amplitude of point sources as a function of frequency,'
           + 'i.e. n = spectral index and amplitude A(v)=v^n.')
opts,args = o.parse_args(sys.argv[1:])


# Get frequency(ies) or frequency range
# parse frequency information from filename
if '/' in opts.uvdata:
    filename = opts.uvdata.split('/')[-1]
else:
    filename = opts.uvdata
freqs_arr = np.array(map(float, filename.split('_')[1].strip('MHz').split('-')))
if not "zenith" in filename:
    freq_res = float(filename.split('_')[2].split('.npy')[0].strip('MHz'))
else:
    freq_res = float(filename.split('_')[2].strip('MHz'))
freqs = np.round(np.arange(freqs_arr.min(), freqs_arr.max() + freq_res, freq_res), decimals=3)
freqs = freqs[np.where(freqs <= freqs_arr[1])]

nfreqs = len(freqs)

#
# if ',' in opts.freq:
#     freqs = np.sort(map(float, opts.freq.split(',')))
# elif '-' in opts.freq:
#     freqs = np.sort(map(float, opts.freq.split('-')))
#     # freqs = np.linspace(freqs[0], freqs[1], int((freqs[1] - freqs[0])/opts.freq_res) + 1)
#     freqs = np.arange(freqs[0], freqs[1] + opts.freq_res, opts.freq_res)
# elif not (',' in opts.freq and '-' in opts.freq):
#     freqs = [float(opts.freq)]
#
nfreqs = len(freqs)
print 'Frequency range: %s' %', '.join(map(str, [freqs.min(), freqs.max()]))
print 'Frequency resolution: %.3f' %freq_res


## ----------------- Construct sky ----------------- ##
print 'Constructing sky...'
pixel_num = 31
ls = np.linspace(-1, 1, pixel_num)
ms = np.linspace(-1, 1, pixel_num)
pixel_side_length = np.diff(ls)[0]/2.
nlm_pix = ls.size*ms.size
L, M = np.meshgrid(ls, ms)
extent_lm = [ls.min(), ls.max(), ms.min(), ms.max()]

if len(ls) % 2 == 0:
    mid_l = len(ls)/2
else:
    mid_l = int(len(ls)/2. - 0.5)
if len(ms) % 2 == 0:
    mid_m = len(ms)/2
else:
    mid_m = int(len(ms)/2. - 0.5)

#Make source catalog
if opts.zenith_source:
    nsources = 1
    grid_pos = np.zeros((nsources, 2), dtype = int)
    true_pos = np.zeros((nsources, 2))
    for i in range(grid_pos.shape[0]):
        grid_pos[i, 0] = mid_l
        grid_pos[i, 1] = mid_m
        true_pos[i, 0] = ls[grid_pos[i, 0]] + np.random.uniform(low=-pixel_side_length,
                                                                                          high=pixel_side_length)
        true_pos[i, 1] = ms[grid_pos[i, 1]] + np.random.uniform(low=-pixel_side_length,
                                                                                            high=pixel_side_length)

elif opts.horizon_source:
    nsources = 1
    grid_pos = np.zeros((nsources, 2), dtype = int)
    true_pos = np.zeros((nsources, 2))
    for i in range(grid_pos.shape[0]):
        grid_pos[i, 0] = mid_l
        grid_pos[i, 1] = 0
        true_pos[i, 0] = ls[grid_pos[i, 0]] + np.random.uniform(low=-pixel_side_length,
                                                                                          high=pixel_side_length)
        true_pos[i, 1] = ms[grid_pos[i, 1]] + np.random.uniform(low=-pixel_side_length,
                                                                                            high=pixel_side_length)

else:
    nsources = opts.nsources
    grid_pos = np.zeros((nsources, 2), dtype = int)
    true_pos = np.zeros((nsources, 2))
    for i in range(grid_pos.shape[0]):
        grid_pos[i, 0] = np.random.randint(0, ls.shape[0])
        grid_pos[i, 1] = np.random.randint(0, ms.shape[0])
        true_pos[i, 0] = ls[grid_pos[i, 0]] + np.random.uniform(low=-pixel_side_length,
                                                                                          high=pixel_side_length)
        true_pos[i, 1] = ms[grid_pos[i, 1]] + np.random.uniform(low=-pixel_side_length,
                                                                                            high=pixel_side_length)

# Make sky matrix
I = np.zeros((nfreqs, pixel_num, pixel_num))
for i, freq in enumerate(freqs):
    if not opts.spec_index == 0.0:
        I[i, grid_pos[:,0], grid_pos[:,1]] = 1./(1 + (freq - freqs.min())**opts.spec_index)
    else:
        I[i, grid_pos[:,0], grid_pos[:,1]] = 1.



## ----------------- Construct uv-plane ----------------- ##
print 'Loading uvdata...'
uvs = np.load(opts.uvdata)
# uvs only contains unique (u,v) assuming perfect degeneracy
# uvs has shape (nfreqs, ntimes, nblns, 3)
us = uvs[:, 0, :, 0]
vs = uvs[:, 0, :, 1]
nuv_pix = us.shape[1]*vs.shape[1]



## ----------------- Analytic solution comparison ----------------- ##
# Point source, Flat beam
Vs_func = lambda u,l,v,m: np.exp(-2*np.pi*1j*(u*l+ v*m))
Vs_analytic = np.zeros((nfreqs, us.shape[1]), dtype=complex)
for i in range(nfreqs):
    for j in range(uvs.shape[2]):
        if not opts.grid_pos:
            Vs_analytic[i, j] = np.sum(Vs_func(uvs[i, 0, j, 0], true_pos[:, 0],
                                                               uvs[i, 0, j, 1], true_pos[:, 1]))
        else:
            Vs_analytic[i, j] = np.sum(Vs_func(uvs[i, 0, j, 0], ls[grid_pos[:, 0]],
                                                               uvs[i, 0, j, 1], ms[grid_pos[:, 1]]))

# FT along frequency axis
Vs_aft = np.zeros_like(Vs_analytic)
for i in range(Vs_analytic.shape[1]):
    Vs_aft[:, i] = np.fft.fftshift(np.fft.fft(Vs_analytic[:, i]))
# Vs_aft = np.fft.fftshift(np.fft.fft(Vs_analytic, axis=0))
# Could need to use fft2 or fftn
ft_freqs = np.fft.fftshift(np.fft.fftfreq(nfreqs, d=1.e6*np.mean(np.diff(freqs))))
ft_freqs /= 1.e-9 # ns

rs = np.round(np.sqrt(uvs[0, 0, :, 0]**2 + uvs[0, 0, :, 1]**2), decimals=8)
rs_unique = np.unique(rs)
blns_unique = rs_unique*3.e8/(freqs[0]*1e6)
pspec_analytic = np.zeros((nfreqs, rs_unique.size))
for i in range(nfreqs):
    rs = np.round(np.sqrt(uvs[i, 0, :, 0]**2 + uvs[i, 0, :, 1]**2), decimals=8)
    rs_unique = np.unique(rs)
    for j,r in enumerate(rs_unique):
        r_inds = np.where(rs == r)[0]
        pspec_analytic[i, j] = np.mean(np.abs(Vs_aft[i, r_inds])**2)



## ----------------- Plotting ----------------- ##
fig = figure(figsize = (15,6))
gs = gridspec.GridSpec(2,3, height_ratios=[1,10], width_ratios=[3,1,3])
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
ttl = master_ax.set_title(opts.uvdata)
ttl.set_position([0.5, 0.5])

skyax = subplot(gs[1,0])
# if opts.log_scale:
#     vmin = np.log10(np.abs(I[0])).min()
#     vmax = np.log10(np.abs(I[0])).max()
#     myim = skyax.imshow(np.log10(np.abs(I[0])),
#                          interpolation='nearest',
#                          origin='lower',
#                          extent=extent_lm,
#                          aspect=aspect,
#                          vmin=vmin, vmax=vmax)
#     skyax.set_title('Log Sky', size=16)
# else:
#     vmin = np.abs(I[0]).min()
#     vmax = np.abs(I[0]).max()
#     skyim = skyax.imshow(np.abs(I[0]),
#                          interpolation='nearest',
#                          origin='lower',
#                          extent=extent_lm,
#                          aspect=aspect,
#                          vmin=vmin, vmax=vmax)
#     skyax.set_title('Sky', size=16)

# # Overlay sources at true positions
skyax.scatter(ls[grid_pos[:, 0]], ms[grid_pos[:, 1]], c='b', marker='o', label='Grid Position')
skyim = skyax.scatter(true_pos[:, 0], true_pos[:, 1], c='r', marker='.', label='Floating Position')
# skyax.legend(loc='best')
skyax.set_xlabel('l', size=16)
skyax.set_ylabel('m', size=16)
skyax.set_xlim([ls.min(), ls.max()])
skyax.set_ylim([ms.min(), ms.max()])

# Shrink current axis by 20%
box = skyax.get_position()
skyax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
skyax.legend(loc='center left', bbox_to_anchor=(1, 0.5))


anax = subplot(gs[1,-1])
if opts.log_scale:
    vmin = np.log10(np.abs(pspec_analytic)).min()
    vmax = np.log10(np.abs(pspec_analytic)).max()
    anim = anax.imshow(np.log10(np.abs(pspec_analytic)),
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

cb = fig.colorbar(anim, ax=anax)#, format='%.1f')
cb.ax.tick_params(labelsize=16)
anax.set_xlabel(r'$\left| \vec{b} \right|\ $ [m]', size=16)
anax.set_ylabel(r'$\tau$ [ns]', size=16)

# Plot tau_max line
tau_max = blns_unique/3.e8*1e9
anax.plot(blns_unique, tau_max, 'w--', lw=2)
anax.plot(blns_unique, -tau_max, 'w--', lw=2)

gs.tight_layout(fig)

show()
#
# filename = 'delay_pspec_'+ '_'.join(opts.uvdata.strip('.npy').split('_')[1:])
# savefig(filename+'.png')
