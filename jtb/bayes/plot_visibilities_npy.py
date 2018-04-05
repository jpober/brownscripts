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
us = uvs[:, 0, :, 0]
vs = uvs[:, 0, :, 1]
I = data_dic['sky']
lms = data_dic['lm']


nfreqs = Vs.shape[0]
# extent_uv = [uvs[:, 0].min(), uvs[:, 0].max(), uvs[:, 1].min(), uvs[:, 1].max()]
print Vs.shape

if '/' in opts.visdata:
    opts.visdata = opts.visdata.split('/')[-1]

freqs = np.array(map(float, opts.visdata.split('_')[1].strip('MHz').split('-')))
if not "zenith" in opts.visdata:
    freq_res = float(opts.visdata.split('_')[2].split('.')[0].strip('MHz'))
else:
    freq_res = float(opts.visdata.split('_')[2].strip('MHz'))
freqs = np.round(np.arange(freqs.min(), freqs.max() + freq_res, freq_res), decimals=3)



## ----------------- Plotting ----------------- ##
nrows = int(np.sqrt(nfreqs))
ncols = int(np.ceil(float(nfreqs) / nrows))
nempty = nrows*ncols - nfreqs

fontsize = 16
figsize = (8, 6)

fig1 = figure(1, figsize=figsize)
fig2 = figure(2, figsize=figsize)
figs = [fig1, fig2]
gs = gridspec.GridSpec(nrows, ncols)

aspect = 'auto'
dmin, dmax = None, None

# Assign range for colorbar range if dmin and dmax not specified
if dmin is None and dmax is None:
    if opts.log_scale:
        if opts.plot_type == 'abs':
            dmin = np.ma.min(np.log10(np.abs(Vs)))
            dmax = np.ma.max(np.log10(np.abs(Vs)))
        elif opts.plot_type == 'real':
            dmin = np.ma.min(np.log10(np.real(Vs)))
            dmax = np.ma.max(np.log10(np.real(Vs)))
        elif opts.plot_type == 'imag':
            dmin = np.ma.min(np.log10(np.imag(Vs)))
            dmax = np.ma.max(np.log10(np.imag(Vs)))
        else:
            print 'Cannot have opts.log_scale with plot_type= \'phase\''
    else:
        if opts.plot_type == 'abs':
            dmin = np.ma.min(np.abs(Vs))
            dmax = np.ma.max(np.abs(Vs))
        elif opts.plot_type == 'real':
            dmin = np.ma.min(np.real(Vs))
            dmax = np.ma.max(np.real(Vs))
        elif opts.plot_type == 'imag':
            dmin = np.ma.min(np.imag(Vs))
            dmax = np.ma.max(np.imag(Vs))
        else:
            dmin = -np.pi
            dmax = np.pi

for fig in figs:
    plot_ind = [0, 0]
    freq_ind = 0

    # Add master labels for all x and y axes
    master_ax = fig.add_subplot(gs[:,:])
    if fig == fig1:
        master_ax.set_xlabel('l', size=fontsize, labelpad=25)
        master_ax.set_ylabel('m', size=fontsize, labelpad=35)
    else:
        master_ax.set_xlabel('u', size=fontsize, labelpad=25)
        master_ax.set_ylabel('v', size=fontsize, labelpad=35)

    master_ax.set_xticks([])
    master_ax.set_yticks([])
    master_ax.spines['top'].set_visible(False)
    master_ax.spines['right'].set_visible(False)
    master_ax.spines['bottom'].set_visible(False)
    master_ax.spines['left'].set_visible(False)
    if fig == fig1:
        ttl = master_ax.set_title(opts.visdata + ', Sky Realizations')
    else:
        if opts.log_scale:
            ttl = master_ax.set_title(opts.visdata + ', Log Visibilities')
        else:
            ttl = master_ax.set_title(opts.visdata + ', Visibilities')
    ttl.set_position([0.5, 1.05])

    while plot_ind[0] < nrows and freq_ind < nfreqs:
        while plot_ind[1] < ncols and freq_ind < nfreqs:
            ax = fig.add_subplot(gs[plot_ind[0], plot_ind[1]])
            if fig == fig1:
                im = ax.imshow(I[freq_ind], origin='center',
                                        aspect=aspect, interpolation='nearest')
                im.set_extent([-1, 1, -1, 1])
                im.set_clim(I.min(), I.max())
            else:
                if len(Vs.shape) > 2:
                    # univorm u,v sampling
                    if opts.plot_type == 'abs':
                        if opts.log_scale:
                            im = ax.imshow(np.log10(np.abs(Vs[freq_ind])), origin='center',
                                                    aspect=aspect, interpolation='nearest')
                        else:
                            im = ax.imshow(np.abs(Vs[freq_ind]), origin='center',
                                                    aspect=aspect, interpolation='nearest')
                    elif opts.plot_type == 'phase':
                        im = ax.imshow(np.angle(Vs[freq_ind]), origin='center',
                                                aspect=aspect, interpolation='nearest')
                    elif opts.plot_type == 'real':
                        if opts.log_scale:
                            im = ax.imshow(np.log10(np.real(Vs[freq_ind])), origin='center',
                                                    aspect=aspect, interpolation='nearest')
                        else:
                            im = ax.imshow(np.real(Vs[freq_ind]), origin='center',
                                                    aspect=aspect, interpolation='nearest')

                    # force clim on all subplots
                    im.set_clim(dmin, dmax)
                    im.set_extent(extent_uv)

                else:
                    # floating point u,v sampling
                    if opts.plot_type == 'abs':
                        if opts.log_scale:
                            im = ax.scatter(us[freq_ind], vs[freq_ind],
                                                    c = np.log10(np.abs(Vs[freq_ind])))
                        else:
                            im = ax.scatter(us[freq_ind], vs[freq_ind],
                                                    c = np.abs(Vs[freq_ind]))
                    elif opts.plot_type == 'phase':
                        im = ax.scatter(us[freq_ind], vs[freq_ind],
                                                c = np.angle(Vs[freq_ind]))
                    elif opts.plot_type == 'real':
                        if opts.log_scale:
                            im = ax.scatter(us[freq_ind], vs[freq_ind],
                                                    c = np.log10(np.real(Vs[freq_ind])))
                        else:
                            im = ax.scatter(us[freq_ind], vs[freq_ind],
                                                    c = np.real(Vs[freq_ind]))

                    # force clim on all subplots
                    im.set_clim(dmin, dmax)

            # remove axis labels from interior plots
            if nrows*ncols == nfreqs:
                if not plot_ind[0] == nrows-1:
                    ax.set_xticks([])
            else:
                if plot_ind[0] < nrows-2:
                    ax.set_xticks([])
                elif plot_ind[0] == nrows-2 and plot_ind[1] <= ncols-nempty-1:
                    ax.set_xticks([])

            if not plot_ind[1] == 0:
                ax.set_yticks([])

            # Add title with frequency
            axttl = ax.set_title(str(freqs[freq_ind]) + ' MHz', color='w')
            axttl.set_position([0.5, 0.05])
            # im.set_extent(extent_uv)
            # ax.set_xlabel('u', size=fontsize)
            # ax.set_ylabel('v', size=fontsize)
            # ax.tick_params(axis='both', labelsize=fontsize)
            #
            # # Add colorbar to each subplot
            # cb = fig.colorbar(im, ax=ax, format='%.1f')
            # cb.ax.tick_params(labelsize=fontsize)

            plot_ind[1] += 1
            freq_ind += 1

        # Reset column counter and update row counter
        plot_ind[1] = 0
        plot_ind[0] += 1

    gs.tight_layout(fig)

    # Append master colorbar
    gs.update(wspace=0.0, hspace=0.0, right=0.85)
    cbar_ax = fig.add_axes([0.875, 0.15, 0.02, 0.75])
    cb = fig.colorbar(im, cax=cbar_ax)
    if fig == fig2:
        cb.set_clim(dmin,dmax)
    else:
        cb.set_clim(I.min(), I.max())

    cb.ax.tick_params(labelsize=fontsize)

show()
