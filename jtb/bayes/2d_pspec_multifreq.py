import numpy as np
import matplotlib.gridspec as gridspec
import time, optparse, sys, os

from matplotlib.pyplot import *
from mpl_toolkits.axes_grid1 import make_axes_locatable

## ----------------- Option Parser ----------------- ##
o = optparse.OptionParser()
o.add_option('--uvdata',
    type = str,
    help = 'Filename for input uvw binary .npy file.')
o.add_option('--freq',
    type=str,
    help='Frequency(ies) of observation in MHz.')
o.add_option('--freq_res',
    type = float,
    default = 1.,
    help = 'Channel width in MHz if --freq passed with \'-\'.')
o.add_option('--plot_type',
    type = str,
    default = 'abs',
    help = 'Plot type, can be abs, phase, or real.')
o.add_option('--write',
    action='store_true',
    help='If passed, write visibility array as .npy file.')
opts,args = o.parse_args(sys.argv[1:])

# Get frequency(ies) or frequency range
if ',' in opts.freq:
    freqs = np.sort(map(float, opts.freq.split(',')))
elif '-' in opts.freq:
    freqs = np.sort(map(float, opts.freq.split('-')))
    # freqs = np.linspace(freqs[0], freqs[1], int((freqs[1] - freqs[0])/opts.freq_res) + 1)
    freqs = np.arange(freqs[0], freqs[1] + opts.freq_res, opts.freq_res)
elif not (',' in opts.freq and '-' in opts.freq):
    freqs = [float(opts.freq)]

nfreqs = len(freqs)

## ----------------- Construct sky ----------------- ##
print 'Constructing sky...'
pixel_num = 150
ls = np.cos(np.linspace(-np.pi, np.pi, pixel_num))
ms = np.sin(np.linspace(-np.pi, np.pi, pixel_num))
pixel_side_length = np.diff(ls)[0]
npix = ls.size*ms.size
extent_lm = [ls.min(), ls.max(), ms.min(), ms.max()]
L, M = np.meshgrid(ls, ms)

#Make source catalog
nsources = 150
grid_pos = np.zeros((nsources,2), dtype = int)
true_pos = np.zeros((nsources,2))
for i in range(grid_pos.shape[0]):
    grid_pos[i, 0] = np.random.randint(0, ls.shape[0])
    grid_pos[i, 1] = np.random.randint(0, ms.shape[0])
    true_pos[i, 0] = ls[grid_pos[i, 0]] + np.random.uniform(low=-pixel_side_length, high=pixel_side_length)
    true_pos[i, 1] = ms[grid_pos[i, 1]] + np.random.uniform(low=-pixel_side_length, high=pixel_side_length)

# Make sky matrix
I = np.zeros((nfreqs, pixel_num, pixel_num))
I[:, grid_pos[:,0], grid_pos[:,1]] = 1.0

# Construct position vectors for npix pixels
I_vec = I.flatten()
ls_vec, ms_vec = np.zeros(0), np.zeros(0)
for l in ls:
    for m in ms:
        ls_vec = np.append(ls_vec, l)
        ms_vec = np.append(ms_vec, m)


## ----------------- Construct uv-plane ----------------- ##
print 'Constructing uv grid...'
us = np.fft.fftshift(np.fft.fftfreq(len(ls), d=np.mean(np.diff(ls))))
vs = np.fft.fftshift(np.fft.fftfreq(len(ms), d=np.mean(np.diff(ms))))
nvispix = us.size*vs.size
extent_uv = [us.min(), us.max(), vs.min(), vs.max()]

# Construct position vectors for nvispix pixels
us_vec, vs_vec = np.zeros(0), np.zeros(0)
for u in us:
    for v in vs:
        us_vec = np.append(us_vec, u)
        vs_vec = np.append(vs_vec, v)

print 'Constructing DFT matrix...'

start = time.time()
# Construct DFT matrix using outer product
DFT_mat = np.exp(-1j*2*np.pi*(np.outer(us_vec,ls_vec) + np.outer(vs_vec, ms_vec)))
# DFT_mat = np.zeros((DFT.shape[0]*nfreqs, DFT.shape[1]*nfreqs))
DFT = np.kron(np.eye(nfreqs), DFT_mat)
stop = time.time()
mytime = stop - start

print 'DFT constructed in %.1f s' %mytime

print 'Computing DFT...'

start = time.time()
Vs_vec = np.dot(DFT, I_vec)
stop = time.time()
mytime = stop - start

print 'DFT completed in %.1f s' %mytime

# Reshape visibilities
Vs = np.reshape(Vs_vec, (nfreqs, us.size, vs.size))

print Vs.shape

if opts.write:
    filename = 'visdata_%sMHz_%sMHz.npy' %(opts.freq, opts.freq_res)
    print 'Writing ' + filename + ' ...\n'
    np.save(filename, Vs)




## ----------------- Plotting ----------------- ##
nrows = int(np.sqrt(nfreqs))
ncols = int(np.ceil(float(nfreqs) / nrows))
nempty = nrows*ncols - nfreqs

fontsize = 16

# Plot sky for reference
figure(1)
imshow(I[0], origin='center', aspect='auto',
            interpolation='nearest', extent=extent_lm)
colorbar()
xlabel('l', size=fontsize)
ylabel('m', size=fontsize)
tick_params(axis='both', labelsize=fontsize)

fig = figure(2)
gs = gridspec.GridSpec(nrows, ncols)

plot_ind = [0, 0]
freq_ind = 0


while plot_ind[0] < nrows and freq_ind < nfreqs:
    while plot_ind[1] < ncols and freq_ind < nfreqs:
        ax = fig.add_subplot(gs[plot_ind[0], plot_ind[1]])
        if opts.plot_type == 'abs':
            im = ax.imshow(np.abs(Vs[freq_ind]), origin='center',
                                    aspect='auto', interpolation='nearest')
        elif opts.plot_type == 'phase':
            im = ax.imshow(np.angle(Vs[freq_ind]), origin='center',
                                    aspect='auto', interpolation='nearest')
        elif opts.plot_type == 'real':
            im = ax.imshow(np.real(Vs[freq_ind]), origin='center',
                                    aspect='auto', interpolation='nearest')
        ax.set_title(str(freqs[freq_ind]))
        im.set_extent(extent_uv)
        ax.set_xlabel('u', size=fontsize)
        ax.set_ylabel('v', size=fontsize)
        ax.tick_params(axis='both', labelsize=fontsize)

        # Add colorbar to each subplot
        cb = fig.colorbar(im, ax=ax, format='%.1f')
        cb.ax.tick_params(labelsize=fontsize)

        plot_ind[1] += 1

    # Reset column counter and update row counter
    plot_ind[1] = 0
    plot_ind[0] += 1

gs.tight_layout(fig)
show()

# start = time.time()
# npVs_0 = np.fft.fftshift(np.fft.fft2(I[0]))
# stop = time.time()
# nptime = stop - start
