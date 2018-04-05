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

o.add_option('--spec_index',
    type = float,
    default = '0.0',
    help = 'Spectral index for amplitude of point sources as a function of frequency,'
           + 'i.e. n = spectral index and amplitude A(v)=v^n.')

o.add_option('--zenith_source',
    action = 'store_true',
    help = 'If passed, only place one source at zenith.')

o.add_option('--horizon_source',
    action = 'store_true',
    help = 'If passed, only place one source at the horizon.')

o.add_option('--nsources',
    type = int,
    default = 1,
    help= 'Number of sources to add to image')

o.add_option('--long_bl_ew',
    type = float,
    default = 84,
    help = 'Longest East-West baseline for array configuration in meters. '
               +
               'Defaults to 84m for HERA-37 compact hex.')

o.add_option('--long_bl_ns',
    type = float,
    default = 72.75,
    help = 'Longest North-South baseline for array configuration in meters. '
               +
               'Defaults to 84m for HERA-37 compact hex.')

# o.add_option('--plot',
#     action = 'store_true',
#     help = 'If passed, show plot of simulated visibilities.')

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

## ----------------- Construct uv-plane ----------------- ##
if not opts.uvdata:
    print 'Please supply a .npy file with uv data.'
    # print 'Constructing uv grid...'
    # nuv_pix = 31
    # u_max = opts.long_bl_ew*(freqs.max()*1e6)/3.e8
    # v_max = opts.long_bl_ns*(freqs.max()*1e6)/3.e8
    # us = np.zeros((nfreqs, nuv_pix))
    # us[] = np.linspace(-u_max, u_max, nuv_pix)
    # vs = np.linspace(-v_max, v_max, nuv_pix)
    # # us = np.fft.fftshift(np.fft.fftfreq(31, d=np.mean(np.diff(ls))))
    # # vs = np.fft.fftshift(np.fft.fftfreq(31, d=np.mean(np.diff(ms))))
    # nvispix = us.size*vs.size
else:
    print 'Loading uvdata...'
    uvs = np.load(opts.uvdata)
    # uvs only contains unique (u,v) assuming perfect degeneracy
    # uvs has shape (nfreqs, ntimes, nblns, 3)
    us = uvs[:, 0, :, 0]
    vs = uvs[:, 0, :, 1]
    nuv_pix = us.shape[1]*vs.shape[1]
    # us and vs have shape (nfreqs, nblns, 1)

# Get max u and v for setting sky resolution
u_max = np.abs(us).max()
v_max = np.abs(vs).max()


# print 'nuv_pix: %d' %nuv_pix
# # Construct position vectors for nvispix pixels
# us_vec = np.zeros((nfreqs, nuv_pix))
# vs_vec = np.zeros((nfreqs, nuv_pix))
# for i in range(len(freqs)):
#     nuv_ind = 0
#     for j in range(us.shape[1]):
#         for k in range(vs.shape[1]):
#             us_vec[i, nuv_ind] = us[i, j]
#             vs_vec[i, nuv_ind] = vs[i, k]
#             nuv_ind += 1


## ----------------- Construct sky ----------------- ##
print 'Constructing sky...'
pixel_num = 31
# Sky resolution set by max u and v
uv_max = np.max([u_max, v_max])
if not pixel_num%2 == 0:
    l_max = (pixel_num-1)/(2*uv_max)
    m_max = (pixel_num-1)/(2*uv_max)
else:
    l_max = pixel_num/(2*uv_max)
    m_max = pixel_num/(2*uv_max)
ls = np.linspace(-l_max, l_max, pixel_num)
ms = np.linspace(-m_max, m_max, pixel_num)
pixel_side_length = np.diff(ls)[0]
nlm_pix = ls.size*ms.size
L, M = np.meshgrid(ls, ms)
extent_lm = [ls.min(), ls.max(), ms.min(), ms.max()]

#Make source catalog
nsources = opts.nsources
if len(ls) % 2 == 0:
    mid_l = len(ls)/2
else:
    mid_l = int(len(ls)/2. - 0.5)
if len(ms) % 2 == 0:
    mid_m = len(ms)/2
else:
    mid_m = int(len(ms)/2. - 0.5)

if opts.zenith_source:
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
        true_pos[i, 1] = ms[grid_pos[i, 1]] + np.random.uniform(low=0,
                                                                                            high=pixel_side_length)

else:
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

# Construct position vectors for npix pixels
I_vec = I.flatten()
ls_vec, ms_vec = np.zeros(0), np.zeros(0)
for l in ls:
    for m in ms:
        ls_vec = np.append(ls_vec, l)
        ms_vec = np.append(ms_vec, m)



## ----------------- DFT Things ----------------- ##
# print 'Constructing DFT matrix...'

# Construct DFT matrix using outer product
# DFT_mat = np.exp(-1j*2*np.pi*(np.outer(us_vec,ls_vec) + np.outer(vs_vec, ms_vec)))
# print 'Size of DFT_mat in GB: ' + str(sys.getsizeof(DFT_mat)/1.e9)
# DFT_mat = np.zeros((DFT.shape[0]*nfreqs, DFT.shape[1]*nfreqs))
# DFT = np.kron(np.eye(nfreqs), DFT_mat)

# print 'DFT constructed in %.1f s' %mytime

# Construct DFT using u,v per frequency
# DFT = np.zeros((nfreqs*nuv_pix, nfreqs*nlm_pix), dtype=complex)
# for i in range(nfreqs):
#     DFT[i*nuv_pix:(i+1)*nuv_pix, i*nlm_pix:(i+1)*nlm_pix] = (
#             np.exp(-1j*2*np.pi*(np.outer(us_vec[i], ls_vec) + np.outer(vs_vec[i], ms_vec)))
#             )


# Perform DFT per frequency using a for loop
print "Performing DFT via for loop..."
Vs = np.zeros((nfreqs, us.shape[1]), dtype=complex)
start = time.time()

for i in range(nfreqs):
    DFT = np.exp(-1j*2*np.pi*(np.outer(us[i], ls_vec) + np.outer(vs[i], ms_vec)))
    Vs[i] = np.dot(DFT, I[i].flatten()) #.reshape((us.shape[1], vs.shape[1]))

stop = time.time()
mytime = stop - start

print 'DFT completed in %.1f s' %mytime
# print 'Size of DFT in GB: ' + str(sys.getsizeof(DFT)/1.e9)

# print 'Computing DFT...'
#
# start = time.time()
# Vs_vec = np.dot(DFT, I_vec)
# stop = time.time()
# mytime = stop - start
#
# print 'DFT completed in %.1f s' %mytime

# # Reshape visibilities
# Vs = np.reshape(Vs_vec, (nfreqs, us.size, vs.size))

print 'Visibilities array shape: %s' %str(Vs.shape)

if opts.write:
    if os.path.exists('./sim_vis/'):
        filename = 'sim_vis/visdata_%sMHz_%sMHz' %(opts.freq, opts.freq_res)
    else:
        filename = 'visdata_%sMHz_%sMHz' %(opts.freq, opts.freq_res)
    if not opts.spec_index == 0.0:
        filename += '_spec-ind_%.1f' %opts.spec_index
    if opts.zenith_source:
        filename += '_zenith'
    elif opts.horizon_source:
        filename += '_horizon'
    if opts.nsources > 1:
        filename += '_%dsources' %opts.nsources
    print 'Writing ' + filename + '.npy ...\n'
    out_dic = {}
    out_dic['vis'] = Vs
    out_dic['uv'] = uvs
    out_dic['sky'] = I
    out_dic['lm'] = np.stack((ls, ms)).T
    out_dic['catalog_grid'] = grid_pos
    out_dic['catalog_true'] = true_pos
    np.save(filename + '.npy', out_dic)


# if opts.plot:
#     ## ----------------- Plotting ----------------- ##
#     nrows = int(np.sqrt(nfreqs))
#     ncols = int(np.ceil(float(nfreqs) / nrows))
#     nempty = nrows*ncols - nfreqs
#
#     fontsize = 16
#
#     # Plot sky for reference
#     figure(1)
#     imshow(I[0], origin='center', aspect='auto',
#                 interpolation='nearest', extent=extent_lm)
#     colorbar()
#     xlabel('l', size=fontsize)
#     ylabel('m', size=fontsize)
#     tick_params(axis='both', labelsize=fontsize)
#
#     fig = figure(2)
#     gs = gridspec.GridSpec(nrows, ncols)
#
#     plot_ind = [0, 0]
#     freq_ind = 0
#
#     while plot_ind[0] < nrows and freq_ind < nfreqs:
#         while plot_ind[1] < ncols and freq_ind < nfreqs:
#             ax = fig.add_subplot(gs[plot_ind[0], plot_ind[1]])
#             if opts.plot_type == 'abs':
#                 im = ax.imshow(np.abs(Vs[freq_ind]), origin='center',
#                                         aspect='equal', interpolation='nearest')
#             elif opts.plot_type == 'phase':
#                 im = ax.imshow(np.angle(Vs[freq_ind]), origin='center',
#                                         aspect='equal', interpolation='nearest')
#             elif opts.plot_type == 'real':
#                 im = ax.imshow(np.real(Vs[freq_ind]), origin='center',
#                                         aspect='equal', interpolation='nearest')
#             ax.set_title(str(freqs[freq_ind]))
#             im.set_extent([us[i].min(), us[i].max(), vs[i].min(), vs[i].max()])
#             ax.set_xlabel('u', size=fontsize)
#             ax.set_ylabel('v', size=fontsize)
#             ax.tick_params(axis='both', labelsize=fontsize)
#
#             # Add colorbar to each subplot
#             cb = fig.colorbar(im, ax=ax, format='%.1f')
#             cb.ax.tick_params(labelsize=fontsize)
#
#             plot_ind[1] += 1
#
#         # Reset column counter and update row counter
#         plot_ind[1] = 0
#         plot_ind[0] += 1
#
#     gs.tight_layout(fig)
#     show()
