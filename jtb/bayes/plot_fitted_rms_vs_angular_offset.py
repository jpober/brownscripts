import numpy as np
import matplotlib.gridspec as gridspec
import time, optparse, sys, os

from scipy.optimize import curve_fit
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

# Gaussian fitting function
def Gaussian(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


# Constants
C = 3.e8 # meters per second
NFREQS = 1
NPIX_SIDE = 31
NPIX = NPIX_SIDE**2

# Noise injection params
RMS = 1.e-6
N_inv = np.eye(NPIX)/RMS**2


# Construct l,m grid
FOV = np.deg2rad(10) # 10 degree FOV
ls = np.linspace(-FOV/2, FOV/2, NPIX_SIDE)
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

# Set up angular offsets, arrays to store fitted RMS values
angular_offsets = np.linspace(0, np.mean(np.diff(ls)), 10)
fitted_RMS = np.zeros_like(angular_offsets)
fitted_RMS_err = np.zeros_like(angular_offsets)

for offset_ind, angular_offset in enumerate(angular_offsets):
    true_pos = np.zeros((1, 2))
    true_pos[0, 0] = angular_offset

    # Use analytical solution to get visibilities using true positions
    Vs = np.zeros((NFREQS, NPIX), dtype=complex)
    for i in range(NFREQS):
        for j in range(us_vec.size):
            Vs[i, j] = np.sum(Vs_func(us_vec[j], true_pos[:, 0],
                                                   vs_vec[j], true_pos[:, 1]))
    Vs = Vs.reshape((NFREQS, NPIX_SIDE, NPIX_SIDE))

    # Construct solution using analytic solution for maximum likelihood
    # Assumes a Gaussian log likelihood function
    # Requires noise injection into data (visibilities above)
    a = np.zeros_like(Vs)
    Vs_maxL = np.zeros_like(a)

    # Create data from visibilities with injected Gaussian noise
    d = np.copy(Vs)

    for i in range(NFREQS):
        d_flat = d[i].flatten()
        Vs_flat = Vs[i].flatten()

        half_ind = int(us_vec.size/2.) + 1
        for j,[u,v] in enumerate(np.stack((us_vec, vs_vec), axis=1)[:half_ind]):
            neg_ind = np.where(np.logical_and(us_vec == -u, vs_vec == -v))[0][0]
            d_flat[[j, neg_ind]] += np.random.normal(0, RMS, 1)

        d[i] = d_flat.reshape([NPIX_SIDE]*2)


    # Iterate over frequencies
    print 'Entering for loop...'
    for i in range(NFREQS):
        DFT = np.exp(-1j*2*np.pi*(np.outer(us_vec, ls_vec)
                            +
                            np.outer(vs_vec, ms_vec)))
        inv_part = np.linalg.inv(np.dot(np.dot(DFT.conj().T, N_inv), DFT))
        right_part = np.dot(np.dot(DFT.conj().T, N_inv), d[i].flatten())
        a[i] = np.dot(inv_part, right_part).reshape((NPIX_SIDE, NPIX_SIDE))
        Vs_maxL[i] = np.dot(DFT, a[i].flatten()).reshape((NPIX_SIDE, NPIX_SIDE))
    print 'For loop finished...'

    # Compute fitted RMS
    diff_data = np.abs(Vs) - np.abs(Vs_maxL)
    counts, bins, patches = hist(diff_data[0].flatten(), bins=50)
    bin_width = np.mean(np.diff(bins))
    fit_xs = bins[:-1] + bin_width/2
    guess_params = [np.max(counts), 0.0, RMS]
    # fit_params: 0, amplitude; 1, mean; 2, std dev
    fit_params, fit_cov = curve_fit(Gaussian, fit_xs, counts, p0=guess_params)

    # Store fitted RMS and fit error
    fitted_RMS[offset_ind] = fit_params[2]
    fitted_RMS_err[offset_ind] = fit_cov[-1, -1]

# Plotting
plot_xs = angular_offsets/np.mean(np.diff(ls))
errorbar(plot_xs, fitted_RMS/RMS, yerr=fitted_RMS_err/RMS)
xlabel('Pixel Offset', size=16)
ylabel('Fitted RMS/$10^{-5}$', size=16)
ylim([(fitted_RMS/RMS).min(), (fitted_RMS/RMS).max()])
tick_params(which='both', labelsize=16)
tight_layout()
show()
