import numpy as np
import matplotlib.gridspec as gridspec
import time, optparse, sys, os
import astropy.cosmology as cosmo

from matplotlib.pyplot import *

## ----------------- Option Parser ----------------- ##
o = optparse.OptionParser()

o.add_option('--data',
    type = str,
    default = '',
    help = 'File path for previously generated rms data via this script.')

o.add_option('--k_vals',
    type = str,
    help = 'Filepath for file with k_vals sampled.')

o.add_option('--jy',
    action = 'store_true',
    help = 'If passed, compute dimensionless power spectrum based on data with input units of Jansky.')

opts,args = o.parse_args(sys.argv[1:])
print o.values

def comoving_distance(z):
    return cosmo.comoving_distance(z).value

def H(z):
    return cosmo.H(z).value

def redshift(nu):
    return v21/nu - 1

# Set up cosmology and parameters
cosmo = cosmo.FlatLambdaCDM(H0=100, Om0=0.308, Tcmb0=2.75, Ob0=0.0222)
H0 = cosmo.H0.value
c = 3.0e8 # m/s
v21 = 1.428e9 # Hz
FWHM = np.deg2rad(10) #radians
k_B = 1.38064852e-23

# Get k_vals
k_vals = np.genfromtxt(opts.k_vals) # Mpc/h

# Read in data
# data_dic = np.load(opts.data).item()
# data_array = data_dic['data_array'] # Visibilities in Janskys (ies?)
# freq_array = data_dic['freq_array'].squeeze()
# sigma_in = 17.096315755646163
# sigma_in = 18.191828652894863
sigma_in = 17
sigma_in *= 1./81
# sigma_in = np.std(data_array, axis=0).mean()
freq_array = np.arange(150, 160)
center_redshift = v21/freq_array.mean() - 1
B = freq_array[-1] - freq_array[0] # bandwidth
nf = freq_array.size

# Cosmology stuff
X = comoving_distance(center_redshift)
Y = comoving_distance(redshift(freq_array[0])) - comoving_distance(redshift(freq_array[-1]))

if opts.jy:
    # Convert Jy to K
    sigma_K_sq = (1.e-26*(c/freq_array.mean())**2/(2*k_B)*sigma_in)**2
else:
    sigma_K_sq = sigma_in**2
# Apply cosmology terms
sigma_K_sq = X**2*Y/(FWHM*B)*sigma_K_sq
# sigma_K_sq /= 512**2

if opts.jy:
    print 'Input [Jy]: %.1e' %sigma_in
else:
    print 'Input [K]: %.1e' %sigma_in
print 'Output [K^2]: %.1e' %sigma_K_sq

# Calculate dimensionless power spectrum estimates
dm_ps = (k_vals)**3/(2*np.pi)*sigma_K_sq

# Plotting
figure()
semilogy(np.log10(k_vals), dm_ps, 'o-')
xlabel(r'$\log(k)\ \left[h^{-1}\ \rm{Mpc}\right]$')
ylabel(r'$\Delta^2$')
title(str(sigma_in))
show()

sys.exit()

# Read in data
data_dic = np.load(opts.data).item()
vis_array = data_dic['data_array'] # Visibilities in Janskys (ies?)
uvw_array = data_dic['uvw_array'] # Phased uvw coordinates in meters
freq_array = data_dic['freq_array'].squeeze()
center_redshift = v21/freq_array.mean() - 1
B = freq_array[-1] - freq_array[0] # bandwidth
nf = freq_array.size

# Only keep first half of data (second half is complex conjugate)
half_ind = uvw_array.shape[0]/2
vis_array = vis_array[:half_ind]
uvw_array = uvw_array[:half_ind]
uvw_array *= freq_array.mean()/c # units of wavelengths

# Frequency axis
eta_array = np.fft.fftshift(np.fft.fftfreq(nf, d=np.diff(freq_array).mean()))
eta_pos_array = eta_array[eta_array >= 0]
k_pars = 2*np.pi*v21*H0*H(center_redshift)*eta_pos_array
k_pars /= c*(1 + center_redshift)**2

# Get k_perps corresponding to uvw array
u_mags = np.round(np.sqrt(uvw_array[:, 0]**2 + uvw_array[:, 1]**2), decimals=0)
unique_u_mags, u_mag_counts = np.unique(u_mags, return_counts=True)
Dz_center = comoving_distance(center_redshift) # Mpc/h
k_perps = 2*np.pi*unique_u_mags/Dz_center

# Power spectrum params
X = Dz_center
# Do you just use the whole bandwidth, then since each k_par
# is a combination of frequency ranges which have different comoving separations?
Y = comoving_distance(redshift(freq_array[0])) - comoving_distance(redshift(freq_array[-1]))

# Frequency FT data
vis_data_nudft = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(vis_array), axis=1))

# Multiply data by complex conjugate to get P(k_par, k_perp)
PS = vis_data_nudft*vis_data_nudft.conjugate()
PS *= X**2*Y/(FWHM*B)
PS_avg = np.zeros((k_perps.size, k_pars.size))
# for i, k_par in enumerate(k_pars):
#     PS[:, i]

# Make dimensionless power spectrum
K_perp, K_par = np.meshgrid(k_perps, k_pars)
ks = np.sqrt(K_par**2 + K_perp**2)
in_k_bins = np.zeros(())















sys.exit()


## ----------------- Plotting ----------------- ##
from matplotlib.pyplot import *
