#! /usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import aipy as a, numpy as n
from matplotlib import pylab as p
import capo as C
import sys, optparse, re, os
import glob

### Plot the |k| power spectra from npz files, each a different color set of points.


o=optparse.OptionParser()
o.add_option('--beam', action='store_true',
                help='scale data by beam square instead of beam')
o.add_option('--cov', action='store', type='float', default=None,
            help='scale factor for signal loss in covariance removal')
o.add_option('--med', action='store_true',
            help='correction factor for using median statistics instead of the mean')
o.add_option('--show', action='store_true',
            help='Show the plot')
opts,args = o.parse_args(sys.argv[1:])
print args

def mean_temp(z):
    return 28. * ((1.+z)/10.)**.5 # mK

#fig=p.figure(figsize=(12,7.2))
fig=p.figure()

colors = 'kbcm' * 10
 
#tau_h = 100 + 15. #in ns
#k_h = C.pspec.dk_deta(C.pspec.f2z(.151)) * tau_h
#p.subplot(121)
#p.vlines(k_h, -1e7, 1e8, linestyles='--', linewidth=1.5)
#p.vlines(-k_h, -1e7, 1e8, linestyles='--', linewidth=1.5)
##p.gca().set_yscale('log', nonposy='clip')
#p.xlabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$', fontsize='large')
#p.ylabel(r'$P(k)[\ {\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$',fontsize='large')
#p.ylim(-.6e7,1.75e7)
##p.ylim(-.6e7,1.4e14)
##p.ylim(1e5,5e16)
#p.grid()


tau_h = 100 + 15. #in ns
k_h = C.pspec.dk_deta(C.pspec.f2z(.151)) * tau_h

p.subplot(111)

#theo_noise = noise_level(freq=freq)

afreqs=[]

maxpk = 0.0

for j,filename in enumerate(args):
    print 'Reading', filename
    f = n.load(filename)
    afreqs.append(f['freq'])
    ks,pk,err = f['k'], f['k3pk'], f['k3err']
    maxpk = n.max([maxpk,n.max(pk)])
    p.errorbar(ks,pk,yerr=err, fmt=colors[j]+'.',label=os.path.basename(filename))
    #k0 = n.argmin(n.abs(ks))

max_y = 0.5e2*maxpk
p.vlines(k_h, -1e7, max_y, linestyles='--', linewidth=1.5)

freq=n.average(afreqs)

#Plot theoretical EoR curve
filename = glob.glob('/users/alanman/capo/ael/lidz_mcquinn_k3pk/*7.32.dat')[0]
#filename = glob.glob('/users/alanman/capo/ael/lidz_mcquinn_k3pk/power_21cm_z11.46.dat')[0]
print 'Reading', filename
d = n.array([map(float, L.split()) for L in open(filename).readlines()])
ks, pk = d[:,0], d[:,1]
#z_file = float(re_z.match(os.path.basename(filename)).groups()[0])
#z = C.pspec.f2z(.151)
z = C.pspec.f2z(freq)    #Average frequency
print 'Redshift = ',z
k3pk = ks**3 / (2*n.pi**2) * pk
p.plot(ks, k3pk * mean_temp(z)**2, 'm-',label='Lidz EoR')
#p.plot(n.array(kpl_pos), 2*n.array(kpl_pos)**3*theo_noise/(2*n.pi**2), 'c--')
p.gca().set_yscale('log', nonposy='clip')
p.xlabel(r'$k\ [h\ {\rm Mpc}^{-1}]$', fontsize='large')
p.ylabel(r'$k^3/2\pi^2\ P(k)\ [{\rm mK}^2]$', fontsize='large')
#p.ylim(1e0,1e7)
p.ylim(1e0,max_y)
p.xlim(0, 0.6)
#p.xlim(0, 3.0)
p.grid()
p.legend(numpoints=1)
p.show()
p.savefig('pspec.png')

sys.exit()

#p.subplot(121)
#p.xlabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$')
#p.ylabel(r'$P(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$')
#p.ylim(-1e2,1e5)
#p.grid()
#p.subplot(122)
#p.xlabel(r'$k\ [h\ {\rm Mpc}^{-1}]$')
#p.ylabel(r'$k^3/2\pi^2\ P(k)\ [{\rm mK}^2]$')
#p.ylim(-1e3,1e3)
#p.xlim(0, 0.6)
#p.grid()
