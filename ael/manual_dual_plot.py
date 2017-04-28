#! /usr/bin/env python
import matplotlib
#matplotlib.use('Agg')
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
o.add_option('--length', default=14.,
            help='Baseline length in meters.')
o.add_option('-o','--ofilename',
            help='Output file name')
o.add_option('-l','--left', action='store_true',
            help='Show only left plot')
o.add_option('-r', '--right', action='store_true',
            help='Show only right plot')
o.add_option('--lin', action='store_true',
            help='Plot on a linear scale')


opts,args = o.parse_args(sys.argv[1:])

def mean_temp(z):
    return 28. * ((1.+z)/10.)**.5 # mK

if opts.ofilename is None:
	#opts.ofilename = ".".join(args[0].split('.')[:-1])    #Drop npz extension, replace with png
	#opts.ofilename = opts.ofilename+".png"
        opts.ofilename = 'pspec.png'

#fig=p.figure(figsize=(12,7.2))
fig=p.figure()

#colors = 'bgrcmyk' * 10
colors = ['black', 'grey', 'darkgreen', 'lime', 'mediumblue', 'darkviolet', 'maroon', 'red', 'rosybrown' ] 
 

#tau_h = 100 + 15. #in ns
#k_h = C.pspec.dk_deta(C.pspec.f2z(.151)) * tau_h

tau_h = float(opts.length)*(100.)/(a.const.len_ns)
k_h = C.pspec.dk_deta(C.pspec.f2z(.151)) * tau_h

p.subplot(111)

#theo_noise = noise_level(freq=freq)

afreqs=[]

maxpk3 = 0.0
maxpk  = 0.0
min_y  = 1e-3


for j,filename in enumerate(args):
    print 'Reading', filename
    f = n.load(filename)
    afreqs.append(f['freq'])

    try:
        ks,k3pk,k3err = f['k'], f['k3pk'], f['k3err']
        kpl,pk,err = f['kpl'], f['pk'], f['err']
    except:
        try:
            ks,k3pk,k3err = f['k'], f['pIv_fold'], f['pIv_fold_err']
            kpl,pk,err = f['kpl'], f['pIv'], f['pIv_err']
            k3pk *= (ks**3/(2*n.pi))      #New version of capo doesn't save k3pk with the k3 factor, so apply here.
            k3err *= (ks**3/(2*n.pi))
        except:
            raise AttributeError


    maxpk3 = n.max([maxpk3,n.max(k3pk)])
    maxpk  = n.max([maxpk ,n.max(pk)])

    max_y3 = 0.5e2*maxpk3; max_y = 0.5e2*maxpk

    if opts.right:
         p.errorbar(ks,k3pk,yerr=k3err, color=colors[j], fmt='.',label=os.path.basename(filename))
         p.vlines(k_h, -1e7, max_y, linestyles='--', linewidth=1.5)
         if j == len(args)-1:
           p.xlabel(r'$k\ [h\ {\rm Mpc}^{-1}]$', fontsize='large')
           p.ylabel(r'$k^3/2\pi^2\ P(k)\ [{\rm mK}^2]$', fontsize='large')
         if opts.lin:   p.ylim(-max_y3,max_y3)
         else:          p.ylim(min_y,max_y3)
         p.xlim(0, 0.6)
    elif opts.left:
         p.errorbar(kpl,pk,yerr=err, color=colors[j], fmt='.',label=os.path.basename(filename))
         p.vlines(k_h, -1e7, max_y, linestyles='--', linewidth=1.5)
         p.vlines(-k_h, -1e7, max_y, linestyles='--', linewidth=1.5)
         if j == len(args)-1:
           p.xlabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$', fontsize='large')
           p.ylabel(r'$P(k)[\ {\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$',fontsize='large')
         if opts.lin:   p.ylim(-max_y,max_y)
         else:          p.ylim(min_y,max_y)
         p.xlim(-0.6, 0.6)
    else:
         p.subplot(121)
         p.errorbar(kpl,pk,yerr=err, color=colors[j], fmt='.')
         p.vlines(k_h, -1e7, max_y, linestyles='--', linewidth=1.5)
         p.vlines(-k_h, -1e7, max_y, linestyles='--', linewidth=1.5)
         if j == len(args)-1:
           if opts.lin:   p.ylim(-max_y,max_y)
           else:          p.ylim(min_y,max_y)
           p.xlim(-0.6, 0.6)
           p.xlabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$', fontsize='large')
           p.ylabel(r'$P(k)[\ {\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$',fontsize='large')
         p.subplot(122)
         p.errorbar(ks,k3pk,yerr=k3err, color=colors[j], fmt='.',label=os.path.basename(filename))
         p.vlines(k_h, -1e7, max_y, linestyles='--', linewidth=1.5)
         if j == len(args)-1:
           if opts.lin:   p.ylim(-max_y3,max_y3)
           else:          p.ylim(min_y,max_y3)
           p.xlim(0, 0.6)
           p.xlabel(r'$k\ [h\ {\rm Mpc}^{-1}]$', fontsize='large')
           p.ylabel(r'$k^3/2\pi^2\ P(k)\ [{\rm mK}^2]$', fontsize='large')


if opts.right or opts.left:
       if not opts.lin: p.gca().set_yscale('log', nonposy='clip')
       p.grid()
else:
       p.subplot(121)
       p.grid()
       if not opts.lin: p.gca().set_yscale('log', nonposy='clip')
       p.subplot(122)
       p.grid()
       if not opts.lin: p.gca().set_yscale('log', nonposy='clip')

freq=n.average(afreqs)

#TODO --- Enable this.  Plot theoretical EoR curve
#fiilename = glob.glob('/users/alanman/capo/ael/lidz_mcquinn_k3pk/*7.32.dat')[0]
#print 'Reading', filename
#d = n.array([map(float, L.split()) for L in open(filename).readlines()])
#ks, pk = d[:,0], d[:,1]
z = C.pspec.f2z(freq)    #Average frequency
print 'Redshift = ',z
#k3pk = ks**3 / (2*n.pi**2) * pk
#
#p.plot(ks, k3pk * mean_temp(z)**2, 'm-',label='Lidz EoR')
##p.plot(n.array(kpl_pos), 2*n.array(kpl_pos)**3*theo_noise/(2*n.pi**2), 'c--')
#p.gca().set_yscale('log', nonposy='clip')
##p.ylim(1e0,1e7)

p.legend(numpoints=1, prop={'size':12}, loc=4)
p.show()
#p.savefig(opts.ofilename)

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

