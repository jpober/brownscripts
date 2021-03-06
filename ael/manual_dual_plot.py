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
o.add_option('--lidz', action='store_true', default=False,
            help='Add the Lidz EoR curve.')
o.add_option('--fhd', action='store_true', default=False,
            help='Add the FHD EoR curve')
o.add_option('--ref_eor', action='store_true', default=False,
            help="Load an EoR pspec curve from the 21cmFAST results.")
o.add_option('--flat', action='store_true', default=False,
            help='Add the FHD Flat EoR curve')

opts,args = o.parse_args(sys.argv[1:])

c = 3e8

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


afreqs=[]

maxpk3 = 0.0
maxpk  = 0.0
max_x  = 0.0
min_y  = 1e-6


f = n.load(args[0])
freq = n.mean(f['freq'])   # Convert to GHz
z = C.pspec.f2z(freq)

if opts.ref_eor:
          path="/gpfs/data/jpober/alanman/BubbleCube/TiledHpxCubes/pspecs_21cmFAST/"
          files = os.listdir(path)
          Z = [ f.split("_")[3] for f in files ]
          Z = n.array([ float(t[1:]) for t in Z])
          ind = n.argmin(n.abs(Z - z))
          filename =  path + files[ind]
          print 'Reading', filename
          eor_label="21cmFAST pspec"
          d = n.array([map(float, L.split()) for L in open(filename).readlines()])
          ks_eor, k3pk_eor = d[:,0]/0.7, d[:,1]*(1./mean_temp(z)**2)*0.5 
          print n.mean(k3pk_eor)
          print 'Redshift = ',z
          u = opts.length*(freq*1e9)/c
          kperp = C.cosmo_units.u2kperp(u,z)
          pk_eor = 1./(ks_eor**3 / (2*n.pi**2)) * k3pk_eor
          kpl_eor = n.sqrt(ks_eor**2 - kperp**2)
          kpl_eor_fold = n.concatenate((-kpl_eor[::-1], kpl_eor))
          pk_eor_fold  = n.concatenate((pk_eor[::-1],pk_eor))

if opts.flat:
          freq = 0.151
          z  = C.pspec.f2z(freq)    #Average frequency
          filename = glob.glob('/users/alanman/capo/ael/FHD_flat_model.dat')[0]
          print 'Reading', filename
          d = n.array([map(float, L.split()) for L in open(filename).readlines()])
          eor_label="Flat Spectrum"
          ks_eor, pk_eor = d[:,0], d[:,1]*(1./mean_temp(z)**2)*(0.7)**3*0.5

if opts.fhd:
          filename = glob.glob('/users/alanman/capo/ael/FHD_EoR_model.dat')[0]
          print 'Reading', filename
          d = n.array([map(float, L.split()) for L in open(filename).readlines()])
          eor_label="Furlanetto EoR"
          ks_eor, pk_eor = d[:,0], d[:,1]*(1./mean_temp(z)**2)*(0.7)**3

if opts.lidz:
          filename = glob.glob('/users/alanman/capo/ael/lidz_mcquinn_k3pk/*8.34.dat')[0]
          print 'Reading', filename
          d = n.array([map(float, L.split()) for L in open(filename).readlines()])
          eor_label = "Lidz EoR"
          ks_eor, pk_eor = d[:,0], d[:,1]

ref_curve = (opts.lidz + opts.fhd + opts.flat + opts.ref_eor) > 0

if ref_curve and not opts.ref_eor:
    print 'Redshift = ',z
    u = opts.length*(freq*1e9)/c
    kperp = C.cosmo_units.u2kperp(u,z)
    k3pk_eor = ks_eor**3 / (2*n.pi**2) * pk_eor
    kpl_eor = n.sqrt(ks_eor**2 - kperp**2)
    kpl_eor_fold = n.concatenate((-kpl_eor[::-1], kpl_eor))
    pk_eor_fold  = n.concatenate((pk_eor[::-1],pk_eor))


n_files = len(args)
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

    print n.mean(k3pk)

    max_x  = n.max([max_x, n.max(kpl)])
    maxpk3 = n.max([maxpk3,n.max(k3pk)])
    maxpk  = n.max([maxpk ,n.max(pk)])

    max_y3 = 0.5e2*maxpk3; max_y = 0.5e2*maxpk

    if opts.right:
         p.errorbar(ks,k3pk,yerr=k3err, color=colors[j], fmt='.',label=os.path.basename(filename))
         p.vlines(k_h, -1e7, max_y, linestyles='--', linewidth=1.5)
         if ref_curve and j == n_files -1 : p.plot(kpl_eor_fold,pk_eor_fold * mean_temp(z)**2, label=eor_label)
         if j == len(args)-1:
           p.xlabel(r'$k\ [h\ {\rm Mpc}^{-1}]$', fontsize='large')
           p.ylabel(r'$k^3/2\pi^2\ P(k)\ [{\rm mK}^2]$', fontsize='large')
         if opts.lin:   p.ylim(-max_y3,max_y3)
         else:          p.ylim(min_y,max_y3)
         p.xlim(0, max_x*1.1)
    elif opts.left:
         p.errorbar(kpl,pk,yerr=err, color=colors[j], fmt='.',label=os.path.basename(filename))
         p.vlines(k_h, -1e7, max_y, linestyles='--', linewidth=1.5)
         p.vlines(-k_h, -1e7, max_y, linestyles='--', linewidth=1.5)
         if j == len(args)-1:
           p.xlabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$', fontsize='large')
           p.ylabel(r'$P(k)[\ {\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$',fontsize='large')
         #if ref_curve and j == n_files - 1:   p.plot(ks_eor, k3pk_eor * mean_temp(z)**2, 'm-',label=eor_label)
         if ref_curve and j == n_files -1 : p.plot(kpl_eor_fold,pk_eor_fold * mean_temp(z)**2, label=eor_label)
         if opts.lin:   p.ylim(-max_y,max_y)
         else:          p.ylim(min_y,max_y)
         p.xlim(-max_x*1.1, max_x*1.1)
    else:
         p.subplot(121)
         p.errorbar(kpl,pk,yerr=err, color=colors[j], fmt='.')
         p.vlines(k_h, -1e7, max_y, linestyles='--', linewidth=1.5)
         p.vlines(-k_h, -1e7, max_y, linestyles='--', linewidth=1.5)
         if ref_curve and j == n_files -1 : p.plot(kpl_eor_fold,pk_eor_fold * mean_temp(z)**2, label=eor_label)
         if j == len(args)-1:
           if opts.lin:   p.ylim(-max_y,max_y)
           else:          p.ylim(min_y,max_y)
           p.xlim(-max_x*1.1, max_x*1.1)
           p.xlabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$', fontsize='large')
           p.ylabel(r'$P(k)[\ {\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$',fontsize='large')
         p.subplot(122)
         p.errorbar(ks,k3pk,yerr=k3err, color=colors[j], fmt='.',label=os.path.basename(filename))
         p.vlines(k_h, -1e7, max_y, linestyles='--', linewidth=1.5)
         if ref_curve and j == n_files - 1:   p.plot(ks_eor, k3pk_eor * mean_temp(z)**2, 'm-',label=eor_label)
         if j == len(args)-1:
           if opts.lin:   p.ylim(-max_y3,max_y3)
           else:          p.ylim(min_y,max_y3)
           p.xlim(0, max_x*1.1)
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
p.legend(numpoints=1, prop={'size':12}, loc=4)


p.show()

