import numpy as n
import pylab as pl
import capo as C
import aipy as a
from glob import glob

folded=True
#dirty = n.load('pspec_PG_Step_95_115_I.npz')
#res = n.load('pspec_PG_StepFS_95_115_I.npz')
#resW = n.load('pspec_PG_StepFS_95_115_I.npz')
files = glob('pspec_PG_Step_*.npz')

files = n.sort(files)
print files
def loadNPZ(NPZfile,fold=True):
    f = n.load(NPZfile)
    if fold:
        pk = n.abs(f['k3pk'])
        err = n.abs(f['k3err'])
        k = f['k']
    else:
        pk = f['pk']
        err = f['err']
        k = f['kpl']
    return pk,err,k

for i in files:
    if i == files[0]:
        pk,err,k = loadNPZ(i,True)
    else:
        pk_i,err_i,k_i = loadNPZ(i,True)
        pk = n.vstack((pk,pk_i))
        err = n.vstack((err,err_i))
        k = n.vstack((k,k_i))

chans=203.
uv = a.miriad.UV('/users/jkerriga/Pzen.2456242.30605.uvcRREcACOTUcHPA')
aa = a.cal.get_aa('psa6240_FHD', uv['sdf'], uv['sfreq'], uv['nchan'])
filters = C.dspec.wedge_width_by_bl(aa, uv['sdf'], uv['nchan'], offset=15.0)
bins = filters[(41,49)]

fig = pl.figure()
ax = fig.add_subplot(111)

tau_h = 100 + 15. #in ns
k_h = C.pspec.dk_deta(C.pspec.f2z(0.16)) * tau_h

pl.vlines(k_h, -1e7, 1e13, linestyles='--', linewidth=1.5)
#pl.vlines(-k_h, -1e7, 1e13, linestyles='--', linewidth=1.5)
pl.axvspan(-k_h,k_h,color='red',alpha=0.5)

ax.set_yscale("log",nonposx='noclip')
mA = n.argmax(n.sum(pk,1))
for j in range(len(files)):
    print pk[j,:]
    ax.errorbar(k[j,:],pk[j,:]/pk[0,:],yerr=err[j,:]/pk[0,:],fmt='o')
#ax.errorbar(k_d,n.abs(pk_d)/pk_d,yerr=err_d/pk_d,fmt='ko')
#ax.errorbar(k_rW,n.abs(pk_rW)/pk_d,yerr=err_rW/pk_d,fmt='r.')
#ax.errorbar(k_r,n.abs(pk_r)/pk_d,yerr=err_r/pk_d,fmt='b.')
#ax.set_ylim(n.min(n.append(pk_d,pk_r))+err_r.min(),n.max(n.append(pk_r,pk_d))+err_d.max())
ax.set_ylim(1e-8,1)
ax.set_xticklabels([])

pl.xlabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$', fontsize='large')
pl.ylabel(r'$P(k)[\ {\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$',fontsize='large')
#pl.title('z = '+str(round(1.42/dirty['freq'] - 1,2)))

#tau_h = 100 + 15. #in ns
#k_h = C.pspec.dk_deta(C.pspec.f2z(dirty['freq'])) * tau_h

pl.vlines(k_h, -1e7, 1e13, linestyles='--', linewidth=1.5)
pl.vlines(-k_h, -1e7, 1e13, linestyles='--', linewidth=1.5)

#ax2 = fig.add_subplot(212)
#err = (pk_r/pk_d)*n.sqrt((err_r/pk_r)**2 +(err_d/pk_d)**2)#/n.sqrt(100)
#ax2.errorbar(k_r,pk_rW/pk_d,yerr=err,fmt='-+')
#pl.hlines(1,-0.5,0.5,linestyles='--')
#ax2.set_xticklabels(k_d)
#ax2.set_ylim(0,1.1)
#pl.vlines(k_h, 0, 1.1, linestyles='--', linewidth=1.5)
#pl.vlines(-k_h, 0, 1.1, linestyles='--', linewidth=1.5)
#pl.axvspan(-k_h,k_h,color='red',alpha=0.5)
#ax2.set_yscale("log",nonposx='clip')
#ax2.errorbar(k_d,n.abs(pk_d),yerr=err_d,fmt='ko')
#ax2.errorbar(k_r,n.abs(pk_r),yerr=err_r,fmt='.')
#ax2.set_ylim(n.min(n.append(pk_r,pk_d))+err_r.min(),n.max(n.append(pk_d,pk_r))+err_d.max())

#pl.xlabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$', fontsize='large')
#pl.ylabel(r'$P(k)[\ {\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$',fontsize='large')
#pl.title('z = '+str(round(1.42/dirty['freq'] - 1,2)))


#pl.semilogy(k_d,pk_d,'.')
pl.show()
#print dirty.keys()
