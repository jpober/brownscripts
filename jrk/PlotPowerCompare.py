import numpy as n
import pylab as pl
import capo as C

dirty = n.load('pspec_Analysis_DirtyWideFRF_95_115_I.npz')
resW = n.load('pspec_Analysis_ResidualWideFRF_95_115_I.npz')
res = n.load('old_dirtywidefrf.npz')

k_d = dirty['k']
k_rW = resW['k']
k_r = res['k']

pk_d = dirty['k3pk']
pk_rW = resW['k3pk']
pk_r = res['k3pk']

err_d = dirty['k3err']
err_r = res['k3err']
err_rW = resW['k3err']

fig = pl.figure()
ax = fig.add_subplot(111)

tau_h = 100 + 15. #in ns
k_h = C.pspec.dk_deta(C.pspec.f2z(dirty['freq'])) * tau_h

pl.vlines(k_h, -1e7, 1e13, linestyles='--', linewidth=1.5)
pl.vlines(-k_h, -1e7, 1e13, linestyles='--', linewidth=1.5)

ax.set_yscale("log",nonposx='clip')
ax.errorbar(k_d,n.abs(pk_d),yerr=err_d,fmt='ko')
ax.errorbar(k_rW,n.abs(pk_rW),yerr=err_rW,fmt='b.')
ax.errorbar(k_r,n.abs(pk_r),yerr=err_r,fmt='r.')
ax.set_ylim(n.min(n.append(pk_d,pk_r))+err_r.min(),n.max(n.append(pk_r,pk_d))+err_d.max())
pl.xlabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$', fontsize='large')
pl.ylabel(r'$P(k)[\ {\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$',fontsize='large')
pl.title('z = '+str(round(1.42/dirty['freq'] - 1,2)))

#ax2 = fig.add_subplot(212)

tau_h = 100 + 15. #in ns
k_h = C.pspec.dk_deta(C.pspec.f2z(dirty['freq'])) * tau_h

pl.vlines(k_h, -1e7, 1e13, linestyles='--', linewidth=1.5)
pl.vlines(-k_h, -1e7, 1e13, linestyles='--', linewidth=1.5)

#ax2.set_yscale("log",nonposx='clip')
#ax2.errorbar(k_d,n.abs(pk_d),yerr=err_d,fmt='ko')
#ax2.errorbar(k_r,n.abs(pk_r),yerr=err_r,fmt='.')
#ax2.set_ylim(n.min(n.append(pk_r,pk_d))+err_r.min(),n.max(n.append(pk_d,pk_r))+err_d.max())

#pl.xlabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$', fontsize='large')
#pl.ylabel(r'$P(k)[\ {\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$',fontsize='large')
#pl.title('z = '+str(round(1.42/dirty['freq'] - 1,2)))


#pl.semilogy(k_d,pk_d,'.')
pl.show()
print dirty.keys()
