import numpy as n
import pylab as pl
import capo as C
dirtyHZ = n.load('pspec_FGSub_DirtyWide_95_115_I.npz')
resHZ = n.load('pspec_FGSub_GResidualBHWide_95_115_I.npz')
dirty = n.load('pspec_FGSub_DirtyDelay_95_115_I.npz')
res = n.load('pspec_FGSub_ResidualDelayWide_95_115_I.npz')
reskais = n.load('pspec_FGSub_ResidualDelay_kaiser_95_115_I.npz')
#master = n.load('pspec_8Day-MASTER_95_115_I.npz')
print dirty['freq']

k_d = dirty['k']
k_r = res['k']
k_dhz = dirtyHZ['k']
k_rhz = resHZ['k']
#k_rk = reskais['k']
#k_m = master['kpl']


pk_d = dirty['k3pk']
pk_r = res['k3pk']
pk_dhz = dirtyHZ['k3pk']
pk_rhz = resHZ['k3pk']
#pk_rk = reskais['k3pk']
#pk_m = master['pk']

err_d = dirty['k3err']
err_r = res['k3err']
err_dhz = dirtyHZ['k3err']
err_rhz = resHZ['k3err']
#err_rk = reskais['k3err']
#err_m = master['err']
fig = pl.figure()
ax = fig.add_subplot(211)

tau_h = 100 + 15. #in ns
k_h = C.pspec.dk_deta(C.pspec.f2z(dirty['freq'])) * tau_h

pl.vlines(k_h, -1e7, 1e13, linestyles='--', linewidth=1.5)
pl.vlines(-k_h, -1e7, 1e13, linestyles='--', linewidth=1.5)

ax.set_yscale("log",nonposx='clip')
ax.errorbar(k_d,n.abs(pk_d),yerr=err_d,fmt='ko')
ax.errorbar(k_r,n.abs(pk_r),yerr=err_r,fmt='.')
ax.errorbar(k_dhz,n.abs(pk_dhz),yerr=err_dhz,fmt='^')
#ax.errorbar(k_rk,n.abs(pk_rk),yerr=err_rk,fmt='k^')
ax.set_ylim(n.min(n.append(pk_d,pk_r))+err_r.min(),n.max(n.append(pk_r,pk_d))+err_d.max())
#ax.errorbar(k_m,pk_m,yerr=err_m,fmt='--o')
pl.xlabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$', fontsize='large')
pl.ylabel(r'$P(k)[\ {\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$',fontsize='large')
pl.title('z = '+str(round(1.42/dirty['freq'] - 1,2)))

ax2 = fig.add_subplot(212)

tau_h = 100 + 15. #in ns
k_h = C.pspec.dk_deta(C.pspec.f2z(dirtyHZ['freq'])) * tau_h

pl.vlines(k_h, -1e7, 1e13, linestyles='--', linewidth=1.5)
pl.vlines(-k_h, -1e7, 1e13, linestyles='--', linewidth=1.5)

ax2.set_yscale("log",nonposx='clip')
ax2.errorbar(k_dhz,n.abs(pk_dhz),yerr=err_dhz,fmt='ko')
ax2.errorbar(k_rhz,n.abs(pk_rhz),yerr=err_rhz,fmt='.')
ax2.set_ylim(n.min(n.append(pk_rhz,pk_dhz))+err_rhz.min(),n.max(n.append(pk_dhz,pk_rhz))+err_dhz.max())
#ax.errorbar(k_m,pk_m,yerr=err_m,fmt='--o')
pl.xlabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$', fontsize='large')
pl.ylabel(r'$P(k)[\ {\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$',fontsize='large')
pl.title('z = '+str(round(1.42/dirtyHZ['freq'] - 1,2)))


#pl.semilogy(k_d,pk_d,'.')
pl.show()
print dirty.keys()
