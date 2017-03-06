import numpy as n
import matplotlib
matplotlib.use('Agg')
import pylab as pl
import capo as C
import aipy as a
from glob import glob


folded=True
dir = 'PG_Ratio'
npz = glob('./'+dir+'/*/I/*.npz')
print npz
chanlist = []
pspecs = {}
for i in npz:
    x = n.load(i)
    chans = i.split('/')[2]
    pspecs[chans+'_pk'] = x['pk']
    pspecs[chans+'_kpl'] = x['kpl']
    pspecs[chans+'_err'] = x['err']
    pspecs[chans+'_freq'] = x['freq']
    chanlist.append(chans)

#print pspecs['50_70_pk']

chans=203.
uv = a.miriad.UV('../../Pzen.2456242.30605.uvcRREcACOTUcHPA')
aa = a.cal.get_aa('psa6240_FHD', uv['sdf'], uv['sfreq'], uv['nchan'])
filters = C.dspec.wedge_width_by_bl(aa, uv['sdf'], uv['nchan'], offset=15.0)
bins = filters[(41,49)]

fig = pl.figure()
ax = fig.add_subplot(111)

tau_h = 100 + 15. #in ns
k_h = C.pspec.dk_deta(C.pspec.f2z(pspecs[chanlist[0]+'_freq'])) * tau_h

for c in chanlist:
    pspecs[c+'_pkfold'] = [pspecs[c+'_pk'][9]]
    pspecs[c+'_errfold'] = [pspecs[c+'_err'][9]]
    for i in range(10):

        pspecs[c+'_pkfold'].append(n.mean((pspecs[c+'_pk'][10-i],pspecs[c+'_pk'][10+i])))
        pspecs[c+'_errfold'].append(n.mean((pspecs[c+'_err'][10-i],pspecs[c+'_err'][10+i])))
        
#pl.vlines(k_h, -1e7, 1e13, linestyles='--', linewidth=1.5)
#pl.vlines(-k_h, -1e7, 1e13, linestyles='--', linewidth=1.5)
pl.axvspan(-k_h,k_h,color='red',alpha=0.5)
ax.set_yscale("linear")

for c in chanlist:
    print c
    if folded==True:
        ax.errorbar(pspecs[c+'_kpl'][10:21],pspecs[c+'_pkfold'],yerr=pspecs[c+'_errfold'],fmt='o',label='z = '+str(round(1.42/pspecs[c+'_freq'] - 1,1)))
        ax.set_xlim(-0.05,n.max(pspecs[chanlist[0]+'_kpl'])+0.05)
    else:
        ax.errorbar(pspecs[c+'_kpl'],pspecs[c+'_pk'],yerr=pspecs[c+'_err'],fmt='o',label='z = '+str(round(1.42/pspecs[c+'_freq'] - 1,1)))
        ax.set_xlim(-1*(n.max(pspecs[chanlist[0]+'_kpl'])+0.05),n.max(pspecs[chanlist[0]+'_kpl'])+0.05)

#ax.set_xlim(-0.05,n.max(pspecs[chanlist[0]+'_kpl'])+0.05)
ax.legend(bbox_to_anchor=(1.1,0.4))#loc=4)
pl.hlines(1,-0.6,0.6,linestyles='--')
#ax.errorbar(k_rW,n.abs(pk_rW),yerr=err_rW,fmt='r.')
#ax.errorbar(k_r,n.abs(pk_r),yerr=err_r,fmt='b.')
#ax.set_ylim(n.min(n.append(pk_d,pk_r))+err_r.min(),n.max(n.append(pk_r,pk_d))+err_d.max())
#ax.set_xlim(-0.05,n.max(k_rW)+0.05)
#ax.set_xticklabels([])

pl.xlabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$', fontsize='large')
pl.ylabel(r'$P_{R}(k)/P_{D}(k)$',fontsize='large')
#pl.title('z = '+str(round(1.42/dirty['freq'] - 1,2)))

#tau_h = 100 + 15. #in ns
#k_h = C.pspec.dk_deta(C.pspec.f2z(dirty['freq'])) * tau_h

#pl.vlines(k_h, -1e7, 1e13, linestyles='--', linewidth=1.5)
#pl.vlines(-k_h, -1e7, 1e13, linestyles='--', linewidth=1.5)
#rho = (n.mean(pk_r*pk_d)-pk_r*pk_d)/(err_d*err_r)#n.corrcoef(pk_d,pk_r)
#print rho
#ax2 = fig.add_subplot(212)


#ax2.errorbar(k_r,pk_r,yerr=err_r,fmt='+')
#ax2.errorbar(k_d,n.abs(pk_d),yerr=err_d,fmt='o')
#ax2.errorbar(k_rW,pk_rW,yerr=err_rW,fmt='xr')
#ax2.set_ylim(n.min(pk_d)+n.max(err_d),n.max(pk_d)+n.max(err_d))
#pl.hlines(1,-0.5,0.5,linestyles='--')
#ax2.set_xticklabels(k_d)

#pl.vlines(k_h, 0, 1.1, linestyles='--', linewidth=1.5)
#pl.vlines(-k_h, 0, 1.1, linestyles='--', linewidth=1.5)
#pl.axvspan(-k_h,k_h,color='red',alpha=0.5)



#pl.xlabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$', fontsize='large')
#pl.ylabel(r'$P(k)[\ {\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$',fontsize='large')
#pl.title('z = '+str(round(1.42/dirty['freq'] - 1,2)))


#pl.semilogy(k_d,pk_d,'.')
pl.savefig('RatioPSPEC.png')
pl.show()

