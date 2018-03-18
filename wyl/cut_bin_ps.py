import numpy as np, sys
from scipy.io.idl import readsav
from plot_vis import *
#power = sys.argv[1]
#power_diff = sys.argv[2]
def ps_cut(power_2d,power_diff,bin_file=None,thresh=1e7,kperp_min_cut=0.003,kperp_max_cut=0.05,kpara_min_cut=0.11,kpara_max_cut=1.08,coarse_band_extent=2):
    d = readsav(power_2d)
    diff = readsav(power_diff)
    kperp = d['kperp_edges']
    kx = diff['kx_mpc']
    ky = diff['ky_mpc']
    kz = diff['kz_mpc']
    h = diff['hubble_param']
    k = np.zeros((kx.size,kx.size))
    for ii in range(kx.size):
        for jj in range(ky.size):
            k[jj][ii]=np.sqrt(kx[ii]**2+ky[jj]**2)
    k=k.flatten()
#dat = organize_power(d,'power_3d')
#dif = organize_power(diff,'power_diff')
    if bin_file is None:
        print "find pixels using model"
        q = np.logical_and(d['power']<thresh,d['power']>0)
    else: 
        fbin = readsav(bin_file)
        q = np.array(fbin['bin_1to2d']).astype(bool)
    r = np.zeros((kz.size,kx.size*ky.size),dtype=bool)
    n = np.digitize(k,kperp)
    q[:,np.where(kperp[:-1] > kperp_max_cut)[0]] = False
    q[:,np.where(kperp[1:] < kperp_min_cut)[0]] = False
    q[np.where(kz > kpara_max_cut)[0]] = False
    q[np.where(kz < kpara_min_cut)[0]] = False
    q[:,12] = False
    q[:,13] = False
    q[:,21] = False
    for ii in range(coarse_band_extent):
        q[24+ii::24] = False
        q[24-ii::24] = False
    for nn in range(kperp.size-1):  
        if kperp[nn] > kperp_max_cut: continue
        if kperp[nn] < kperp_min_cut: continue          
        for zi in range(kz.size): 
            if kz[zi] < kpara_min_cut: continue
            if kz[zi] > kpara_max_cut: continue
            if q[zi][nn]:             
                ind0 = np.where(n==nn+1)
                r[zi][ind0[0]] = True
#    return r
    dif = diff['power_diff']*(h**3)
    wgt = diff['weight_diff']
    #wgt /= np.sum(wgt)
    r = r.reshape(dif.shape)   
    ind = np.where(r)
    cut = dif[ind]#*(wgt[ind]>10.*np.mean(wgt))
    #n0ind = np.nonzero(cut)
    kzpower = np.divide(np.sum(dif*wgt*r,axis=(1,2)),np.sum(wgt*r,axis=(1,2)))
    return q,d['power']*(h**3), cut, wgt[ind], kzpower, kz/h

def check_2to3d(power_2d,power_diff,bin_file):
    d = readsav(power_2d)
    diff = readsav(power_diff)
    fbin = readsav(bin_file)
    q = np.array(fbin['bin_1to2d']).astype(bool)
    r0 = np.array(fbin['bin_arr_3d']).astype(bool)
    kperp = d['kperp_edges']
    kx = diff['kx_mpc']
    ky = diff['ky_mpc']
    kz = diff['kz_mpc']
    k = np.zeros((kx.size,kx.size))
    for ii in range(kx.size):
        for jj in range(ky.size):
            k[jj][ii]=np.sqrt(kx[ii]**2+ky[jj]**2)
    k=k.flatten()
    r = np.zeros((kz.size,kx.size*ky.size),dtype=bool)
    n = np.digitize(k,kperp)
    for nn in range(kperp.size-1):
        for zi in range(kz.size):
            if q[zi][nn]:
                ind0 = np.where(n==nn+1)
                r[zi][ind0[0]] = True
    return r0,r


def histbin(data,wgts,bins):
    w = wgts/np.sum(wgts)
    mean = np.sum(data*w)
    d = (data-mean)*w/np.mean(w) + mean             
    plt.hist(data,bins=bins)     
    plt.yscale("log",linthresh=1e-5)
    plt.xlabel('diff power')      
    plt.ylabel('counts')            
    plt.grid(True)              
    plt.show() 
#bins = np.linspace(-1e8,1e8,1000)
#histbin(cut,bins)

def statdata(cut,wgt):
    w = wgt/np.sum(wgt)
    cutmean = np.sum(w*cut)
    cutmeanvar = np.sum(((cut-cutmean)*w)**2)
    cutmeanstd = np.sqrt(cutmeanvar)
    return cutmean, cutmeanstd, cutmean/cutmeanstd

def plot_ps(p,q,clim=None):
    kper_min = 0.0032069491162178276
    kper_max = 1.2186406641627747
    kpar_min = 0.00537899
    kpar_max = 0.30424922974941621
    fig = plt.figure()
    p1 = fig.add_subplot(1,1,1)
    i1=p1.imshow((p-q*1e10)[::-1],aspect='equal',cmap='jet',extent=(kpar_min,kpar_max,kper_min,kper_max),clim=clim,norm=LogNorm())
    p1.set_xscale('log',linthresh=1e-3)
    p1.set_yscale('log',linthresh=1e-3)
    plt.colorbar(i1)
    #p2 = fig.add_subplot(1,2,2)
    #i2=p2.imshow(q[::-1][8:],aspect='equal',cmap='jet',extent=(kpar_min,kpar_max,kper_min,kper_max))
    #p2.set_xscale('log') 
    #p2.set_yscale('log') 
    #plt.colorbar(i2)
    plt.show()

def kzpower_plot(k1,k2,k3,kz=None):
    if kz is None: kz = np.linspace(0,1.22,k1.size)
    kzres = kz[1] - kz[0]
    kz_bin,k1_bin,k2_bin,k3_bin = [],[],[],[]
    for ii in range(kz.size):
        kz_bin.append(kz[ii]-kzres/2)
        k1_bin.append(k1[ii])
        k2_bin.append(k2[ii])
        k3_bin.append(k3[ii])
        kz_bin.append(kz[ii]+kzres/2)
        k1_bin.append(k1[ii])
        k2_bin.append(k2[ii])
        k3_bin.append(k3[ii])
    kz_bin = np.array(kz_bin)
    k1_bin = np.array(k1_bin)/1e6
    k2_bin = np.array(k2_bin)/1e6
    k3_bin = np.array(k3_bin)/1e6
    #ind = np.where(kz_bin>0)
    p1,=plt.plot(kz_bin,k1_bin,label='FHD - OF')
    p2,=plt.plot(kz_bin,k2_bin,label='FHD - FO')
    p3,=plt.plot(kz_bin,k3_bin,label='FO - OF')
    plt.xlabel('$k_\parallel (hMpc^{-1})$')
    plt.ylabel('$P_k (10^6 mK^2h^{-3}Mpc^3)$')
    plt.xscale('log', linthresh=1e-3)
    plt.xlim(1e-1,np.max(kz_bin))
    plt.ylim(-0.6,2.2)
    ind0 = np.argmin(np.abs(kz_bin-0.7-kzres))
    ind1 = np.argmin(np.abs(kz_bin-0.7+2*kzres))
    plt.axvline(x=kz_bin[ind0],color='black',linewidth=0.5,linestyle='-.')
    plt.axvline(x=kz_bin[ind1],color='black',linewidth=0.5,linestyle='-.')
    plt.fill_between(kz_bin,-0.7,2.3,where=np.isnan(k1_bin),color='silver',alpha=0.8)
    below_wedge_ind = np.where(np.isnan(k1_bin)==False)[0][0]
    wedge_cut = np.zeros(kz_bin.shape,dtype=bool)
    wedge_cut[:below_wedge_ind] = True
    plt.fill_between(kz_bin,-0.7,2.3,where=wedge_cut,color='grey',alpha=0.8)
    
#    if log_scale: plt.yscale('symlog',linthresh=1e-7)
  #  plt.ticklabel_format(axis='y', style='sci',scilimits=(-2,7))
    plt.grid(True,axis='y')
    plt.legend(handles=[p1,p2,p3],loc=1)
    plt.show()
    
