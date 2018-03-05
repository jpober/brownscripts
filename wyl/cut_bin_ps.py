import numpy as np, sys
from scipy.io.idl import readsav
from plot_vis import *
#power = sys.argv[1]
#power_diff = sys.argv[2]
def ps_cut(power_2d,power_diff,thresh=1e7,kperp_min_cut=0.003,kperp_max_cut=0.05,kpara_min_cut=0.11,kpara_max_cut=1.08):
    d = readsav(power_2d)
    diff = readsav(power_diff)
    kperp = d['kperp_edges']
    kx = diff['kx_mpc']
    ky = diff['ky_mpc']
    kz = diff['kz_mpc']
    k = np.zeros((kx.size,kx.size))
    for ii in range(kx.size):
        for jj in range(ky.size):
            k[jj][ii]=np.sqrt(kx[ii]**2+ky[jj]**2)
    k=k.flatten()
#dat = organize_power(d,'power_3d')
#dif = organize_power(diff,'power_diff')
    q = np.logical_and(d['power']<thresh,d['power']>0)
    r = np.zeros((kz.size,kx.size*ky.size),dtype=bool)
    n = np.digitize(k,kperp)
    q[:,np.where(kperp[:-1] > kperp_max_cut)[0]] = False
    q[:,np.where(kperp[1:] < kperp_min_cut)[0]] = False
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
    dif = diff['power_diff']
    wgt = diff['weight_diff']
    #wgt /= np.sum(wgt)
    r = r.reshape(dif.shape)   
    ind = np.where(r)
    cut = dif[ind]#*(wgt[ind]>10.*np.mean(wgt))
    #n0ind = np.nonzero(cut)
    return q,d['power'], cut, wgt[ind]

def histbin(data,bins):              
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
    fig = plt.figure()
    p1 = fig.add_subplot(1,2,1)
    i1=p1.imshow(p[::-1],aspect='equal',cmap='jet',extent=(0.005,0.30424922974941621,0.003,1.22),clim=clim,norm=LogNorm())
    p1.set_xscale('log',linthresh=1e-3)
    p1.set_yscale('log',linthresh=1e-3)
    plt.colorbar(i1)
    p2 = fig.add_subplot(1,2,2)
    i2=p2.imshow(q[::-1],aspect='equal',cmap='jet',extent=(0.005,0.30424922974941621,0.003,1.22))
    p2.set_xscale('log',linthresh=1e-3) 
    p2.set_yscale('log',linthresh=1e-3) 
    plt.colorbar(i2)
    plt.show()

