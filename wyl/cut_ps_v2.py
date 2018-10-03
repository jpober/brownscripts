import numpy as np, sys, glob
from scipy.io.idl import readsav
from plot_vis import *
import scipy.special as sp

def erf(x): return sp.erf(x)

def erfinv(y): return sp.erfinv(y)

def rebin(bin_file,kperp_min_cut=0.0,kperp_max_cut=0.1,kpara_min_cut=0.12,kpara_max_cut=1.2,coarse_band_extent=4,cut_vmodes=True):
    fbin = readsav(bin_file)
    q = np.array(fbin['bin_1to2d']).astype(bool)
    kx = fbin['kx_mpc']
    ky = fbin['ky_mpc']
    kz = fbin['kz_mpc']
    kperp = fbin['kperp_edges']
    h = fbin['hubble_param']
    k = np.zeros((ky.size,kx.size))
    for ii in range(kx.size):
        for jj in range(ky.size):
            k[jj][ii]=np.sqrt(kx[ii]**2+ky[jj]**2)
    k=k.flatten()
    r = np.zeros((kz.size,kx.size*ky.size),dtype=bool)
    n = np.digitize(k,kperp)
    q[:,np.where(kperp[:-1] > kperp_max_cut)[0]] = False
    q[:,np.where(kperp[1:] < kperp_min_cut)[0]] = False
    q[np.where(kz > kpara_max_cut)[0]] = False
    q[np.where(kz < kpara_min_cut)[0]] = False
    if cut_vmodes:
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
    r = r.reshape(kz.size,ky.size,kx.size)
    bin_dict = {"2d_mask":q,"3d_mask":r,"kx":kx/h,"ky":ky/h,"kz":kz/h,"h":h,"k":k/h,"kperp":kperp/h}
    return bin_dict

def get_diff_power(dif_fn, bin_dict):
    diff = readsav(dif_fn)
    h = bin_dict["h"]
    dm = bin_dict["3d_mask"]
    dp = diff['power_diff']*(h**3)
    dw = diff['weight_diff']
    cut = dp[dm]
    wgt = dw[dm]
    p_kpar = np.divide(np.sum(dp*dw*dm,axis=(1,2)),np.sum(dw*dm,axis=(1,2)))
    pkpt = np.sum(dp*dw*dm,axis=0).flatten()
    wkpt = np.sum(dw*dm,axis=0).flatten()
    n = np.digitize(bin_dict["k"],bin_dict["kperp"])
    p_kperp = []
    for ii in range(1,bin_dict["kperp"].size):
        ind = np.where(n==ii)
        p_kperp.append(np.divide(np.sum(pkpt[ind]),np.sum(wkpt[ind])))
    p_kperp = np.array(p_kperp)
    power_dict = {"p_kpar": p_kpar, "p_kperp":p_kperp, "power": cut, "weight": wgt}
    return power_dict

def bin_3d(kx,ky,kz):
    k = np.zeros(kz.size,ky.size,kx.size)
    for ix in range(kx.size):
        for iy in range(ky.size):
            for iz in range(kz.size):
                k[ix][iy][iz] = np.sqrt(kx[ix]**2+ky[iy]**2+kz[iz]**2)
    return k

def get_3d_ps(powerfn, bin_dict, slots_1d):
    d = readsav(powerfn)
    p = np.copy(d['power_3d'])
    w = np.copy(d['weights_3d'])
    kx = d['kx_mpc']
    ky = d['ky_mpc']
    kz = d['kz_mpc']
    p *= bin_dict['3d_mask']
    w *= bin_dict['3d_mask']
    p1d = np.zeros(slots_1d.size-1)
    w1d = np.zeros(slots_1d.size-1)
    for ix in range(kx.size):
        for iy in range(ky.size):
            for iz in range(kz.size):
                k = np.sqrt(kx[ix]**2+ky[iy]**2+kz[iz]**2)
                ind = np.where(slots_1d<k)[0][-1]
                try:
                    p1d[ind] += p[iz,iy,ix]*w[iz,iy,ix]
                    w1d[ind] += w[iz,iy,ix]*(p[iz,iy,ix]!=0)
                except: pass
    ind = np.where(w1d!=0)[0]
    p1d[ind] /= w1d[ind]
    p1d *= 2
    return p1d

def histbin(data,wgts,bins):
    w = wgts/np.sum(wgts)
    counts,edges,_=plt.hist(data,weights=w,bins=bins)     
    plt.yscale("log",linthresh=1e-5)
    plt.xlabel('diff power')      
    plt.ylabel('counts')            
    plt.grid(True)             
    plt.show() 
    return counts,edges

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
    plt.show()

def binx(x,y):
    x_res = x[1] - x[0]
    xbin = np.linspace(x[0]-x_res/2,x[-1]+x_res/2,2*x.size+1)
    xbool = np.zeros(xbin.shape,dtype=bool)
    ind = np.where(np.isnan(y))
    xbool[1::2][ind] = True
    xbool[0::2][ind] = True
    xbool[2::2][ind] = True
    return xbin, xbool

def xtostep(x):
    x_res = (x[1]-x[0])/2
    x_bin = np.zeros((2*x.size))
    x_bin[0::2] = x-x_res
    x_bin[1::2] = x+x_res
    return x_bin

def ytostep(y):
    y_bin = np.zeros((2*y.size))
    y_bin[0::2] = y
    y_bin[1::2] = y
    return y_bin


def kpar_plot(k1,k2,k3,kz):
    p1,=plt.step(kz,k1/1e6,where='mid',label='FHD - OF')
    p2,=plt.step(kz,k2/1e6,where='mid',label='FHD - FO')
    p3,=plt.step(kz,k3/1e6,where='mid',label='FO - OF')
    plt.xlabel('$k_\parallel (hMpc^{-1})$')
    plt.ylabel('$P_k (10^6 mK^2h^{-3}Mpc^3)$')
    plt.xscale('log', linthresh=1e-3)
    plt.xlim(1e-1,np.max(kz))
    plt.ylim(-0.6,2.2)
    plt.axvline(x=0.70,color='black',linewidth=0.8,linestyle='-.')
    plt.axvline(x=1.1,color='black',linewidth=0.8,linestyle='--')
    kbin, kbool = binx(kz,k1)
    plt.fill_between(kbin,-0.7,2.3,where=kbool,color='silver',alpha=0.8)
    below_wedge_ind = np.where(kbool==False)[0][0]
    wedge_cut = np.zeros(kbin.shape,dtype=bool)
    wedge_cut[:below_wedge_ind] = True
    plt.fill_between(kbin,-0.7,2.3,where=wedge_cut,color='grey',alpha=0.8)
    plt.grid(True,axis='y')
    plt.legend(handles=[p1,p2,p3],loc=1)
    plt.show()

def kperp_plot(k1,k2,k3,kperp):
    kp = 0.5*(kperp[1:]+kperp[:-1])
    p1,=plt.step(kp,k1/1e6,where='mid',label='FHD - OF')
    p2,=plt.step(kp,k2/1e6,where='mid',label='FHD - FO')
    p3,=plt.step(kp,k3/1e6,where='mid',label='FO - OF')
    plt.xlabel('$k_\perp (hMpc^{-1})$')
    plt.ylabel('$P_k (10^6 mK^2h^{-3}Mpc^3)$')
    plt.xlim(0.01,0.06)
    plt.ylim(-1.2,1)
    plt.xscale('log', linthresh=1e-3)
    kbin,kbool = binx(kp,k1)
    plt.fill_between(kbin,-1.2,1,where=kbool,color='silver',alpha=0.8)
    plt.grid(True,axis='y')
    plt.legend(handles=[p1,p2,p3],loc=1)
    plt.show()

def power1dplot(fn):
    lgd_hand = []
    fsp = fn.split('/')[-1].split('_') 
    tp = fsp[9]
    pol = fsp[10]
    d = readsav(fn)
    h = d['hubble_param']
    p = d['power']*(h**3)
    n = d['noise']*(h**3)
    s = (h**3)/np.sqrt(d['weights'])
    k_edges = d['k_edges']/h
    k = k_edges[1:]/2+k_edges[:-1]/2
    p = p*k**3/2/np.pi**2
    n = n*k**3/2/np.pi**2
    s = s*k**3/2/np.pi**2
    plt.step(k,p,where='mid',label='measured power',c='black')
    plt.step(k,s,where='mid',label='1 sigma thermal noise',linestyle='--',c='black')
    k_bin = xtostep(k)
    p_bin = ytostep(p)
    s_bin = ytostep(s)
    pkup = np.sqrt(2)*s_bin*erfinv(0.977-(1-0.977)*erf(p_bin/s_bin/np.sqrt(2))) + p_bin
    #ind0 = np.where(pkup<0)
    #pkup[ind0] = 2*s_bin[ind0]
    #pup = np.sqrt(2)*s_bin*scipy.special.erfinv(0.977-(1-0.977)*scipy.special.erf(p_bin/s_bin/np.sqrt(2)))
    plt.fill_between(k_bin,p_bin-2*s_bin,pkup,color='silver',alpha=0.8)
    plt.plot(k_bin,pkup,label='2 sigma upper limit',c='indigo')
    plt.ylabel('$\Delta^2(mK^2)$')
    plt.xlabel('$k(hMpc^{-1})$')
    plt.grid(True,axis='y')
    plt.legend(loc=2)
    plt.xlim(1e-1,np.max(k))
    plt.ylim(1e2,1e7)
    plt.xscale('log')
    plt.yscale('log')
    plt.show()

def get_1d_limit(fn):
    d = readsav(fn)
    h = d['hubble_param']
    p = d['power']*(h**3)
    s = (h**3)/np.sqrt(d['weights'])
    k_edges = d['k_edges']/h
    k = k_edges[1:]/2+k_edges[:-1]/2
    p = p*k**3/2/np.pi**2
    s = s*k**3/2/np.pi**2
    pkup = np.sqrt(2)*s*erfinv(0.977-(1-0.977)*erf(p/s/np.sqrt(2))) + p
    return k, pkup

def exratio(fi):
    #fi=obs+'_cubeXX__even_odd_joint_bh_res_'+pol+'_averemove_swbh_dencorr_no_horizon_wedge_kperplambda5-50_1dkpower.idlsave'
    d=readsav(fi)
    p=d['power']
    n=d['noise']
    k=d['k_edges']
    k=0.5*(k[1:]+k[:-1])
    h=d['hubble_param']
    k /= h
    sec=np.logical_and(k>0.2,k<0.3)
    return np.mean(p[sec])/np.mean(n[sec])*np.sqrt(np.sum(sec)-1)

def flagr(s):
    sc=np.ma.masked_array(s,np.zeros((len(s)),dtype=bool))
    for ii in range(500):
        m=np.mean(sc)
        q=np.std(sc)
        ind=np.argmax(sc)
        if sc[ind]<m+3*q: break
        else: sc.mask[ind]=True
    print str(ii)+" observations flagged."
    return sc

def metric_wrap(fpath='/users/wl42/data/wl42/FHD_out/fhd_int_PhaseII/ps/data/1d_binning/'):
    fx = glob.glob(fpath+'116*xx*_averemove_swbh_dencorr_no_horizon_wedge_kperplambda5-50_1dkpower.idlsave')
    fy = glob.glob(fpath+'116*yy*_averemove_swbh_dencorr_no_horizon_wedge_kperplambda5-50_1dkpower.idlsave')
    fx.sort()
    fy.sort()
    obs, rx, ry = [], [], []
    for fi in fx: obs.append(int(fi.split('/')[-1].split('_')[0]))
    obs = np.array(obs)
    day = np.int32(obs/86164.1)
    sid = obs%86164.1
    for fi in fx:
        r = exratio(fi)
        rx.append(r)
    for fi in fy:
        r = exratio(fi)
        ry.append(r)
    rx = flagr(rx)
    ry = flagr(ry)
    return { 'obs': obs, 'day': day, 'sid': sid, 'rx': rx, 'ry': ry }

def scatter_snr(day,sid,r,flag=False):
    for d in np.unique(day):
        ind = np.where(day==d)[0]
        if flag: plt.scatter(sid[ind],r[ind])
        else: plt.scatter(sid[ind],r.data[ind])
    plt.axvline(x=24214.35,color='black',linewidth=0.8,linestyle='-.')
    plt.axvline(x=26292.00,color='black',linewidth=0.8,linestyle='-.')
    plt.axvline(x=28392.95,color='black',linewidth=0.8,linestyle='-.')
    plt.axvline(x=30322.70,color='black',linewidth=0.8,linestyle='-.')
    plt.axvline(x=32141.35,color='black',linewidth=0.8,linestyle='-.')
    plt.axvline(x=33951.70,color='black',linewidth=0.8,linestyle='-.')
    plt.axvline(x=35810.80,color='black',linewidth=0.8,linestyle='-.')
    plt.grid(True,axis='y')
    plt.xlim(np.min(sid)-56,np.max(sid)+56)
    plt.xlabel('lst(s)')
    plt.ylabel('power/noise')
    plt.show()

