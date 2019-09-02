import numpy as np, sys, glob, copy, matplotlib
from scipy.io.idl import readsav
from plot_vis import *
import scipy.special as sp

colorsets = [plt.cm.Set1, plt.cm.Set2, plt.cm.Set3]

def erf(x): return sp.erf(x)

def erfinv(y): return sp.erfinv(y)

def getcolor(n):
    ii = n%29
    if ii<9: return colorsets[0](ii)
    elif ii<17: return colorsets[1](ii-9)
    else: return colorsets[2](ii-17)

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
        q[:,29] = False
        q[:,30] = False
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
    try:
       dp = diff['power_diff']*(h**3)
       dw = diff['weight_diff']
    except:
       dp = diff['power_3d']*(h**3)
       dw = diff['weights_3d']
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

def get_3d_ps(powerfn, bin_dict, slots):
    d = readsav(powerfn)
    h = d['hubble_param']
    slots_1d = slots/h
    try:
        p = np.copy(d['power_3d'])*(h**3)
        w = np.copy(d['weights_3d'])
    except:
        p = np.copy(d['power_diff'])*(h**3)
        w = np.copy(d['weight_diff'])
    kx = d['kx_mpc']/h
    ky = d['ky_mpc']/h
    kz = d['kz_mpc']/h
    p *= bin_dict['3d_mask']
    w *= bin_dict['3d_mask']
    p1d = np.zeros(slots_1d.size-1)
    w1d = np.zeros(slots_1d.size-1)
    nx = kx.size
    ny = ky.size
    nz = kz.size
    def fill(n):
        ix = n%nx
        iy = (n/nx)%ny
        iz = n/nx/ny
        k = np.sqrt(kx[ix]**2+ky[iy]**2+kz[iz]**2)
        ind = np.where(slots_1d<k)[0][-1]
        try:
            p1d[ind] += p[iz,iy,ix]*w[iz,iy,ix]
            w1d[ind] += w[iz,iy,ix]*(p[iz,iy,ix]!=0)
        except: pass
    map(fill, np.arange(nx*ny*nz))
#    for ix in range(kx.size):
#        for iy in range(ky.size):
#            for iz in range(kz.size):
#                k = np.sqrt(kx[ix]**2+ky[iy]**2+kz[iz]**2)
#                ind = np.where(slots_1d<k)[0][-1]
#                try:
#                    p1d[ind] += p[iz,iy,ix]*w[iz,iy,ix]
#                    w1d[ind] += w[iz,iy,ix]*(p[iz,iy,ix]!=0)
#                except: pass
    ind = np.where(w1d!=0)[0]
    p1d[ind] /= w1d[ind]
    kp = 0.5*(slots_1d[1:]+slots_1d[:-1]) 
    p1d *= 2
    p1d *= kp**3/2/np.pi**2
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

def power1dplot(fn, colornum=0, kmeasure = False, show = False, lgd=False, fiducial=True, label=''):
    lgd_hand = []
    fsp = fn.split('/')[-1].split('_') 
    tp = fsp[9]
    pol = fsp[10]
    d = readsav(fn)
    deor = readsav('/users/wl42/IDL/FHD/catalog_data/eor_power_1d.idlsave')
    h = d['hubble_param']
    p = d['power']*(h**3)
    #n = d['noise']*(h**3)
    s = (h**3)/np.sqrt(d['weights'])
    k_edges = d['k_edges']/h
    k = k_edges[1:]/2+k_edges[:-1]/2
    keor = deor['k_centers']/h
    peor = deor['power']*(h**3)
    p = p*k**3/2/np.pi**2
    #n = n*k**3/2/np.pi**2
    s = s*k**3/2/np.pi**2
    peor = peor*keor**3/2/np.pi**2
    if fiducial: plt.step(keor,peor,where='mid',label='fiducial theory',c='r')
    if kmeasure: plt.step(k,p,where='mid',label='measured power',linestyle=':',c=getcolor(colornum))
    plt.step(k,s,where='mid',label=label+'thermal noise',linestyle='--',c=getcolor(colornum))
    k_bin = xtostep(k)
    p_bin = ytostep(p)
    s_bin = ytostep(s)
    pkup = np.sqrt(2)*s_bin*erfinv(0.977-(1-0.977)*erf(p_bin/s_bin/np.sqrt(2))) + p_bin
    #ind0 = np.where(pkup<0)
    #pkup[ind0] = 2*s_bin[ind0]
    #pup = np.sqrt(2)*s_bin*scipy.special.erfinv(0.977-(1-0.977)*scipy.special.erf(p_bin/s_bin/np.sqrt(2)))
    plt.fill_between(k_bin,p_bin-2*s_bin,pkup,color='silver',alpha=0.8)
    plt.plot(k_bin,pkup,label=label+'2 sigma upper limit',c=getcolor(colornum))
    plt.ylabel('$\Delta^2$ ($\mathrm{mK^2}$)')
    plt.xlabel('$k$ ($h\mathrm{Mpc^{-1}}$)')
    plt.grid(True,axis='y')
    if lgd: plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left")
    plt.xlim(1e-1,np.max(k))
    plt.ylim(1e0,1e8)
    plt.xscale('log')
    plt.yscale('log')
    if show: plt.show()

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
    return k, p, pkup, s

def exratio(fi):
    #fi=obs+'_cubeXX__even_odd_joint_bh_res_'+pol+'_averemove_swbh_dencorr_no_horizon_wedge_kperplambda5-50_1dkpower.idlsave'
    d=readsav(fi)
    p=d['power']
    n=d['noise']
    k=d['k_edges']
    k=0.5*(k[1:]+k[:-1])
    h=d['hubble_param']
    k /= h
    sec=np.logical_and(k>0.2,k<0.35)
    return np.mean(p[sec])/np.mean(n[sec])*np.sqrt(np.sum(sec)-1)

def flagr(s):
    sc=np.ma.masked_array(s,np.zeros((len(s)),dtype=bool))
    for ii in range(sc.size):
        m=np.mean(sc)
        q=np.std(sc)
        ind=np.argmax(sc)
        if sc[ind]<m+3*q: break
        else: sc.mask[ind]=True
    print str(ii)+" observations flagged."
    return sc

lstfit = np.array([ 7.29211507e-05, -2.40663844e+00])

def metric_wrap(fpath='/users/wl42/data/wl42/FHD_out/fhd_int_PhaseII/ps/data/1d_binning/',filename=None):
    if filename is None:
        fx = glob.glob(fpath+'116*xx*no_horizon*1dkpower.idlsave')
        fy = glob.glob(fpath+'116*yy*no_horizon*1dkpower.idlsave')
    else:
        fx, fy = [], []
        for line in open(filename,'rb'):
            o = line.strip()
            fx.append(fpath+o+'_cubeXX__even_odd_joint_bh_res_xx_averemove_swbh_dencorr_no_horizon_wedge_kperplambda5-50_1dkpower.idlsave')
            fy.append(fpath+o+'_cubeXX__even_odd_joint_bh_res_yy_averemove_swbh_dencorr_no_horizon_wedge_kperplambda5-50_1dkpower.idlsave')
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
    z=lstfit
    return { 'obs': obs, 'day': day, 'sid': z[0]*sid+z[1], 'rx': rx, 'ry': ry }

def scatter_snr(dic,flag=False,mk='x',pointobslist=None,xex=None,yex=None,box=None,lw=0.8): 
    day = dic['day']
    sid = dic['sid']
    rx = dic['rx']
    ry = dic['ry']
    if xex:
        obsflag = []
        for line in open(xex,'rb'):
            obsflag.append(int(line.strip()))
        for o in obsflag:
            try:
                ind=np.where(dic['obs']==o)[0][0]
                rx[ind] = np.nan
            except: pass
    if yex:
        obsflag = []
        for line in open(yex,'rb'):
            obsflag.append(int(line.strip()))
        for o in obsflag:
            try:
                ind=np.where(dic['obs']==o)[0][0]
                ry[ind] = np.nan
            except: pass
    fig,axs=plt.subplots(2,1,sharex=True)
    fig.subplots_adjust(hspace=0.03)
    z=lstfit
    sid0 = [np.min(sid)-56*z[0], np.max(sid)+56*z[0]]
    lx = np.nanmean(rx)-3*np.nanstd(rx)
    ux = np.nanmean(rx)+3*np.nanstd(rx)
    ly = np.nanmean(ry)-3*np.nanstd(ry)
    uy = np.nanmean(ry)+3*np.nanstd(ry)
    if box:
        axs[0].plot(box,[lx,lx],linewidth=lw,c='k')
        axs[0].plot(box,[ux,ux],linewidth=lw,c='k')
        axs[0].plot([box[0],box[0]],[lx,ux],linewidth=lw,c='k')
        axs[0].plot([box[1],box[1]],[lx,ux],linewidth=lw,c='k')
        axs[1].plot(box,[ly,ly],linewidth=lw,c='k')
        axs[1].plot(box,[uy,uy],linewidth=lw,c='k')
        axs[1].plot([box[0],box[0]],[ly,uy],linewidth=lw,c='k')
        axs[1].plot([box[1],box[1]],[ly,uy],linewidth=lw,c='k')
    for d in np.unique(day):
        ind = np.where(day==d)[0]
        if flag: 
	    axs[0].scatter(sid[ind],rx[ind],marker=mk)
	    axs[1].scatter(sid[ind],ry[ind],marker=mk)
        else: 
	    axs[0].scatter(sid[ind],rx.data[ind],marker=mk)
	    axs[1].scatter(sid[ind],ry.data[ind],marker=mk)
    if pointobslist is None:
        axs[0].axvline(x=z[0]*24214.35+z[1],color='black',linewidth=0.8,linestyle='-.')
        axs[0].axvline(x=z[0]*26292.00+z[1],color='black',linewidth=0.8,linestyle='-.')
        axs[0].axvline(x=z[0]*28392.95+z[1],color='black',linewidth=0.8,linestyle='-.')
        axs[0].axvline(x=z[0]*30322.70+z[1],color='black',linewidth=0.8,linestyle='-.')
        axs[0].axvline(x=z[0]*32141.35+z[1],color='black',linewidth=0.8,linestyle='-.')
        axs[0].axvline(x=z[0]*33951.70+z[1],color='black',linewidth=0.8,linestyle='-.')
        axs[0].axvline(x=z[0]*35810.80+z[1],color='black',linewidth=0.8,linestyle='-.')
        axs[1].axvline(x=z[0]*24214.35+z[1],color='black',linewidth=0.8,linestyle='-.')
        axs[1].axvline(x=z[0]*26292.00+z[1],color='black',linewidth=0.8,linestyle='-.')
        axs[1].axvline(x=z[0]*28392.95+z[1],color='black',linewidth=0.8,linestyle='-.')
        axs[1].axvline(x=z[0]*30322.70+z[1],color='black',linewidth=0.8,linestyle='-.')
        axs[1].axvline(x=z[0]*32141.35+z[1],color='black',linewidth=0.8,linestyle='-.')
        axs[1].axvline(x=z[0]*33951.70+z[1],color='black',linewidth=0.8,linestyle='-.')
        axs[1].axvline(x=z[0]*35810.80+z[1],color='black',linewidth=0.8,linestyle='-.')
    else:
        #please make sure pointobslist is in lst order
        obslst = None
        for fi in pointobslist:
            l = []
            for line in open(fi, 'rb'):
		l.append(int(line.strip()))
            if obslst:
                s = 0.5*(obslst+l[0]%86164.1)
                axs[0].axvline(x=z[0]*s+z[1],color='black',linewidth=0.8,linestyle='-.')
                axs[1].axvline(x=z[0]*s+z[1],color='black',linewidth=0.8,linestyle='-.')
            obslst = l[-1]%86164.1
    #plt.grid(True,axis='y')
    axs[0].set_xlim(np.min(sid)-56*z[0],np.max(sid)+56*z[0])
    axs[1].set_xlim(np.min(sid)-56*z[0],np.max(sid)+56*z[0])
    axs[1].set_xlabel('lst(rad)')
    axs[0].set_ylabel('E-W')
    axs[1].set_ylabel('N-S')
    axs[0].set_yscale('symlog',linthreshy=1.)
    axs[1].set_yscale('symlog',linthreshy=1.)
    axs[0].set_ylim(-5,250)
    axs[1].set_ylim(-5,250)
    #plt.ylabel('power/noise')
    plt.show()

def plotdiffomnical(sz=(10,10),cl='C3'):
    fig,axs=plt.subplots(3,2,sharex=True,sharey=True,figsize=sz,gridspec_kw={'wspace': 0})
    deor = readsav('/users/wl42/IDL/FHD/catalog_data/eor_power_1d.idlsave')
    keor = deor['k_centers']
    peor = deor['power']
    z=['7.1','6.8','6.5']
    pol=['E-W','N-S']
    linpol=['xx','yy']
    chs=['0-95','48-143','96-191']
    fhd='/users/wl42/data/wl42/FHD_out/fhd_ps_sky/fhd_ps__sky_minus_red/data/1d_binning/'
    for ii in range(6):
        ix=ii/2
        iy=ii%2
        ch=chs[ix]
        _,_,pup,_=get_1d_limit('/users/wl42/data/wl42/FHD_out/fhd_ps_sky/ps/data/1d_binning/Combined_obs_cath0_cubeXX__even_odd_joint_fullimg_ch'+ch+'_res_xx_averemove_swbh_dencorr_no_120deg_wedge_cbw3_kperplambda12-50_kpar0.1-10_1dkpower.idlsave')
        ind = np.where(np.logical_or(np.isnan(pup),np.isinf(pup)))[0]
        fn = fhd+'Combined_obs_cath0_cubeXX__even_odd_joint_fullimg_ch'+ch+'_fullimg_ch'+ch+'_averemove_swbh_res_'+linpol[iy]+'_dencorr_no_horizon_wedge_1dkpower.idlsave'
        d = readsav(fn)
        h = d['hubble_param']
        k0 = keor/h
        p0 = peor*keor**3/2/np.pi**2 
        k,p,_,s=get_1d_limit(fn)
        p[ind] = np.nan
        p[np.where(p==0)] = np.nan
        axs[ix][iy].set_title('z='+z[ix]+' '+pol[iy])
        axs[ix][iy].set_xscale('log')
        axs[ix][iy].set_yscale('log')
        if iy==0: axs[ix][iy].set_ylabel('$\Delta^2$ ($\mathrm{mK^2}$)')
        if ix==2: axs[ix][iy].set_xlabel('$k$ ($h$ $\mathrm{Mpc^{-1}}$)')
        m1,m2,m3=None,None,None
        if ii==0: m1 = 'Power Difference'
        if ii==2: m2 = 'Fiducial EoR'
        if ii==4: m3 = '1-$\sigma$ Noise'
        axs[ix][iy].step(k,np.abs(p),where='mid',linestyle=':',c=cl)
        p[np.where(p<=0)] = np.nan
        axs[ix][iy].step(k,p,where='mid',c=cl,label=m1)    
        axs[ix][iy].step(k0,p0,where='mid',c='indigo',label=m2)
        axs[ix][iy].fill_between(xtostep(k),0,ytostep(s),color='silver',alpha=0.3,label=m3)
        axs[ix][iy].set_xlim(0.165,1.7)
        axs[ix][iy].set_ylim(10,1e6)
        axs[ix][iy].axvline(x=0.70,color='k',linewidth=0.8,linestyle='--')
        axs[ix][iy].axvline(x=1.05,color='k',linewidth=0.8,linestyle='-.')
        axs[ix][iy].grid(axis='y')
        if iy==0: axs[ix][iy].legend()
    plt.show()

def plotdiffbp(sz=(10,10)):
    fig,axs=plt.subplots(3,2,sharex=True,sharey=True,figsize=sz,gridspec_kw={'wspace': 0})
    z=['7.1','6.8','6.5']
    pol=['E-W','N-S']
    linpol=['xx','yy']
    chs=['0-95','48-143','96-191']
    fhd1='/users/wl42/data/wl42/FHD_out/fhd_ps_obsolete/fhd_ps__obsolete_minus_auto/data/1d_binning/'
    fhd2='/users/wl42/data/wl42/FHD_out/fhd_ps_obsolete/fhd_ps__obsolete_minus_wlcal/data/1d_binning/'
    for ii in range(6):
        ix=ii/2
        iy=ii%2
        ch=chs[ix]
        _,_,pup,_=get_1d_limit('/users/wl42/data/wl42/FHD_out/fhd_ps_sky/ps/data/1d_binning/Combined_obs_cath0_cubeXX__even_odd_joint_fullimg_ch'+ch+'_res_xx_averemove_swbh_dencorr_no_120deg_wedge_cbw3_kperplambda12-50_kpar0.1-10_1dkpower.idlsave')
        ind = np.where(np.logical_or(np.isnan(pup),np.isinf(pup)))[0]
        fn1 = fhd1+'Combined_obs_cath0_cubeXX__even_odd_joint_fullimg_ch'+ch+'_fullimg_ch'+ch+'_averemove_swbh_res_'+linpol[iy]+'_dencorr_no_horizon_wedge_1dkpower.idlsave'
        fn2 = fhd2+'Combined_obs_cath0_cubeXX__even_odd_joint_fullimg_ch'+ch+'_fullimg_ch'+ch+'_averemove_swbh_res_'+linpol[iy]+'_dencorr_no_horizon_wedge_1dkpower.idlsave'
        k1,p1,_,_=get_1d_limit(fn1)
        k2,p2,_,s=get_1d_limit(fn2)
        p1[ind] = np.nan
        p1[np.where(p1==0)] = np.nan
        p2[ind] = np.nan
        p2[np.where(p2==0)] = np.nan
        axs[ix][iy].set_title('z='+z[ix]+' '+pol[iy])
        axs[ix][iy].set_xscale('log')
        axs[ix][iy].set_yscale('log')
        if iy==0: axs[ix][iy].set_ylabel('$\Delta^2$ ($\mathrm{mK^2}$)')
        if ix==2: axs[ix][iy].set_xlabel('$k$ ($h$ $\mathrm{Mpc^{-1}}$)')
        m1,m2,m3=None,None,None
        if ii==0: m1 = 'Barry 2019'
        if ii==2: m2 = 'This Work'
        if ii==4: m3 = '1-$\sigma$ Noise'
        axs[ix][iy].step(k1,np.abs(p1),where='mid',linestyle=':',c='C9')
        p1[np.where(p1<=0)] = np.nan
        axs[ix][iy].step(k1,p1,where='mid',c='C9',label=m1)
        axs[ix][iy].step(k2,np.abs(p2),where='mid',linestyle=':',c='C3')
        p2[np.where(p2<=0)] = np.nan
        axs[ix][iy].step(k2,p2,where='mid',c='C3',label=m2)
        axs[ix][iy].fill_between(xtostep(k2),0,ytostep(s),color='silver',alpha=0.3,label=m3)
        axs[ix][iy].set_xlim(0.165,1.7)
        #axs[ix][iy].set_ylim(1,1e6)
        axs[ix][iy].axvline(x=0.70,color='k',linewidth=0.8,linestyle='--')
        axs[ix][iy].axvline(x=1.05,color='k',linewidth=0.8,linestyle='-.')
        axs[ix][iy].grid(axis='y')
        if iy==0: axs[ix][iy].legend()
    plt.show()
        
def plotlimit(f1,f2,sz=(10,10)):
    fig,axs=plt.subplots(3,2,sharex=True,sharey=True,figsize=sz,gridspec_kw={'wspace': 0})  
    deor = readsav('/users/wl42/IDL/FHD/catalog_data/eor_power_1d.idlsave')
    keor = deor['k_centers']
    peor = deor['power'] 
    z=['7.1','6.8','6.5'] 
    pol=['E-W','N-S']
    for ii in range(6): 
        ix=ii/2
        iy=ii%2
        fb=f1[ii]
        fl=f2[ii] 
        d=readsav(fl) 
        h=d['hubble_param']
        k0=keor/h
        p0=peor*keor**3/2/np.pi**2 
        kb,pb,pbup,sb=get_1d_limit(fb)
        kl,pl,plup,sl=get_1d_limit(fl)
        pl[np.where(pl==0)]=np.nan
        axs[ix][iy].set_xlim(0.1,2) 
        axs[ix][iy].set_ylim(10,1e7)
        axs[ix][iy].set_title('z='+z[ix]+' '+pol[iy])
        axs[ix][iy].set_xscale('log')
        axs[ix][iy].set_yscale('log') 
        m1,m2,m3=None,None,None 
        if iy==0: axs[ix][iy].set_ylabel('$\Delta^2$ ($\mathrm{mK^2}$)')
        if ix==2: axs[ix][iy].set_xlabel('$k$ ($h$ $\mathrm{Mpc^{-1}}$)')
        if ii==1: 
            m1='Fiducial Theory'
            m2='2016 2-$\sigma$ Upper Limit' 
            m3='2016 Noise Level' 
        axs[ix][iy].step(k0,p0,where='mid',c='r',label=m1) 
        axs[ix][iy].step(kb,pbup,where='mid',c='c',label=m2)
        axs[ix][iy].step(kb,sb,where='mid',c='c',linestyle='--',label=m3)
        l1,l2,l3=None,None,None
        if ii==0:
            l1='Measured Power' 
            l2='Noise Level'
            l3='2-$\sigma$ Upper Limit'
        axs[ix][iy].step(kl,pl,where='mid',c='k',label=l1) 
        axs[ix][iy].step(kl,sl,where='mid',c='k',linestyle='--',label=l2) 
        axs[ix][iy].step(kl,plup,where='mid',c='indigo',label=l3)
        axs[ix][iy].fill_between(xtostep(kl),ytostep(pl)-2*ytostep(sl),ytostep(plup),color='silver',alpha=0.8)
        axs[ix][iy].grid(axis='y') 
    plt.subplots_adjust(top=0.94,bottom=0.12,left=0.07,right=0.965,hspace=0.21,wspace=0.005)
    axs[0][0].legend(loc=2,ncol=3,fontsize='x-small')
    axs[0][1].legend(loc=2,ncol=3,fontsize='x-small')
    plt.show()

def diffcalup():
    crt = [1.0272699078325964,1.030787170292735]
    pref='/users/wl42/data/wl42/FHD_out/fhd_ps_'
    f1='/ps/data/1d_binning/Combined_obs_cath0_cubeXX__even_odd_joint_fullimg_ch0-191_res_'
    f2='_averemove_swbh_dencorr_no_120deg_wedge_cbw3_kperplambda12-50_kpar0.1-10_1dkpower.idlsave'
    def gfn(s,p): 
        k,p,pup,s=get_1d_limit(pref+s+f1+p+f2)
        return k,p,pup,s
    fl = ['red','sky','obsolete','auto']
    pl = [[],[]]
    polarizations = ['xx','yy']
    for f in fl:
         for i in range(2):
             k,_,p,s=gfn(f, polarizations[i])
             if f=='auto':
                 p *= crt[i]
                 s *= crt[i]
             pl[i].append(p)
    pols = {'xx': 'E-W', 'yy': 'N-S'}
    lgds = [['Hybrid Calibration',None],['Diff Beardsley 2016',None],['Diff Barry 2019',None],[None,'1-$\sigma$ Thermal Noise'],[None,'Fiducial EoR']]
    fig,axs=plt.subplots(1,2,sharex=True,sharey=True) 
    d0 = readsav('/users/wl42/IDL/FHD/catalog_data/eor_power_1d.idlsave')
    k0, f0 = d0['k_centers'], d0['power']
    f0 = f0*(k0**3) / (2*np.pi*np.pi)
    k0 = k0/0.71
    for ii in range(2):
        pol = polarizations[ii]
        ps1 = pl[ii][1] - pl[ii][0]
        ps2 = pl[ii][2] - pl[ii][1]
        ps3 = pl[ii][3] - pl[ii][1]
        axs[ii].step(k,np.abs(ps1),where='mid',linestyle=':',c='C2')
        ps1[np.where(ps1<=0)] = np.nan
        axs[ii].step(k,ps1,where='mid',c='C2',label=lgds[0][ii])
        axs[ii].step(k,np.abs(ps2),where='mid',linestyle=':',c='C3')
        ps2[np.where(ps2<=0)] = np.nan
        axs[ii].step(k,ps2,where='mid',c='C3',label=lgds[1][ii])
        axs[ii].step(k,np.abs(ps3),where='mid',linestyle=':',c='C9')
        ps3[np.where(ps3<=0)] = np.nan
        axs[ii].step(k,ps3,where='mid',c='C9',label=lgds[2][ii])
        axs[ii].fill_between(xtostep(k),0,ytostep(s),color='silver',alpha=0.5,label=lgds[3][ii])
        axs[ii].step(k0,f0,where='mid',c='indigo',label=lgds[4][ii])
        axs[ii].set_xlim(0.1,2)
        axs[ii].set_ylim(1,1e7)
        axs[ii].set_title(pols[pol])
        axs[ii].set_xlabel('$k$ ($h\mathrm{Mpc^{-1}}$)')
        axs[ii].set_xscale('log')
        axs[ii].set_yscale('log')
        axs[ii].legend()
        axs[ii].axvline(x=0.70,color='k',linewidth=0.8,linestyle='--')
        axs[ii].axvline(x=1.05,color='k',linewidth=0.8,linestyle='-.')
        axs[ii].grid(axis='y')
    axs[0].set_ylabel('$\Delta^2$ ($\mathrm{mK^2}$)')
    plt.tight_layout()
    plt.show()

def diffcalps(ch='0-191'):
    _,_,pup,_=get_1d_limit('/users/wl42/data/wl42/FHD_out/fhd_ps_obsolete/ps/data/1d_binning/Combined_obs_cath0_cubeXX__even_odd_joint_fullimg_ch'+ch+'_res_xx_averemove_swbh_dencorr_no_120deg_wedge_cbw3_kperplambda12-50_kpar0.1-10_1dkpower.idlsave')
    ind = np.where(np.logical_or(np.isnan(pup),np.isinf(pup)))[0]
    p1='/users/wl42/data/wl42/FHD_out/fhd_ps_sky/fhd_ps__sky_minus_red/data/1d_binning/'
    p2='/users/wl42/data/wl42/FHD_out/fhd_ps_obsolete/fhd_ps__obsolete_minus_wlcal/data/1d_binning/'
    p3='/users/wl42/data/wl42/FHD_out/fhd_ps_obsolete/fhd_ps__obsolete_minus_auto/data/1d_binning/'
    fig,axs=plt.subplots(1,2,sharex=True,sharey=True,gridspec_kw={'wspace': 0})
    pols = {'xx': 'E-W', 'yy': 'N-S'}
    lgds = [['Diff no Hybrid Cal',None],['This work',None],['Barry 2019',None],[None,'1-$\sigma$ Thermal Noise'],[None,'Fiducial EoR']]
    polarizations = ['xx', 'yy']
    d0 = readsav('/users/wl42/IDL/FHD/catalog_data/eor_power_1d.idlsave')
    k0, f0 = d0['k_centers'], d0['power']
    f0 = f0*(k0**3) / (2*np.pi*np.pi)
    k0 = k0/0.71
    for ii in range(2):
        pol = polarizations[ii]
        #k,ps1,_,s=get_1d_limit(p1+'Combined_obs_cath0_cubeXX__even_odd_joint_fullimg_ch'+ch+'_fullimg_ch'+ch+'_averemove_swbh_res_'+pol+'_dencorr_no_horizon_wedge_1dkpower.idlsave')
        k,ps2,_,s=get_1d_limit(p2+'Combined_obs_cath0_cubeXX__even_odd_joint_fullimg_ch'+ch+'_fullimg_ch'+ch+'_averemove_swbh_res_'+pol+'_dencorr_no_horizon_wedge_1dkpower.idlsave') 
        k,ps3,_,_=get_1d_limit(p3+'Combined_obs_cath0_cubeXX__even_odd_joint_fullimg_ch'+ch+'_fullimg_ch'+ch+'_averemove_swbh_res_'+pol+'_dencorr_no_horizon_wedge_1dkpower.idlsave')
        c1='C2'
        c2='C3'
        c3='C9'
        #print ps1.shape, ps2.shape,ps3.shape,pup.shape
        #ps1[ind] = np.nan
        ps2[ind] = np.nan
        ps3[ind] = np.nan
        #ps1[np.where(ps1==0)] = np.nan
        ps2[np.where(ps2==0)] = np.nan
        ps3[np.where(ps3==0)] = np.nan
        #ps2[np.where(np.isnan(ps1))] = np.nan
        #ps3[np.where(np.isnan(ps1))] = np.nan
        #axs[ii].step(k,np.abs(ps1),where='mid',linestyle=':',c=c1)
        #ps1[np.where(ps1<=0)]=np.nan
        #axs[ii].step(k,ps1,where='mid',c=c1,label=lgds[0][ii])
        axs[ii].step(k,np.abs(ps2),where='mid',linestyle=':',c=c2)
        ps2[np.where(ps2<=0)]=np.nan
        axs[ii].step(k,ps2,where='mid',c=c2,label=lgds[1][ii])
        axs[ii].step(k,np.abs(ps3),where='mid',linestyle=':',c=c3)
        ps3[np.where(ps3<=0)]=np.nan
        axs[ii].step(k,ps3,where='mid',c=c3,label=lgds[2][ii])
        axs[ii].fill_between(xtostep(k),0,ytostep(s),color='silver',alpha=0.3,label=lgds[3][ii])
        #axs[ii].step(k0,f0,where='mid',c='indigo',label=lgds[4][ii])
        kind=np.where(np.isnan(ps2)==False)[0]
        axs[ii].set_xlim(k[kind[0]],k[kind[-1]])
        ##axs[ii].set_xlim(0.1,2)
        ##axs[ii].set_xticks([0.2,0.4,0.8,1.0])
        axs[ii].get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        #axs[ii].set_xticklabels(['$\mathrm{2 \times 10^{-1})$','$\mathrm{4 \times 10^{-1})$','$\mathrm{8 \times 10^{-1})$','$\mathrm{10^0}$'])
        #axs[ii].set_ylim(1,1e7)
        axs[ii].set_title(pols[pol])
        axs[ii].set_xlabel('$k$ ($h\mathrm{Mpc^{-1}}$)')
        axs[ii].set_xscale('log')
        axs[ii].set_yscale('log')
        axs[ii].legend()
        axs[ii].axvline(x=0.70,color='k',linewidth=0.8,linestyle='--')
        axs[ii].axvline(x=1.05,color='k',linewidth=0.8,linestyle='-.')
        axs[ii].grid(axis='y')
    axs[0].set_ylabel('$\Delta^2$ ($\mathrm{mK^2}$)')
    for ax in axs: ax.label_outer()
    #plt.tight_layout()
    plt.show()

def readrts(fn):
    f=open(fn)
    k,px,py,sx,sy=[],[],[],[],[]
    def conv(st):
        if st[-1].isdigit():
            return float(st)
        else: return np.nan
    for line in f:
        ll=line.strip().split()
        k.append(conv(ll[0]))
        px.append(conv(ll[1]))
        py.append(conv(ll[2]))
        sx.append(conv(ll[3]))
        sy.append(conv(ll[4]))
    def uplim(p,s):
        return np.sqrt(2)*s*erfinv(0.977-(1-0.977)*erf(p/s/np.sqrt(2))) + p
    k=np.array(k)
    py=np.array(py)
    px=np.array(px)
    sx=np.array(sx)/2
    sy=np.array(sy)/2
    return k,uplim(px,sx),uplim(py,sy),sx,sy

def valpipe(lw=2,alpha=0.8):
    #kr,prx,pry,srx,sry=readrts('/users/wl42/data/wl42/OBS0/wenyang_updated.txt')
    kr,prx,pry,srx,sry=readrts('/users/wl42/data/wl42/OBS0/2016_phaseII_RTS_119obs.txt')
    prx[np.where(kr>1.7)]=np.nan
    pry[np.where(kr>1.7)]=np.nan
    srx[np.where(kr>1.7)]=np.nan
    sry[np.where(kr>1.7)]=np.nan
    pr=[prx,pry]
    sr=[srx,sry]
    lgds = [['FHD-eppsilon',None,'1-$\sigma$ Noise',None],[None, 'RTS-CHIPS', None, '1-$\sigma$ Noise']]
    fig,axs=plt.subplots(1,2,sharex=True,sharey=True,gridspec_kw={'wspace': 0})
    pols = {'xx': 'E-W', 'yy': 'N-S'}
    polarizations = ['xx', 'yy']
    c1='#1f77b4'
    c2='#ff7f0e'
    for ii in range(2):
        pol = polarizations[ii]
        kf,_,pf,sf=get_1d_limit('/users/wl42/data/wl42/FHD_out/fhd_ps_wlcal/ps/data/1d_binning/Combined_obs_rtschips2_cubeXX__even_odd_joint_fullimg_ch0-191_res_'+pol+'_averemove_swbh_dencorr_no_120deg_wedge_cbw3_kperplambda12-50_kpar0.14-10_1dkpower.idlsave')
        axs[ii].step(kf,pf,where='mid',label=lgds[ii][0],c=c1,linewidth=lw,alpha=alpha)
        axs[ii].step(kr,pr[ii],where='mid',label=lgds[ii][1],c=c2,linewidth=lw,alpha=alpha)
        axs[ii].step(kf,sf,where='mid',linestyle='--',label=lgds[ii][2],c=c1,linewidth=lw,alpha=alpha)
        axs[ii].step(kr,sr[ii],where='mid',linestyle='--',label=lgds[ii][3],c=c2,linewidth=lw,alpha=alpha)
        axs[ii].set_xlim(0.1,2)
        axs[ii].set_title(pols[pol])
        axs[ii].set_xlabel('$k$ ($h$ $\mathrm{Mpc^{-1}}$)')
        axs[ii].set_xscale('log')
        axs[ii].set_yscale('log')
        axs[ii].grid(axis='y')
        axs[ii].legend()
    axs[0].set_ylabel('$\Delta^2$ ($\mathrm{mK^2}$)')
    plt.tight_layout()
    plt.show()

def P1vsP2(alpha=0.5,c1='C4',c2='C5'):
    p2 = '/users/wl42/data/wl42/FHD_out/fhd_INTEGRATED/ps/data/1d_binning/Combined_obs_ALL_cubeXX__even_odd_joint_fullimg_ch9-126_res_'
    p1 = '/users/wl42/data/wl42/FHD_out/barry2019/ps/data/1d_binning/Combined_obs_btl_noalltv_noocc4_cubeXX__even_odd_joint_fullimg_ch9-126_res_'
    suf = '_averemove_swbh_dencorr_no_120deg_wedge_cbw3_kperplambda12-50_kpar0.15-10_1dkpower.idlsave'
    lgds = [['Barry 2019',None,'1-$\sigma$ Noise',None],[None, 'This Work', None, '1-$\sigma$ Noise']]
    fig, axs = plt.subplots(1,2,sharex=True,sharey=True,gridspec_kw={'wspace': 0})
    pols = {'xx': 'E-W', 'yy': 'N-S'}
    polarizations = ['xx', 'yy']
    for ii in range(2):
        pol = polarizations[ii]
        kb,mb,pb,sb=get_1d_limit(p1+pol+suf)
        kl,ml,pl,sl=get_1d_limit(p2+pol+suf)
        axs[ii].step(kb,pb,where='mid',label=lgds[ii][0],c=c1)
        axs[ii].step(kl,pl,where='mid',label=lgds[ii][1],c=c2)
        axs[ii].step(kb,sb,where='mid',label=lgds[ii][2],c=c1,linestyle='--')
        axs[ii].step(kl,sl,where='mid',label=lgds[ii][3],c=c2,linestyle='--')
        axs[ii].fill_between(xtostep(kb),ytostep(mb-2*sb),ytostep(pb),color=c1,alpha=alpha)
        axs[ii].fill_between(xtostep(kl),ytostep(ml-2*sl),ytostep(pl),color=c2,alpha=alpha)
        axs[ii].set_xlim(0.1,2)
        axs[ii].set_title(pols[pol])
        axs[ii].set_xlabel('$k$ ($h$ $\mathrm{Mpc^{-1}}$)')
        axs[ii].set_xscale('log')
        axs[ii].set_yscale('log')
        axs[ii].grid(axis='y')
        axs[ii].legend()
    axs[0].set_ylabel('$\Delta^2$ ($\mathrm{mK^2}$)')
    plt.tight_layout()
    plt.show()

def aelsim(T=13.515464586088115,df=8e4,fc=182.395e6): #T in mK. df: fq resolution. fc: center fq. Outputs the power in unit of mK^2Mpc^(-3)
    k = 4*np.pi/(12*256**2)
    dcv = 2.84947343e+10
    f0 = 1420405751.7667
    return k*df*dcv*(T**2)*f0/fc/fc

def plotssins(fn):
    d=np.load(fn)
    ins=mp2cal.qltm.INS(ins_arr=d)
    fq=np.linspace(167.035,197.755,768)
    fig=plt.figure(figsize=(10,6))
    p1=fig.add_subplot(3,1,1)
    i1=p1.imshow(ins.ins[:,0,:,0], aspect='auto', extent=(fq[0],fq[-1],len(d)-1,0), cmap='coolwarm')
    plt.colorbar(i1)
    frac_diff = ins.ins / np.ma.median(ins.ins, axis=0) - 1
    frac_diff /= np.std(frac_diff)
    p1.set_title('Sky-Subtracted Incoherent Noise Spectrum (SSINS)')
    p1.set_ylabel('Time Pairs')
    p2=fig.add_subplot(3,1,2)
    i2=p2.imshow(frac_diff[:,0,:,0], aspect='auto', extent=(fq[0],fq[-1],len(d)-1,0), cmap='coolwarm', clim=(-5,5))
    plt.colorbar(i2)
    p2.set_title('Fractional SSINS')
    p2.set_ylabel('Time Pairs')
    ins.outliers_flagging()
    ins.time_flagging()
    ins.coherence_flagging()
    ins.merge_flagging()
    frac_diff = ins.ins / np.ma.median(ins.ins, axis=0) - 1
    frac_diff /= np.std(frac_diff)
    p3=fig.add_subplot(3,1,3) 
    i3=p3.imshow(frac_diff[:,0,:,0], aspect='auto', extent=(fq[0],fq[-1],len(d)-1,0), cmap='coolwarm', clim=(-5,5))
    plt.colorbar(i3)
    p3.set_title('Fractional SSINS, Post Flagging')
    p3.set_xlabel('Frequency (MHz)')
    p3.set_ylabel('Time Pairs')
    plt.tight_layout()
    plt.show()

def plotchissins(fnchi,fnssins):
    dssins = np.load(fnssins)
    dchi=np.load(fnchi)
    fq=np.linspace(167.035,197.755,768)
    fig=plt.figure(figsize=(15,6))
    frac_diff = dssins / np.ma.median(dssins, axis=0) - 1
    s = np.std(frac_diff)
    p1=fig.add_subplot(3,1,1)
    i1=p1.imshow(frac_diff[:,0,:,0], aspect='auto', extent=(fq[0],fq[-1],len(dssins)-1,0), cmap='coolwarm',clim=(-5*s,5*s))
    plt.colorbar(i1)
    p1.set_title('Fractional SSINS')
    p1.set_ylabel('Time Steps')
    chi=mp2cal.qltm.Chisq('xx',dchi)
    chi0=np.ma.masked_array(dchi['chisq'],arr)
    m = np.ma.median(chi0)
    s = np.std(chi0)
    p2=fig.add_subplot(3,1,2)
    i2=p2.imshow(chi0, aspect='auto', extent=(fq[0],fq[-1],len(chi0)-1,0), cmap='coolwarm', clim=(m-5*s,m+5*s))
    plt.colorbar(i2)
    p2.set_title('$\chi^2$')
    p2.set_ylabel('Time Steps')
    chi.ch7det()
    chi1=np.ma.masked_array(dchi['chisq'],chi.chi.mask)
    m = np.ma.median(chi1)
    s = np.std(chi1)   
    p3=fig.add_subplot(3,1,3)
    i3=p3.imshow(chi1, aspect='auto', extent=(fq[0],fq[-1],len(chi1)-1,0), cmap='coolwarm', clim=(m-5*s,m+5*s))
    plt.colorbar(i3)
    p3.set_title('$\chi^2$, Post Flagging')
    p3.set_ylabel('Time Steps')
    p3.set_xlabel('Frequency (MHz)')
    plt.tight_layout()
    plt.show()

def pltsim(fn,T=0,show=True,lgd=None):
    d = readsav(fn)
    k = 0.5*(d['k_edges'][1:]+d['k_edges'][:-1])
    p = d['power']
    plt.plot(k,p,label=lgd)
    plt.xlabel('$k$ ($\mathrm{Mpc^{-1}}$)')
    plt.ylabel('$P_k$ ($mK^2Mpc^3$)')
    if T>0: 
        plt.plot(k, aelsim(T=T)*np.ones_like(k),label='Input EoR')
        print(np.mean(p[1:]/aelsim(T=T)))
    plt.legend()
    plt.xlim(k[5],k[-1])
    plt.yscale('log')
    if show: plt.show()

def plotbps(cable1,cable2,gbp1,gbp2):
    fq=np.linspace(167.055,197.735,768)
    def getb(arr):
        b=arr/pfb
        b[np.where(arr.mask)]=np.nan
        return b.data
    fig,axs=plt.subplots(4,1,sharex=True,sharey=True)
    for k in cable1.keys():
        axs[0].plot(fq,getb(np.mean(cable1[k],axis=0)),label=str(k)+' m')
        axs[0].legend(ncol=4)
    axs[0].set_title('Per cable type averged $|\hat{g}_i|$')
    axs[0].set_xticklabels([])
    lb = '$\hat{g}_i$'
    for ii in range(gbp1.shape[0]):
        axs[1].plot(fq,getb(gbp1[ii]),c='silver',label=lb)
        lb = None
    axs[1].plot(fq,getb(np.mean(gbp1,axis=0)),c='k',label='$<\hat{g}_i>_i$')
    axs[1].set_title('Globally averaged $|\hat{g}_i|$')
    axs[1].set_xticklabels([])
    axs[1].legend(ncol=2)
    for k in cable2.keys():
        axs[2].plot(fq,getb(np.mean(cable2[k],axis=0)),label=str(k)+' m')
        axs[2].legend(ncol=4)
    axs[2].set_title('Per cable type averaged $|\hat{g}_i|$/$A^o_i$')
    axs[2].set_xticklabels([])
    lb = '$\hat{g}_i$/$A^o_i$'
    for ii in range(gbp2.shape[0]):
        axs[3].plot(fq,getb(gbp2[ii]),c='silver',label=lb)
        lb = None
    axs[3].plot(fq,getb(np.mean(gbp2,axis=0)),c='k',label='$B_o$')
    axs[3].set_xlim(fq[0],fq[-1])
    axs[3].set_ylim(0.8,1.2)
    axs[3].set_title('Globally averaged $|\hat{g}_i|$/$A^o_i$')
    axs[3].set_xlabel('Frequency (MHz)')
    axs[3].legend(ncol=2)
    plt.show()
