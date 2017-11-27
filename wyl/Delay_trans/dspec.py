import aipy as a, numpy as n
import capo as C
import numpy as np #lets start the transition away from n
def wedge_width_by_bl(aa, sdf, nchan, offset=0., horizon=1.):
    '''Generate dict of (upper,lower) bounds of delay-domain wedge for 
    baseline (i,j).  Bounds are given in bin number, derived from sdf (the
    delta freq between channels in GHz) and nchan.  offset (in ns) adds to
    the width of the wedge.  horizon is the fraction of the bl length to
    use, with 1 being the horizon and 0 being zenith.'''
    filters = {}
    for i in range(len(aa.ants)):
      for j in range(len(aa.ants)):
        if j < i: continue
        bl_len = aa.get_baseline(i,j)
        bl_len = horizon * n.sqrt(n.dot(bl_len,bl_len)) + offset
        filters[(i,j)] = wedge_width(bl_len, sdf, nchan)
    return filters

def wedge_width(bl_len, sdf, nchan):
    '''Return the (upper,lower) delay bins that geometrically
    correspond to the sky, for bl_len in ns, sdf in GHz, and the number
    of channels nchan.'''
    bin_dly = 1. / (sdf * nchan)
    uthresh, lthresh = bl_len/bin_dly + 1.5, -bl_len/bin_dly - 0.5
    uthresh, lthresh = int(n.ceil(uthresh)), int(n.floor(lthresh))
    return (uthresh,lthresh)
    
def delay_filter(data, wgts, uthresh, lthresh, tol=1e-4, window='none', maxiter=100):
    '''Apply a wideband delay filter to data.  Data are weighted again by wgts,
    windowed, Fourier transformed, and deconvolved allowing clean components
    between lthresh and uthresh.  The mdl, residual, and info are returned in
    frequency domain.'''
    window = a.dsp.gen_window(data.size, window=window)
    _d = n.fft.ifft(data * wgts * window)
    _w = n.fft.ifft(wgts * wgts * window)
    area = n.ones(_d.size, dtype=n.int); area[uthresh:lthresh] = 0
    _d_cl, info = a.deconv.clean(_d, _w, area=area, tol=tol, stop_if_div=False, maxiter=maxiter)
    d_mdl = n.fft.fft(_d_cl)# + info['res'])
    d_res = data - d_mdl * wgts
    return d_mdl, d_res, info

def phs2lstbin(aa, data, i, j, jd, lst_res=C.pspec.LST_RES):
    '''Phase data to the closest lst bin, which is the point that transits zenith
    at the sidereal time corresponding to the center of an lst bin.'''
    if aa.get_jultime() != jd: aa.set_jultime(jd)
    lst = aa.sidereal_time()
    ubin,vbin,lstbin = C.pspec.bin2uv(C.pspec.uv2bin(0,0,lst,lst_res=lst_res), lst_res=lst_res)
    zen = a.phs.RadioFixedBody(lstbin, aa.lat, epoch=aa.date)
    zen.compute(aa)
    data = aa.phs2src(data,zen,i,j)
    # XXX should we conjugate to u >= 0?
    # u = aa.get_uvw(i,j,src=zen)[0].flatten()[-1]
    # if u < 0: data = n.conj(data)
    return data
    

def apply_delay_filter(aa, data, wgts, i, j, wedges, phs2lst=False, jd=None, 
        skip_wgt=0.5, **kwargs):
    uthresh,lthresh = wedges[(i,j)]
    assert(data.ndim in [1,2])
    if data.ndim == 1: # apply only once and return
        if phs2lst: data = phs2lstbin(aa, data, i, j, jd, **kwargs)
        return delay_filter(data, wgts, uthresh, lthresh, **kwargs)
    else: # apply to each spectrum 
        dmdl, dres, info = n.zeros_like(data), n.zeros_like(data), {}
        for t in xrange(data.shape[0]): # loop through time dimension
            dt,wt = data[t], wgts[t]
            if n.average(wt) < skip_wgt: continue
            if phs2lst: dt = phs2lstbin(aa, dt, i, j, jd[t], **kwargs)
            dtm, dtr, i = delay_filter(dt, wt, uthresh, lthresh, **kwargs)
            dmdl[t] = dtm; dres[t] = dtr; info.update(i)
        return dmdl, dres, info
def delayfiltercov(C,horizon_bins=5,eig_cut_dnr=2):
    #delay filter a spectral covariance matrix
    #horizon_bins = distance delay=0 to be retained, ie the size of the wedge in bins
    # eig_cut_dnr = retain eigenvalues with a dynamic range of  median(dnr)*eig_cut_dnr 
    # where dnr is max(dspec eigenvector)/mean(abs(dpsec eigenvector outside horizon))    
    #
    # returns filtered_covariance,matching_projection matrix
    S,V = np.linalg.eig(C)
    dV = np.fft.ifft(V,axis=0)
    #calculate eigenvalue cut, selecting only those eigenvectors with strong delay spectrum signals
    dnr = np.max(np.abs(dV),axis=0)/np.mean(np.abs(dV)[horizon_bins:-horizon_bins,:],axis=0)
    median_dnr = np.median(dnr)
    eig_cut_dnr *= median_dnr
    S[dnr<eig_cut_dnr] = 0 #apply eigenvalue cut
    #mask outside wedge
    dV[horizon_bins:-horizon_bins,:] = 0 # mask out stuff outside the horizon
    V_filtered = np.fft.fft(dV,axis=0)
    #return filtered covariance and its matching projection matrix
    return np.einsum('ij,j,jk',V_filtered,S,V_filtered.T),np.einsum('ij,j,jk',V_filtered,S!=0,V_filtered.T)

