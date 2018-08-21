import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

def downsample(masked_arr,dt=4,df=16):
    nt,nf,npol = masked_arr.shape
    assert nt%dt == 0 and nf%df==0, "invalid downsampling size"
    temp = masked_arr.reshape(nt/dt,dt,nf/df,df,npol)
    return np.mean(temp,axis=(1,3))

hori_filter = np.array([[ 1, 1, 1],
			[ 0, 0, 0],
			[-1,-1,-1]])
vert_filter = np.array([[ 1, 0,-1],
			[ 1, 0,-1],
			[ 1, 0,-1]])

def zero_pad(X,pad):
    X_pad = np.pad(X, ((pad,pad),(pad,pad),(0,0)),'constant', constant_values=((0,0),(0,0),(0,0)))
    return X_pad

def filt_arr(arr,filter_arr=hori_filter,stride=1,pad=0):
    arr_pad = zero_pad(arr,pad)
    nt,nf,n_pol = arr_pad.shape
    fh,fw = filter_arr.shape
    n_H = int((nt+2*pad-fh)/stride)+1
    n_W = int((nf+2*pad-fw)/stride)+1
    Z = np.zeros((n_H,n_W,n_pol))
    for h in range(n_H):
        for w in range(n_W):
            for p in range(n_pol):
                vert_start = h * stride
                vert_end = h * stride + fh
                horiz_start = w * stride
                horiz_end = w * stride + fw
                zz = arr[vert_start:vert_end,horiz_start:horiz_end,p]*filter_arr
                if np.sum(zz.mask) > 0: Z[h,w,p] = 0
                else: Z[h,w,p] = np.sum(zz)
    return Z

def stackomnichi(obs,path='/users/wl42/data/wl42/OBS0/sol_chi/'):
    dx = np.load(path+str(obs)+'.xx.omni.npz')
    dy = np.load(path+str(obs)+'.yy.omni.npz')
    cx = dx['chisq2']
    cy = dy['chisq2']
    mx = dx['flags']
    my = dy['flags']
    arr = np.zeros(cx.shape+(2,))
    flg = np.zeros(mx.shape+(2,),dtype=bool)
    arr[:,:,0] = cx
    arr[:,:,1] = cy
    flg[:,:,0] = mx
    flg[:,:,1] = my
    return np.ma.masked_array(arr,flg)

def getnoise(obs,path='/users/wl42/data/wl42/OBS1/mwilensky/Unflagged/all_spw/arrs/'):
    d = np.load(path+str(obs)+'_Unflagged_Amp_INS_frac_diff.npym')
    return d[:,0,:,:]

def flag_freqs(obs,path='/users/wl42/data/wl42/OBS1/mwilensky/All/all_spw/arrs/',nsig=6, f_thresh=0.15, t_thresh=0.25):
    d = np.load(path+str(obs)+'_All_Amp_INS_frac_diff_mask.npym')
    f = np.mean(d.mask,axis=(0,1,3))
    t = np.max(np.mean(d.mask,axis=(1,2)),axis=1)
    ind = np.where(f>f_thresh)[0]
    tnd = np.where(t>t_thresh)[0]
    for ii in ind: d.mask[:,:,ii,:] = True
    for tt in tnd: d.mask[tt,:,:,:] = True
    d2 = d[1:]*d[:-1]
    d2 = np.mean(d2,axis=(0,1,3))
    d2 -= np.mean(d2)
    d2 /= np.std(d2)
    fqind = np.where(d2 > nsig)
    return fqind, d2[fqind], tnd

def get_auto(uv):
    a1 = uv.ant_1_array[:uv.Nbls]
    a2 = uv.ant_2_array[:uv.Nbls]
    auto = {}
    ind = np.where(a1 == a2)[0]
    data = uv.data_array.reshape(uv.Ntimes,uv.Nbls,uv.Nspws,uv.Nfreqs,uv.Npols)
    scale = np.mean(data[:,ind,:,0:2].real)
    assert scale > 0, "Auto has to be positive"
    for ii in ind:
        a = str(a1[ii])
        auto[a] = np.zeros((uv.Ntimes,uv.Nspws,uv.Nfreqs,2))
        auto[a][:,0] = np.sqrt(data[:,ii,:,:,0].real/scale)
        auto[a][:,1] = np.sqrt(data[:,ii,:,:,1].real/scale)
        auto[a][np.where(auto[a]==0)] += 1e-10
    return auto

def flag_over_bls(uv):
    return np.sum(np.logical_not(uv.flag_array.reshape(uv.Ntimes,uv.Nbls,uv.Nspws,uv.Nfreqs,uv.Npols)),axis=1)==0

def smooth_over_t(INS,fc=1):
    frac_diff = INS / INS.mean(axis=0) - 1
    #time_ave = np.mean(frac_diff,axis=(1,2))
    #for i1 in range(INS.shape[1]):
    #    for i2 in range(INS.shape[2]):
    #        frac_diff[:,i1,i2,:] -= time_ave.data
    SH = INS.shape
    nc = SH[2]/24
    cf = int(fc*nc)
    m2 = np.zeros(SH)
    x = np.linspace(-cf,cf,2*cf+1)
    window = np.zeros((SH[0],SH[1],2*cf+1,SH[3]))
    for ii in range(2*cf+1): window[:,:,ii,:]=np.exp(-(x[ii]/float(cf))**2)
    for ff in range(SH[2]):
        min_ind = max(0,ff-cf)
        max_ind = min(SH[2],ff+cf+1)
        dm = np.sum(frac_diff[:,:,min_ind:max_ind,:]*window[:,:,cf-(ff-min_ind):cf+(max_ind-ff),:],axis=2)
        dn = np.sum(np.logical_not(frac_diff.mask[:,:,min_ind:max_ind,:])*window[:,:,cf-(ff-min_ind):cf+(max_ind-ff),:],axis=2)+1e-10
        m2[:,:,ff,:] = dm.data/dn
    return m2

def smooth_plotter(INS,fc=1):
    frac_diff = INS / INS.mean(axis=0) - 1
    fracsig = np.std(frac_diff)
    lims = (-5*fracsig, 5*fracsig)
    fig = plt.figure(figsize=(30,15))
    p1 = fig.add_subplot(4,2,1)
    i1 = p1.imshow(frac_diff[:,0,:,0],aspect='auto',cmap='coolwarm',clim=lims)
    plt.colorbar(i1)
    p1.set_title('XX masked')
    p2 = fig.add_subplot(4,2,2)
    i2 = p2.imshow(frac_diff[:,0,:,1],aspect='auto',cmap='coolwarm',clim=lims)
    plt.colorbar(i2)
    p2.set_title('YY masked')
    p3 = fig.add_subplot(4,2,3)
    i3 = p3.imshow(frac_diff[:,0,:,2],aspect='auto',cmap='coolwarm',clim=lims)
    plt.colorbar(i3)
    p3.set_title('XY masked')
    p4 = fig.add_subplot(4,2,4)
    i4 = p4.imshow(frac_diff[:,0,:,3],aspect='auto',cmap='coolwarm',clim=lims)
    plt.colorbar(i4)
    p4.set_title('YX masked')
    frac_diff = smooth_over_t(INS,fc=fc)
    fracsig = np.std(frac_diff)
    lims = (-5*fracsig, 5*fracsig)
    p5 = fig.add_subplot(4,2,5)
    i5 = p5.imshow(frac_diff[:,0,:,0],aspect='auto',cmap='coolwarm',clim=lims)
    plt.colorbar(i5)
    p5.set_title('XX smoothed')
    p6 = fig.add_subplot(4,2,6)
    i6 = p6.imshow(frac_diff[:,0,:,1],aspect='auto',cmap='coolwarm',clim=lims)
    plt.colorbar(i6)
    p6.set_title('YY smoothed')
    p7 = fig.add_subplot(4,2,7)
    i7 = p7.imshow(frac_diff[:,0,:,2],aspect='auto',cmap='coolwarm',clim=lims)
    plt.colorbar(i7)
    p7.set_title('XY smoothed')
    p8 = fig.add_subplot(4,2,8)
    i8 = p8.imshow(frac_diff[:,0,:,3],aspect='auto',cmap='coolwarm',clim=lims)
    plt.colorbar(i8)
    p8.set_title('YX post flagging')
    plt.tight_layout()

def time_flagging(INS,nthresh=1e-5):
    maxiter = INS.shape[0]
    for niter in range(maxiter):
        frac_smt = smooth_over_t(INS,fc=1)
        frac_co = frac_smt[:,:,1:,:]*frac_smt[:,:,:-1,:]
        dt_slice = np.max(np.mean(frac_co,axis=(1,2)),axis=1)
        dt_ind = np.argmax(np.abs(dt_slice))
        if np.abs(dt_slice[dt_ind]) > nthresh: INS.mask[dt_ind] = True
        else: break

def freq_flagging(INS,nsig=6):
    frac_diff = INS / INS.mean(axis=0) - 1
    df_slice = frac_diff[1:]*frac_diff[:-1]
    df_slice = np.mean(df_slice,axis=(0,1,3))
    df_slice -= np.mean(df_slice)
    df_ind = np.where(df_slice>nsig*np.std(df_slice))[0]
    INS.mask[:,:,df_ind,:] = True

def coherence_flagging(INS,fc=0.25,nsig=5,maxiter=10):
    SH = INS.shape
    nc = SH[2]/24
    cf = int(fc*nc)
    m2 = np.zeros(SH)
    x = np.linspace(-cf,cf,2*cf+1)
    window = np.zeros((SH[0],SH[1],2*cf+1,SH[3]))
    for ii in range(2*cf+1): window[:,:,ii,:]=np.exp(-(x[ii]/float(cf))**2)
    for niter in range(maxiter):
        frac_diff = INS / INS.mean(axis=0) - 1
        #time_ave = np.mean(frac_diff,axis=(1,2))
        #for i1 in range(INS.shape[1]):
        #    for i2 in range(INS.shape[2]):
        #        frac_diff[:,i1,i2,:] -= time_ave.data
        for ff in range(SH[2]): 
            min_ind = max(0,ff-cf)
            max_ind = min(SH[2],ff+cf+1)
            dm = np.sum(frac_diff[:,:,min_ind:max_ind,:]*window[:,:,cf-(ff-min_ind):cf+(max_ind-ff),:],axis=2)
            dn = np.sum(np.logical_not(frac_diff.mask[:,:,min_ind:max_ind,:])*window[:,:,cf-(ff-min_ind):cf+(max_ind-ff),:],axis=2)+1e-10
            m2[:,:,ff,:] = dm.data/dn
        sigma = np.std(m2)
        indf = np.where(np.abs(m2)>nsig*sigma)
        if indf[0].size == 0: break
        INS.mask[indf] = True

def extend_flagging(INS, f_thresh=0.5, t_thresh=0.5):
    mt = np.max(np.mean(INS.mask, axis=(1,2)),axis=1)
    indt = np.where(mt>t_thresh)[0]
    mf = np.max(np.mean(INS.mask, axis=(0,1)),axis=1)
    indf = np.where(mf>f_thresh)[0]
    INS.mask[indt,:,:,:] = True
    INS.mask[:,:,indf,:] = True

def frac_diff_plotter(INS,tf=True,ff=True,cf=True,ef=True):
    frac_diff = INS / INS.mean(axis=0) - 1
    fig = plt.figure(figsize=(30,15))
    p1 = fig.add_subplot(4,2,1)
    i1 = p1.imshow(frac_diff[:,0,:,0],aspect='auto',cmap='coolwarm')
    plt.colorbar(i1)
    p1.set_title('XX masked')
    p2 = fig.add_subplot(4,2,2)
    i2 = p2.imshow(frac_diff[:,0,:,1],aspect='auto',cmap='coolwarm')
    plt.colorbar(i2)
    p2.set_title('YY masked')
    p3 = fig.add_subplot(4,2,3)
    i3 = p3.imshow(frac_diff[:,0,:,2],aspect='auto',cmap='coolwarm')
    plt.colorbar(i3)
    p3.set_title('XY masked')
    p4 = fig.add_subplot(4,2,4)
    i4 = p4.imshow(frac_diff[:,0,:,3],aspect='auto',cmap='coolwarm')
    plt.colorbar(i4)
    p4.set_title('YX masked')
    if tf: time_flagging(INS)
    if ff: freq_flagging(INS)
    if cf: coherence_flagging(INS)
    if ef: extend_flagging(INS)
    frac_diff = INS / INS.mean(axis=0) - 1
    p5 = fig.add_subplot(4,2,5)
    i5 = p5.imshow(frac_diff[:,0,:,0],aspect='auto',cmap='coolwarm')
    plt.colorbar(i5)
    p5.set_title('XX post flagging')
    p6 = fig.add_subplot(4,2,6)
    i6 = p6.imshow(frac_diff[:,0,:,1],aspect='auto',cmap='coolwarm')
    plt.colorbar(i6)
    p6.set_title('YY post flagging')
    p7 = fig.add_subplot(4,2,7)
    i7 = p7.imshow(frac_diff[:,0,:,2],aspect='auto',cmap='coolwarm')
    plt.colorbar(i7)
    p7.set_title('XY post flagging')
    p8 = fig.add_subplot(4,2,8)
    i8 = p8.imshow(frac_diff[:,0,:,3],aspect='auto',cmap='coolwarm')
    plt.colorbar(i8)
    p8.set_title('YX post flagging')
    plt.tight_layout()
    #plt.show()

def further_flagging(UV, INS):
    mask1 = np.zeros((UV.Ntimes,UV.Nspws,UV.Nfreqs,UV.Npols),dtype=bool)
    mask2 = np.zeros((UV.Ntimes,UV.Nspws,UV.Nfreqs,UV.Npols),dtype=bool)
    mask1[0] = True
    mask1[1:] = INS.mask
    mask2[-1] = True
    mask2[:-1] = INS.mask
    mask_all = np.logical_and(mask1,mask2)
    for ii in range(UV.Nbls): UV.flag_array[ii::UV.Nbls] = np.logical_or(UV.flag_array[ii::UV.Nbls], mask_all)

