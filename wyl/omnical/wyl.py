import numpy as np, omnical, aipy
import subprocess, datetime, os
from astropy.io import fits
import copy
import heracal
from scipy.io.idl import readsav


def unwrap(arr):
    brr = np.unwrap(arr)
    crr = []
    for ii in range(1,brr.size): crr.append(brr[ii]-brr[ii-1])
    crr = np.unwrap(crr)
    nn = np.round(crr[0]/(2*np.pi))
    crr -= (nn*2.*np.pi)
    drr = np.zeros(brr.shape)+brr[0]
    for ii in range(crr.size): drr[ii+1] += np.sum(crr[:ii+1])
    return drr


def output_mask_array(flag_array):
    invf = 1 - flag_array
    sf = np.sum((np.sum(invf,axis=0)),axis=0).astype(bool)
    st = np.sum((np.sum(invf,axis=1)),axis=1).astype(bool)
    mask_array = 1 - np.outer(st,sf)
    mask_array = mask_array.astype(bool)
    return mask_array


def find_ex_ant(uvdata):
    ex_ant = []
    for ii in uvdata.antenna_numbers:
        if not ii in uvdata.ant_1_array and not ii in uvdata.ant_2_array:
            ex_ant.append(ii)
    return ex_ant


def scale_gains(g0, amp_ave=1.):
    g = copy.deepcopy(g0)
    for p in g.keys():
        amp = 0
        n = 0
        for a in g[p].keys():
            amp += np.abs(g[p][a])
            n += 1
        amp /= n
        q = amp/amp_ave
        inds = np.where(amp!=0)
        for a in g[p].keys(): g[p][a][inds] /= q[inds]
    return g


def uv_wrap_fc(uv,redbls,pols=['xx','yy']):
    wrap_list = []
    a1 = uv.ant_1_array[:uv.Nbls]
    a2 = uv.ant_2_array[:uv.Nbls]
    data = uv.data_array
    flag = uv.flag_array
    for jj in range(uv.Npols):
        pp = aipy.miriad.pol2str[uv.polarization_array[jj]]
        if not pp in pols: continue
        wrap = {}
        wrap['pol'] = pp
        wrap['data'] = {}
        wrap['flag'] = {}
        for ii in range(uv.Nbls):
            if (a1[ii],a2[ii]) in redbls: bl = (a1[ii],a2[ii])
            elif (a2[ii],a1[ii]) in redbls: bl = (a2[ii],a1[ii])
            else: continue
            if not wrap['data'].has_key(bl):
                if bl == (a1[ii],a2[ii]): dat_temp = data[:,0][:,:,jj].reshape(uv.Ntimes,uv.Nbls,uv.Nfreqs)[:,ii]
                else: dat_temp = data[:,0][:,:,jj].reshape(uv.Ntimes,uv.Nbls,uv.Nfreqs)[:,ii].conj()
                flg_temp = flag[:,0][:,:,jj].reshape(uv.Ntimes,uv.Nbls,uv.Nfreqs)[:,ii]
                dat_ma = np.ma.masked_array(dat_temp, mask=flg_temp)
                dat_ma = np.mean(dat_ma,axis=0)
                wrap['data'][bl] = {pp: np.complex64([dat_ma.data])}
                wrap['flag'][bl] = {pp: np.array([dat_ma.mask])}
        wrap_list.append(wrap)
    return wrap_list


def uv_wrap_omni(uv,pols=['xx','yy']):
    data_wrap = {}
    a1 = uv.ant_1_array[:uv.Nbls]
    a2 = uv.ant_2_array[:uv.Nbls]
    data = uv.data_array
    flag = uv.flag_array
    for jj in range(uv.Npols):
        pp = aipy.miriad.pol2str[uv.polarization_array[jj]]
        if not pp in pols: continue
        wrap = {}
        wrap['pol'] = pp
        wrap['data'] = {}
        wrap['flag'] = {}
        wrap['auto'] = {}
        wrap['mask'] = output_mask_array(flag[:,0][:,:,jj].reshape(uv.Ntimes,uv.Nbls,uv.Nfreqs))
        auto_scale = 0
        for ii in range(uv.Nbls):
            if a1[ii] == a2[ii]:
                auto_m = np.ma.masked_array(data[:,0][:,:,jj].reshape(uv.Ntimes,uv.Nbls,uv.Nfreqs)[:,ii].real,mask=wrap['mask'])
                wrap['auto'][a1[ii]] = np.sqrt(np.mean(auto_m,axis=0).data) + 1e-10
                auto_scale += np.nanmean(wrap['auto'][a1[ii]])
            else:
                bl = (a1[ii],a2[ii])
                wrap['data'][bl] = {pp: np.complex64(data[:,0][:,:,jj].reshape(uv.Ntimes,uv.Nbls,uv.Nfreqs)[:,ii])}
                wrap['flag'][bl] = {pp: np.array(flag[:,0][:,:,jj].reshape(uv.Ntimes,uv.Nbls,uv.Nfreqs)[:,ii])}
        auto_scale /= len(wrap['auto'].keys())
        for a in wrap['auto'].keys(): wrap['auto'][a] /= auto_scale
        data_wrap[pp] = wrap
    return data_wrap


def polyfunc(x,z):
    sum = np.zeros((x.size))
    for ii in range(z.size):
        sum *= x
        sum += z[ii]
    return sum


def mwa_bandpass_fit(gains0, auto, tile_info, amp_order=2, phs_order=1, fit_reflection=True):
    gains = copy.deepcopy(gains0)
    fqs = np.linspace(167.075,197.715,384)
    freq = np.arange(384)
    for p in gains.keys():
        for ant in gains[p].keys():
            x = np.where(gains[p][ant]!=0)[0]
            if x.size == 0: continue
            A = np.zeros((384),dtype=np.float)
            for n in range(0,24):
                chunk = np.arange(16*n+1,16*n+15)
                induse = np.where(gains[p][ant][chunk]!=0)
                z1 = np.polyfit(freq[chunk[induse]],np.abs(gains[p][ant][chunk[induse]])/auto[ant][chunk[induse]],amp_order)
                A[chunk[induse]] = auto[ant][chunk[induse]]*polyfunc(freq[chunk[induse]],z1)
            y2 = np.angle(gains[p][ant][x])
            y2 = np.unwrap(y2)
            z2 = np.polyfit(x,y2,phs_order)
            rp = np.zeros((384))
            cable = tile_info[ant]['cable']
            if fit_reflection and cable==150:
                vf = tile_info[ant]['vf']
                t0 = 2*cable/299792458.0/vf*1e6
                rp[x] = y2 - polyfunc(x,z2)
                tau = np.fft.fftfreq(384,(fqs[-1]-fqs[0])/383)
                fftrp = np.fft.fft(rp,n=384)
                inds = np.where(abs(np.abs(tau)-t0)<0.05)
                imax = np.argmax(np.abs(fftrp[inds]))
                ind = np.where(np.abs(tau)==np.abs(tau[inds][imax]))
                mask =np.zeros((384))
                mask[ind] = 1.
                fftrp *= mask
                rp = np.fft.ifft(fftrp)
            gains[p][ant][x] = A[x]*np.exp(1j*polyfunc(x,z2))
            gains[p][ant][x] *= np.exp(1j*rp[x])
    return gains


def poly_bandpass_fit(gains0,fit_order=4):
    gains = copy.deepcopy(gains0)
    for p in gains.keys():
        for a in gains[p].keys():
            g = np.copy(gains[p][a])
            for ff in range(24):
                chunk = np.arange(16*ff+1,16*ff+15)
                z1 = np.polyfit(chunk,g.real[chunk],fit_order)
                z2 = np.polyfit(chunk,g.imag[chunk],fit_order)
                gains[p][a][chunk] = polyfunc(chunk,z1) + 1j*polyfunc(chunk,z2)
    return gains


def amp_bandpass_fit(gains0,fit_order=4):
    gains = copy.deepcopy(gains0)
    for p in gains.keys():
        for a in gains[p].keys():
            g = np.abs(gains[p][a])
            for ff in range(24):
                chunk = np.arange(16*ff+1,16*ff+15)
                z = np.polyfit(chunk,g[chunk],fit_order)
                gains[p][a][chunk] = polyfunc(chunk,z)
    return gains


def ampproj(g_input,g_target):
    amppar = {}
    for p in g_input.keys():
        SH = g_input[p][g_input[p].keys()[0]].shape
        s = np.zeros(SH)
        n = np.zeros(SH)
        for a in g_input[p].keys():
            if not a in g_target[p].keys(): continue
            if np.isnan(np.mean(g_target[p][a])): continue
            if np.isnan(np.mean(g_input[p][a])): continue
            num = np.ones(SH)
            amp_in = np.abs(g_input[p][a])
            amp_ta = np.resize(np.abs(g_target[p][a]),SH)
            ind = np.where(amp_in==0)
            amp_in[ind] = 1.
            amp_ta[ind] = 0.
            num[ind] = 0
            s += amp_ta/amp_in
            n += num
        ind = np.where(n==0)
        n[ind] = 1.
        s[ind] = 0.
        amppar[p] = s/n
    return amppar


def phsproj(g_input,g_target,antpos,EastHex,SouthHex): #only returns slopes
    phspar = {}
    ax1,ax2 = [],[]
    for ii in range(EastHex.shape[0]):
        if ii == 3: continue
        ind_east = EastHex[ii]
        ind_south = SouthHex[ii]
        ax1.append(ind_east)
        ax1.append(ind_south)
    for jj in range(EastHex.shape[1]):
        if jj == 3: continue
        ind_east = EastHex[:,jj]
        ind_south = SouthHex[:,jj]
        ax2.append(ind_east)
        ax2.append(ind_south)
    for p in g_input.keys():
        phspar[p] = {}
        a0 = g_input[p].keys()[0]
        SH = g_input[p][a0].shape
        if len(SH) == 2:
            for a in g_input[p].keys(): g_input[p][a] = np.mean(g_input[p][a],axis=0)
        slp1 = []
        slp2 = []
        for ff in range(0,384):
            if ff%16 in [0,15]:
                slp1.append(0)
                slp2.append(0)
                continue
            #***** East-West direction fit *****#
            slope = []
            for inds in ax1:
                x,tau = [],[]
                for ii in inds:
                    if not ii in g_input[p].keys(): continue
                    if not ii in g_target[p].keys(): continue
                    if np.isnan(g_input[p][ii][ff]): continue
                    if np.isnan(g_target[p][ii][ff]): continue
                    x.append(float(np.argwhere(inds==ii)))
                    tau.append(np.angle(g_target[p][ii][ff]*g_input[p][ii][ff].conj()))
                if len(tau) < 3: continue
                if np.round(x[-1])-np.round(x[0])+1 != len(x): continue
                tau = unwrap(tau)
                z = np.polyfit(x,tau,1)
                slope.append(z[0])
            slope = np.unwrap(slope)
            slp1.append(np.median(slope))
            #***** 60 deg East-South direction fit *****#
            slope = []
            for inds in ax2:
                x,tau = [],[]
                for ii in inds:
                    if not ii in g_input[p].keys(): continue
                    if not ii in g_target[p].keys(): continue
                    if np.isnan(g_input[p][ii][ff]): continue
                    if np.isnan(g_target[p][ii][ff]): continue
                    x.append(float(np.argwhere(inds==ii)))
                    tau.append(np.angle(g_target[p][ii][ff]*g_input[p][ii][ff].conj()))
                if len(tau) < 3: continue
                if np.round(x[-1])-np.round(x[0])+1 != len(x): continue
                tau = unwrap(tau)
                z = np.polyfit(x,tau,1)
                slope.append(z[0])
            slope = np.unwrap(slope)
            slp2.append(np.median(slope))
        phspar[p]['phi1'] = np.array(slp1)
        phspar[p]['phi2'] = np.array(slp2)
    return phspar


def plane_fitting(gains,antpos):
    phspar = {}
    for p in gains.keys():
        phspar[p] = {}
        phix,phiy,offset_east,offset_south = [],[],[],[]
        for f in range(384):
            if f%16 in [0,15]:
                phix.append(0)
                phiy.append(0)
                offset_east.append(0)
                offset_south.append(0)
                continue
            M0 = np.zeros((4,4))
            p0 = np.zeros((4,1))
            for a in gains[p].keys():
                x = antpos[a]['top_x']
                y = antpos[a]['top_y']
                if gains[p][a].ndim == 2:
                    z = np.angle(np.mean(gains[p][a],axis=0)[f])
                else:
                    z = np.angle(gains[p][a][f])
                if 56 < a < 93:
                    M0 += np.array([[x*x, x*y, x , 0 ],
                                    [x*y, y*y, y , 0 ],
                                    [ x ,  y , 1 , 0 ],
                                    [ 0 ,  0 , 0 , 0 ]])
                    p0 += np.array([[z*x],
                                    [z*y],
                                    [ z ],
                                    [ 0 ]])
                if 92 < a < 128:
                    M0 += np.array([[x*x, x*y, 0 , x ],
                                    [x*y, y*y, 0 , y ],
                                    [ 0 ,  0 , 0 , 0 ],
                                    [ x ,  y , 0 , 1 ]])
                    p0 += np.array([[z*x],
                                    [z*y],
                                    [ 0 ],
                                    [ z ]])
            C = np.linalg.inv(M0).dot(p0)
            #Attention: append negative results here
            phix.append(-C[0][0])
            phiy.append(-C[1][0])
            offset_east.append(-C[2][0])
            offset_south.append(-C[3][0])
        phspar[p]['phix'] = np.array(phix)
        phspar[p]['phiy'] = np.array(phiy)
        phspar[p]['offset_east'] = np.array(offset_east)
        phspar[p]['offset_south'] = np.array(offset_south)
    return phspar


def degen_project_OF(gomni,gfhd,antpos,EastHex,SouthHex,v2={}):
    gains = copy.deepcopy(gomni)
    for p in gains.keys():
        ref1 = min(gains[p].keys())
        ref2 = max(gains[p].keys())
        ref_exp1 = np.exp(1j*np.angle(gains[p][ref1]*gfhd[p][ref1].conj()))
        ref_exp2 = np.exp(1j*np.angle(gains[p][ref2]*gfhd[p][ref2].conj()))
        for a in gains[p].keys():
            if a < 93: gains[p][a] /= ref_exp1
            else: gains[p][a] /= ref_exp2
        amppar = ampproj(gains,gfhd)
        phspar = phsproj(gains,gfhd,antpos,EastHex,SouthHex)
        for a in gains[p].keys():
            if a < 93:
                dx = antpos[a]['top_x']-antpos[ref1]['top_x']
                dy = antpos[a]['top_y']-antpos[ref1]['top_y']
            else:
                dx = antpos[a]['top_x']-antpos[ref2]['top_x']
                dy = antpos[a]['top_y']-antpos[ref2]['top_y']
            nx = dx/14.-dy/np.sqrt(3)/14.
            ny = -2*dy/np.sqrt(3)/14.
            proj = amppar[p]*np.exp(1j*(nx*phspar[p]['phi1']+ny*phspar[p]['phi2']))
            gains[p][a] *= proj
        ratio = {p:{}}
        for a in gains[p].keys():
            r = gains[p][a]*gfhd[p][a].conj()
            if np.isnan(np.mean(r)): continue
            ratio[p][a] = r
        phspar2 = plane_fitting(ratio,antpos)
        for a in gains[p].keys():
            dx = antpos[a]['top_x']
            dy = antpos[a]['top_y']
            proj = np.exp(1j*(dx*phspar2[p]['phix']+dy*phspar2[p]['phiy']))
            if a > 92: proj *= np.exp(1j*phspar2[p]['offset_south'])
            else: proj *= np.exp(1j*phspar2[p]['offset_east'])
            gains[p][a] *= proj
        if not v2 == {}:
            pp = p+p
            for bl in v2[pp].keys():
                i,j = bl
                if i < 93: v2[pp][bl] *= (ref_exp1*np.exp(-1j*phspar2[p]['offset_east']))
                else: v2[pp][bl] *= (ref_exp2*np.exp(-1j*phspar2[p]['offset_south']))
                if j < 93: v2[pp][bl] *= (ref_exp1.conj()*np.exp(1j*phspar2[p]['offset_east']))
                else: v2[pp][bl] *= (ref_exp2.conj()*np.exp(1j*phspar2[p]['offset_south']))
                dx = antpos[i]['top_x']-antpos[j]['top_x']
                dy = antpos[i]['top_y']-antpos[j]['top_y']
                nx = dx/14.-dy/np.sqrt(3)/14.
                ny = -2*dy/np.sqrt(3)/14.
                proj = amppar[p]*amppar[p]*np.exp(1j*(nx*phspar[p]['phi1']+ny*phspar[p]['phi2']))*np.exp(1j*(dx*phspar2[p]['phix']+dy*phspar2[p]['phiy']))
                proj = np.resize(proj,v2[pp][bl].shape)
                ind = np.where(proj!=0)
                v2[pp][bl][ind] /= proj[ind]
    return gains


def degen_project_FO(gomni,antpos,v2={}):
    gains = scale_gains(gomni)
    phspar = plane_fitting(gains,antpos)
    for p in gains.keys():
        for a in gains[p].keys():
            dx = antpos[a]['top_x']
            dy = antpos[a]['top_y']
            proj = np.exp(1j*(dx*phspar[p]['phix']+dy*phspar[p]['phiy']))
            if a > 92: proj *= np.exp(1j*phspar[p]['offset_south'])
            else: proj *= np.exp(1j*phspar[p]['offset_east'])
            gains[p][a] *= proj
        if not v2 == {}:
            pp = p+p
            for bl in v2[pp].keys():
                i,j = bl
                if i < 93: v2[pp][bl] *= np.exp(-1j*phspar[p]['offset_east'])
                else: v2[pp][bl] *= np.exp(-1j*phspar[p]['offset_south'])
                if j < 93: v2[pp][bl] *= np.exp(1j*phspar[p]['offset_east'])
                else: v2[pp][bl] *= np.exp(1j*phspar[p]['offset_south'])
                dx = antpos[i]['top_x']-antpos[j]['top_x']
                dy = antpos[i]['top_y']-antpos[j]['top_y']
                proj = np.exp(-1j*(dx*phspar[p]['phix']+dy*phspar[p]['phiy']))
                v2[pp][bl][ind] *= proj
    return gains


def degen_project_simple(g_input,g_target,antpos):
    g_output = copy.deepcopy(g_input)
    amppar = ampproj(g_input,g_target)
    for p in g_output.keys():
        ratio = {p:{}}
        for a in g_output[p].keys():
            r = g_input[p][a]*g_target[p][a].conj()
            if np.isnan(np.mean(r)): continue
            ratio[p][a] = r
        phspar = plane_fitting(ratio,antpos)
        for a in g_input[p].keys():
            dx = antpos[a]['top_x']
            dy = antpos[a]['top_y']
            proj = amppar[p]*np.exp(1j*(dx*phspar[p]['phix']+dy*phspar[p]['phiy']))
            if a > 92: proj *= np.exp(1j*phspar[p]['offset_south'])
            else: proj *= np.exp(1j*phspar[p]['offset_east'])
            g_output[p][a] *= proj
    return g_output


def cal_var_wgt(v,m,w):
    n = np.ma.masked_array(v-m,mask=w,fill_value=0.+0.j)
    var = np.var(n,axis=0).data
    zeros = np.where(var==0)
    var[zeros] = 1.
    inv = 1./var
    inv[zeros] = 0.
    return inv


def pos_to_info(position, pols=['x'], fcal=False, **kwargs):
    nant = position['nant']
    antpos = -np.ones((nant*len(pols),3))
    xmin,ymin = 0,0
    for key in position.keys():
        if key == 'nant': continue
        if position[key]['top_x'] < xmin: xmin = position[key]['top_x']
        if position[key]['top_y'] < ymin: ymin = position[key]['top_y']
    for ant in range(0,nant):
        try:
            x = position[ant]['top_x'] - xmin + 0.1
            y = position[ant]['top_y'] - ymin + 0.1
        except(KeyError): continue
        for z, pol in enumerate(pols):
            z = 2**z
            i = heracal.omni.Antpol(ant,pol,nant)
            antpos[i.val,0],antpos[i.val,1],antpos[i.val,2] = x,y,z
    reds = heracal.omni.compute_reds(nant, pols, antpos[:nant],tol=0.01)
    ex_ants = [heracal.omni.Antpol(i,nant).ant() for i in range(antpos.shape[0]) if antpos[i,0] < 0]
    kwargs['ex_ants'] = kwargs.get('ex_ants',[]) + ex_ants
    reds = heracal.omni.filter_reds(reds, **kwargs)
    if fcal:
        from heracal.firstcal import FirstCalRedundantInfo
        info = FirstCalRedundantInfo(nant)
    else:
        info = heracal.omni.RedundantInfo(nant)
    info.init_from_reds(reds, antpos)
    return info


def cal_reds_from_pos(position,**kwargs):
    nant = position['nant']
    antpos = -np.ones((nant,3))
    xmin = 0
    ymin = 0
    for key in position.keys():
        if key == 'nant': continue
        if position[key]['top_x'] < xmin: xmin = position[key]['top_x']
        if position[key]['top_y'] < ymin: ymin = position[key]['top_y']
    for ant in range(0,nant):
        try:
            x = position[ant]['top_x'] - xmin + 0.1
            y = position[ant]['top_y'] - ymin + 0.1
        except(KeyError): continue
        z = 0
        i = ant
        antpos[i,0],antpos[i,1],antpos[i,2] = x,y,z
    reds = omnical.arrayinfo.compute_reds(antpos,tol=0.01)
    kwargs['ex_ants'] = kwargs.get('ex_ants',[]) + [i for i in range(antpos.shape[0]) if antpos[i,0] < 0]
    reds = omnical.arrayinfo.filter_reds(reds,**kwargs)
    return reds


def get_phase(fqs,tau, offset=False):
    fqs = fqs.reshape(-1,1) #need the extra axis
    if offset:
        delay = tau[0]
        offset = tau[1]
        return np.exp(-1j*(2*np.pi*fqs*delay) - offset)
    else:
        return np.exp(-2j*np.pi*fqs*tau)


def save_gains_fc(s,fqs,outname):
    s2 = {}
    for k,i in s.iteritems():
        if len(i) > 1:
            s2[str(k)] = get_phase(fqs,i,offset=True).T
            s2['d'+str(k)] = i[0]
            s2['o'+str(k)] = i[1]
        else:
            s2[str(k)] = get_phase(fqs,i).T
            s2['d'+str(k)] = i
    np.savez(outname,**s2)


def load_gains_fc(fcfile):
    g0 = {}
    fc = np.load(fcfile)
    for k in fc.keys():
        if k[0].isdigit():
            a = int(k[:-1])
            p = k[-1]
            if not g0.has_key(p): g0[p] = {}
            g0[p][a] = fc[k]
    return g0


def save_gains_omni(filename, meta, gains, vismdl, xtalk):
    d = {}
    metakeys = ['jds','lsts','freqs','history']
    for key in meta:
        if key.startswith('chisq'): d[key] = meta[key] #separate if statements  pending changes to chisqs
        for k in metakeys:
            if key.startswith(k): d[key] = meta[key]
    for pol in gains:
        for ant in gains[pol]:
            d['%d%s' % (ant,pol)] = gains[pol][ant]
    for pol in vismdl:
        for bl in vismdl[pol]:
            d['<%d,%d> %s' % (bl[0],bl[1],pol)] = vismdl[pol][bl]
    for pol in xtalk:
        for bl in xtalk[pol]:
            d['(%d,%d) %s' % (bl[0],bl[1],pol)] = xtalk[pol][bl]
    np.savez(filename,**d)


def load_gains_omni(filename):
    meta, gains, vismdl, xtalk = {}, {}, {}, {}
    def parse_key(k):
        bl,pol = k.split()
        bl = tuple(map(int,bl[1:-1].split(',')))
        return pol,bl
    npz = np.load(filename)
    for k in npz.files:
        if k[0].isdigit():
            pol,ant = k[-1:],int(k[:-1])
            if not gains.has_key(pol): gains[pol] = {}
            gains[pol][ant] = npz[k]
        try: pol,bl = parse_key(k)
        except(ValueError): continue
        if k.startswith('<'):
            if not vismdl.has_key(pol): vismdl[pol] = {}
            vismdl[pol][bl] = npz[k]
        elif k.startswith('('):
            if not xtalk.has_key(pol): xtalk[pol] = {}
            xtalk[pol][bl] = npz[k]
        kws = ['chi','hist','j','l','f']
        for kw in kws:
            for k in [f for f in npz.files if f.startswith(kw)]: meta[k] = npz[k]
    return meta, gains, vismdl, xtalk


def quick_load_gains(filename):
    d = np.load(filename)
    gains = {}
    for k in d.keys():
        if k[0].isdigit():
            p = k[-1]
            if not gains.has_key(p): gains[p] = {}
            a = int(k[:-1])
            gains[p][a] = d[k]
    return gains


def load_gains_fhd(fhdsav):
    fhd_cal = readsav(fhdsav,python_dict=True)
    gfhd = {'x':{},'y':{}}
    for a in range(fhd_cal['cal']['N_TILE'][0]):
        gfhd['x'][a] = fhd_cal['cal']['GAIN'][0][0][a]
        gfhd['y'][a] = fhd_cal['cal']['GAIN'][0][1][a]
    return gfhd


def fill_flags(data,flag,fit_order = 4):
    dout = np.copy(data)
    wgt = np.logical_not(flag)
    SH = data.shape
    time_stack = np.sum(wgt,axis=1)
    for ii in range(SH[0]):
        if time_stack[ii] <= (SH[1]/2 + 1) : continue
        for jj in range(24):
            chunk = np.arange(16*jj+1,16*jj+15)
            ind = np.where(wgt[ii][chunk])
            if ind[0].size == 14: continue
            x = chunk[ind]
            y = dout[ii][chunk][ind]
            z1 = np.polyfit(x,y.real,fit_order)
            z2 = np.polyfit(x,y.imag,fit_order)
            zeros = np.where(flag[ii][chunk])
            d_temp = dout[ii][chunk]
            d_temp[zeros] = (polyfunc(chunk,z1) + 1j*polyfunc(chunk,z2))[zeros]
            dout[ii][chunk] = d_temp
    return dout


def fit_data(data,fit_order=2):
    if data.ndim == 2: d = np.mean(data,axis=0)
    else: d = data
#    fq = np.arange(d.size)
#    zr = np.polyfit(fq,d.real,fit_order)
#    zi = np.polyfit(fq,d.imag,fit_order)
#    fit_data = polyfunc(fq,zr) + 1j*polyfunc(fq,zi)
    fit_data = np.zeros(d.shape,dtype=np.complex64)
    for ii in range(24):
        chunk = np.arange(16*ii+1,16*ii+15)
        dr = d.real[chunk]
        di = d.imag[chunk]
        zr = np.polyfit(chunk,dr,fit_order)
        zi = np.polyfit(chunk,di,fit_order)
        fit_data[chunk] = polyfunc(chunk,zr)+1j*polyfunc(chunk,zi)
    return fit_data



def rough_cal(data,info,pol='xx'): #The data has to be the averaged over time axis
    p = pol[0]
    g0 = {p: {}}
    phi = {}
    reds = info.get_reds()
    reds[0].sort()
    reds[1].sort()
    redbls = reds[0] + reds[1]
    redbls.sort()
    SH = data[reds[0][0]][pol].shape
    gamma0 = fit_data(data[reds[0][0]][pol])
    gamma1 = fit_data(data[reds[1][0]][pol])
    subsetant = info.subsetant
    fixants = (min(subsetant), min(subsetant[np.where(subsetant>92)]))
    for a in fixants: phi[a] = np.zeros(SH)
    while len(redbls) > 0:
        i,j = redbls[0]
        r = (i,j)
        redbls.remove(r)
        if phi.has_key(i) and phi.has_key(j): continue
        elif phi.has_key(i) and not phi.has_key(j):
            if r in reds[0]:
                phi[j] = np.angle(fit_data(data[r][pol])*np.exp(1j*phi[i])*gamma0.conj())
            elif r in reds[1]:
                phi[j] = np.angle(fit_data(data[r][pol])*np.exp(1j*phi[i])*gamma1.conj())
        elif phi.has_key(j) and not phi.has_key(i):
            if r in reds[0]:
                phi[i] = np.angle(fit_data(data[r][pol]).conj()*np.exp(1j*phi[j])*gamma0)
            elif r in reds[1]:
                phi[i] = np.angle(fit_data(data[r][pol]).conj()*np.exp(1j*phi[j])*gamma1)
        else: redbls.append(r)
    if len(phi.keys()) != subsetant.size: raise IOError('Missing antennas')
    for a in phi.keys():
        g0[p][a] = np.exp(-1j*phi[a])
    return g0


def run_omnical(data, info, gains0=None, xtalk=None, maxiter=500, conv=1e-3,
                     stepsize=.3, trust_period=1):
    
    m1,g1,v1 = omnical.calib.logcal(data, info, xtalk=xtalk, gains=gains0,
                                    maxiter=maxiter, conv=conv, stepsize=stepsize,
                                    trust_period=trust_period)
    m2,g2,v2 = omnical.calib.lincal(data, info, xtalk=xtalk, gains=g1, vis=v1,
                                    maxiter=maxiter, conv=conv, stepsize=stepsize,
                                    trust_period=trust_period)

    return m2,g2,v2


def remove_degen_hex(gomni, antpos):
    g2 = copy.deepcopy(gomni)
    for p in g2.keys():
        ref_exp1 = np.exp(-1j*np.angle(g2[p][57]))
        ref_exp2 = np.exp(-1j*np.angle(g2[p][93]))
        for a in g2[p].keys():
            if a < 93: g2[p][a] *= ref_exp1
            else: g2[p][a] *= ref_exp2
        phi58 = g2[p][58]
        phi61 = g2[p][61]
        phi1 = np.angle(phi58)
        phi2 = np.angle(phi61)
        for a in g2[p].keys():
            if a < 93:
                dx = antpos[a]['top_x'] - antpos[57]['top_x']
                dy = antpos[a]['top_y'] - antpos[57]['top_y']
            else:
                dx = antpos[a]['top_x'] - antpos[93]['top_x']
                dy = antpos[a]['top_y'] - antpos[93]['top_y']
            nx = dx/14.-dy/np.sqrt(3)/14.
            ny = -2*dy/np.sqrt(3)/14.
            g2[p][a] *= np.exp(-1j*(phi1*nx+phi2*ny))
    g2 = scale_gains(g2)
    return g2

