import numpy as np, omnical, aipy
import pyuvdata.uvdata as uvd
import subprocess, datetime, os
from astropy.io import fits
from uv_data_only import *
import copy

#def output_mask_array(filename, filetype, flag_array):
def output_mask_array(flag_array):
#    outfn = ''
#    if filetype == 'fhd':
#        for fn in filename:
#            if fn.endswith('_params.sav'): outfn = ('_').join(fn.split('_')[0:-1]) + '_mask.npy'
#    ### there must be a params file, otherwise an error is already raised earlier ###
#    elif filetype == 'uvfits':
#        outfn = '.'.join(filename.split('.')[0:-1]) + '_mask.npy'
#    elif filetype == 'miriad':
#        outfn = '.'.join(filename.split('.')[0:-2]) + '_mask.npy'
    invf = 1 - flag_array
    sf = np.sum((np.sum(invf,axis=0)),axis=0).astype(bool)
    st = np.sum((np.sum(invf,axis=1)),axis=1).astype(bool)
    mask_array = 1 - np.outer(st,sf)
    mask_array = mask_array.astype(bool)
#    np.save(outfn,mask_array)
    return mask_array


def writefits(npzfiles, repopath, ex_ants=[], name_dict={}):
    
    p2pol = {'EE': 'x','NN': 'y','EN': 'cross', 'NE': 'cross'}

    fn0 = npzfiles[0].split('.')
    if len(npzfiles) > 1: fn0[-2] = 'O'
    else: fn0[-2] += 'O'
    fn0[-1] = 'fits'
    outfn = '.'.join(fn0)
    print outfn
    if os.path.exists(outfn):
        print '   %s exists, skipping...' % outfn
        return 0
    githash = subprocess.check_output(['git','rev-parse','HEAD'], cwd=repopath)
    today = datetime.date.today().strftime("Date: %d, %b %Y")
    ori = subprocess.check_output(['git','remote','show','origin'], cwd=repopath)
    ori = ori.split('\n')[1].split(' ')[-1]
    githash = githash.replace('\n','')

    datadict = {}
    ant = []
    for f,filename in enumerate(npzfiles):
        data = np.load(filename)
        for ii, ss in enumerate(data):
            if ss[0].isdigit():
                datadict[ss] = data[ss]
                intss = int(ss[0:-1])
                if not intss in ant:
                    ant.append(intss)
    ant.sort()
    if name_dict == {}: tot = ant + ex_ants
    else: tot = name_dict.keys()
    tot.sort()
    time = data['jds']
    freq = data['freqs']/1e6
    pol = ['EE', 'NN', 'EN', 'NE']
    nt = time.shape[0]
    nf = freq.shape[0]
    na = len(tot)
    nam = []
    for nn in range(0,na):
        try: nam.append(name_dict[tot[nn]])
        except(KeyError): nam.append('ant'+str(tot[nn]))
    datarray = []
    flgarray = []
    for ii in range(0,4):
        dd = []
        fl = []
        for jj in range(0,na):
            try: dd.append(datadict[str(tot[jj])+p2pol[pol[ii]]])
            except(KeyError): dd.append(np.ones((nt,nf)))
            if tot[jj] in ex_ants: fl.append(np.ones((nt,nf),dtype=bool))
            else: fl.append(np.zeros((nt,nf),dtype=bool))
        datarray.append(dd)
        flgarray.append(fl)
    datarray = np.array(datarray)
    datarray = datarray.swapaxes(0,2).swapaxes(1,2).swapaxes(2,3).reshape(4*nt*nf*na)
    flgarray = np.array(flgarray)
    flgarray = flgarray.swapaxes(0,2).swapaxes(1,2).swapaxes(2,3).reshape(4*nt*nf*na)
    tarray = np.resize(time,(4*nf*na,nt)).transpose().reshape(4*nf*nt*na)
    parray = np.array((['EE']*(nf*na)+['NN']*(nf*na)+['EN']*(nf*na)+['NE']*(nf*na))*nt)
    farray = np.array(list(np.resize(freq,(na,nf)).transpose().reshape(na*nf))*4*nt)
    numarray = np.array(tot*4*nt*nf)
    namarray = np.array(nam*4*nt*nf)

    prihdr = fits.Header()
    prihdr['DATE'] = today
    prihdr['ORIGIN'] = ori
    prihdr['HASH'] = githash
    prihdr['PROTOCOL'] = 'Divide uncalibrated data by these gains to obtain calibrated data.'
    prihdr['NTIMES'] = nt
    prihdr['NFREQS'] = nf
    prihdr['NANTS'] = na
    prihdr['NPOLS'] = 4
    prihdu = fits.PrimaryHDU(header=prihdr)
    colnam = fits.Column(name='ANT NAME', format='A10', array=namarray)
    colnum = fits.Column(name='ANT INDEX', format='I',array=numarray)
    colf = fits.Column(name='FREQ (MHZ)', format='E', array=farray)
    colp = fits.Column(name='POL', format='A4', array=parray)
    colt = fits.Column(name='TIME (JD)', format='D', array=tarray)
    coldat = fits.Column(name='GAIN', format='M', array=datarray)
    colflg = fits.Column(name='FLAG', format='L', array=flgarray)
    cols = fits.ColDefs([colnam, colnum, colf, colp, colt, coldat, colflg])
    tbhdu = fits.BinTableHDU.from_columns(cols)
    hdulist = fits.HDUList([prihdu, tbhdu])
    hdulist.writeto(outfn)


def writetxt(npzfiles, repopath, ex_ants=[], name_dict={}):
    
    p2pol = {'EE': 'x','NN': 'y','EN': 'cross', 'NE': 'cross'}  #check the convension
    
    #create output file
    fn0 = npzfiles[0].split('.')
    if len(npzfiles) > 1: fn0.remove[fn0[-2]]
    fn0[-1] = 'txt'
    outfn = '.'.join(fn0)
    if os.path.exists(outfn):
        print '   %s exists, skipping...' % outfn
        return 0
    outfile = open(outfn,'w')
    githash = subprocess.check_output(['git','rev-parse','HEAD'], cwd=repopath)
    today = datetime.date.today().strftime("Date: %d, %b %Y")
    ori = subprocess.check_output(['git','remote','show','origin'], cwd=repopath)
    ori = ori.split('\n')[1].split(' ')[-1]
    outfile.write("# %s\n"%today)
    outfile.write("# Program of origin: %s\n"%ori)
    outfile.write("# Git Hash: %s"%githash)
    outfile.write("# Convention: Divide uncalibrated data by these gains to obtain calibrated data.\n")
    outfile.write("# ANT NAME, ANT INDEX, FREQ (MHZ), POL, TIME (JD), RE(GAIN), IM(GAIN), FLAG\n")
    
    #read gain solutions from npz
    
    datadict = {}
    ant = []
    for f,filename in enumerate(npzfiles):
        data = np.load(filename)
        for ii, ss in enumerate(data):
            if ss[0].isdigit():
                datadict[ss] = data[ss]
                intss = int(ss[0:-1])
                if not intss in ant:
                    ant.append(intss)
    ant.sort()
    if name_dict == {}: tot = ant + ex_ants
    else: tot = name_dict.keys()
    tot.sort()
    time = data['jds']
    freq = data['freqs']/1e6
    pol = ['EE', 'NN', 'EN', 'NE']
    nt = time.shape[0]
    nf = freq.shape[0]
    na = len(tot)
    for tt in range(0, nt):
        for pp in range(0, 4):
            for ff in range(0, nf):
                for iaa in range(0, na):
                    aa = tot[iaa]
                    dfl = 0
                    if aa in ex_ants: dfl=1
                    dt = time[tt]
                    dp = pol[pp]
                    df = freq[ff]
                    stkey = str(aa) + p2pol[pol[pp]]
                    try: da = datadict[stkey][tt][ff]
                    except(KeyError): da = 1.0
                    try: antstr = name_dict[aa]
                    except(KeyError): antstr = 'ant'+str(aa)
                    outfile.write("%s, %d, %f, %s, %.8f, %.8f, %.8f, %d\n"%(antstr,aa,df,dp,dt,da.real,da.imag,dfl))
    outfile.close()

def uv_read_fc(filenames, filetype=None, bl_str=None,antstr='cross',p_list = ['xx','yy']):
    info = {'lsts':[], 'times':[]}
    dat, flg = {},{}
    ginfo = [0,0,0]
    freqarr = []
    ex_ant = []
    if type(filenames) == str: filenames = [filenames]
    for filename in filenames:
        if filetype == 'miriad':
            uvdata = data_miriad()
            uvdata.read_data_only(filename)
        elif filetype == 'uvfits':
            uvdata = data_uvfits()
            uvdata.read_data_only(filename)
        elif filetype == 'fhd':
            uvdata = data_fhd()   #in this case filename should be a list of files
            uvdata.read_data_only(filename)
        else:
            raise IOError('invalid filetype, it should be miriad, uvfits, or fhd')
        Nt = uvdata.Ntimes
        blt = uvdata.Nblts
        nbl = uvdata.Nbls
        nfreq = uvdata.Nfreqs
        
        ex_ant = find_ex_ant(uvdata)
        info['times'] = uvdata.time_array[0]
        info['lsts'] = np.array([])   #uvdata.lst_array[::nbl]
        pol = uvdata.polarization_array
        npol = len(pol)
        data = uvdata.data_array
        flag = uvdata.flag_array
        ant1 = uvdata.ant_1_array[:nbl]
        ant2 = uvdata.ant_2_array[:nbl]
        
#        if not (0 in ant1 or 0 in ant2):          #if the index starts from 1
#            ones = np.ones((len(ant1)))
#            ant1 = ant1 - ones
#            ant2 = ant2 - ones

        freqarr = uvdata.freq_array[0]
        auto = 0
        exconj = 0

        dindex = ant1 - ant2
        for ii in range(0, nbl):
            if ant1[ii] == ant2[ii]: auto += 1

        if 1 in dindex and -1 in dindex: #if both (i,j) and (j,i) are included, use -1 to flag (j,i) (if j>i)
            for ii in range(0,nbl):
                if ant1[ii] > ant2[ii]:
                    ant1[ii]=-1
                    ant2[ii]=-1
                    exconj += 1

        nbl -= (auto + exconj)
        
        nant = int((1+np.sqrt(1+8*nbl))/2)
        bl_list = bl_str.split(',')
        for ii in range(0,uvdata.Nbls):
            if ant1[ii] < 0: continue
            if ant1[ii] == ant2[ii] and antstr == 'cross': continue
            ss1 = str(int(ant1[ii]))+'_'+str(int(ant2[ii]))
            ss2 = str(int(ant2[ii]))+'_'+str(int(ant1[ii]))
            if not (ss1 in bl_list or ss2 in bl_list): continue
            bl = (ant1[ii],ant2[ii])
            if not dat.has_key(bl): dat[bl],flg[bl] = {},{}
            for jj in range(0,npol):
                pp = aipy.miriad.pol2str[pol[jj]]
                if not pp in p_list: continue
                dat_temp = data[:,0][:,:,jj].reshape(uvdata.Ntimes,uvdata.Nbls,uvdata.Nfreqs)[:,ii]
                flg_temp = flag[:,0][:,:,jj].reshape(uvdata.Ntimes,uvdata.Nbls,uvdata.Nfreqs)[:,ii]
                dat_ma = np.ma.masked_array(dat_temp, mask=flg_temp)
                dat_ma = np.mean(dat_ma,axis=0)
                dat[bl][pp] = np.complex64([dat_ma.data])
                flg[bl][pp] = np.array([dat_ma.mask])

        ginfo[0] = nant
        ginfo[1] = 1
        ginfo[2] = nfreq
    return info, dat, flg, ginfo, freqarr, ex_ant

def uv_read_omni(filenames, filetype=None, antstr='cross', p_list = ['xx','yy'], use_model=False):
    ### Now only support reading in one data file once, don't load in multiple obs ids ###
    info = {'lsts':[], 'times':[]}
    ginfo = [0,0,0]
    freqarr = []
    infodict = {}
    #    uvdata=uvd.UVData()
    if type(filenames) == str: filenames = [filenames]
    for filename in filenames:
        if filetype == 'miriad':
            uvdata = data_miriad()
            uvdata.read_data_only(filename)
        elif filetype == 'uvfits':
            uvdata = data_uvfits()
            uvdata.read_data_only(filename)
        elif filetype == 'fhd':
            uvdata = data_fhd()  #in this case filename should be a list of files
            uvdata.read_data_only(filename,use_model=use_model)
        else:
            raise IOError('invalid filetype, it should be miriad, uvfits, or fhd')
        Nt = uvdata.Ntimes
        blt = uvdata.Nblts
        nbl = uvdata.Nbls
        nfreq = uvdata.Nfreqs
        
        info['times'] = uvdata.time_array[::nbl]
        info['lsts'] = np.array([]) #uvdata.lst_array[::nbl]
        pol = uvdata.polarization_array
        npol = len(pol)
        data = uvdata.data_array
        flag = uvdata.flag_array
        ant1 = uvdata.ant_1_array[:nbl]
        ant2 = uvdata.ant_2_array[:nbl]
        
        #        if not (0 in ant1 or 0 in ant2):          #if the index starts from 1
        #            ones = np.ones((len(ant1)))
        #            ant1 = ant1 - ones
        #            ant2 = ant2 - ones
        
        freqarr = uvdata.freq_array[0]
        auto = 0
        exconj = 0
        dindex = ant1 - ant2
        for ii in range(0, nbl):
            if ant1[ii] == ant2[ii]:
                auto += 1
        if 1 in dindex and -1 in dindex: #if both (i,j) and (j,i) are included, use -1 to flag (j,i) (if j>i)
            for ii in range(0,nbl):
                if ant1[ii] > ant2[ii]:
                    ant1[ii]=-1
                    ant2[ii]=-1
                    exconj += 1
        nbl -= (auto + exconj)
        nant = int((1+np.sqrt(1+8*nbl))/2)
        ex_ant = find_ex_ant(uvdata)
        for jj in range(0,npol):
            auto_corr = {}
            pp = aipy.miriad.pol2str[pol[jj]]
            if not pp in p_list: continue
            infodict[pp] = {}
            dat, flg = {},{}
            for ii in range(0,uvdata.Nbls):
                if ant1[ii] < 0: continue
                if ant1[ii] == ant2[ii]:
                    auto_corr[ant1[ii]] = np.sqrt(np.mean((data[:,0][:,:,jj].reshape(uvdata.Ntimes,uvdata.Nbls,uvdata.Nfreqs)[:,ii].real)[1:-2],axis=0)) + 1e-10 # +1e-10 to avoid division error
                    continue
                bl = (ant1[ii],ant2[ii])
                if not dat.has_key(bl): dat[bl],flg[bl] = {},{}
                if not dat[bl].has_key(pp):
                    dat[bl][pp],flg[bl][pp] = [],[]
                dat[bl][pp] = np.complex64(data[:,0][:,:,jj].reshape(uvdata.Ntimes,uvdata.Nbls,uvdata.Nfreqs)[:,ii])
                flg[bl][pp] = np.array(flag[:,0][:,:,jj].reshape(uvdata.Ntimes,uvdata.Nbls,uvdata.Nfreqs)[:,ii])
            scale = 1e-10
            for a in auto_corr.keys():
                ave = np.nanmean(auto_corr[a])
                if np.isnan(ave): continue
                scale += ave
            scale /= len(auto_corr.keys())
            for a in auto_corr.keys():
                auto_corr[a] /= scale
            ginfo[0] = nant
            ginfo[1] = Nt
            ginfo[2] = nfreq
            infodict[pp]['timeinfo'] = info
            infodict[pp]['data'] = dat
            infodict[pp]['flag'] = flg
            infodict[pp]['ginfo'] = ginfo
            infodict[pp]['freqs'] = freqarr
            infodict[pp]['pol'] = pp
            infodict[pp]['ex_ants'] = ex_ant
            infodict[pp]['auto_corr'] = auto_corr
            infodict[pp]['mask'] = output_mask_array(flag[:,0][:,:,jj].reshape(uvdata.Ntimes,uvdata.Nbls,uvdata.Nfreqs))
        infodict['name_dict'] = {}
        for ii in range(0,uvdata.Nants_telescope):
            if not infodict['name_dict'].has_key(uvdata.antenna_numbers[ii]):
                infodict['name_dict'][uvdata.antenna_numbers[ii]] = uvdata.antenna_names[ii]
#        if output_mask:
#            output_mask_array(filename, filetype, flag[:,0][:,:,0].reshape(uvdata.Ntimes,uvdata.Nbls,uvdata.Nfreqs))
    return infodict


def polyfunc(x,z):
    sum = np.zeros((x.size))
    for ii in range(z.size):
        sum *= x
        sum += z[ii]
    return sum


def mwa_bandpass_fit(gains, tile_info, amp_order=2, phs_order=1, band = 'high'):
    if band.lower() == 'high':
        fqs = np.linspace(167.075,197.715,384)
    elif band.lower() == 'low':
        fqs = np.linspace(138.995,169.635,384)
    for p in gains.keys():
        bandpass = {}
        for ant in gains[p].keys():
            cable = tile_info[ant]['cable']
            if not bandpass.has_key(cable): bandpass[cable] = {}
            bandpass[cable][ant] = gains[p][ant]
#            bandpass[cable][ant] = np.mean(gains[p][ant][1:53],axis=0)
            SH = gains[p][ant].shape
        global_bp = {}
        freq = np.arange(384)
        fuse = []
        for ii in range(384):
            if not ii%16 in [0,15]: fuse.append(ii)
        for length in bandpass.keys():
            amp = []
            for ant in bandpass[length].keys():
                normbp = np.abs(bandpass[length][ant])/np.mean(np.abs(bandpass[length][ant][fuse]))
                amp.append(normbp)
            amp = np.array(amp)
            global_bp[length] = np.mean(amp,axis=0)
        residual = {}
        for length in global_bp.keys():
            residual[length] = {}
            for ant in bandpass[length].keys():
                residual[length][ant] = bandpass[length][ant]/global_bp[length]
        fitamp,fitphs = {},{}
        for length in residual.keys():
            for ant in residual[length].keys():
                x = np.array(fuse)
                y1 = np.abs(residual[length][ant][fuse])
                y2 = np.angle(residual[length][ant][fuse])
                y2 = np.unwrap(y2)
                z1 = np.polyfit(x,y1,amp_order)
                z2 = np.polyfit(x,y2,phs_order)
                fitamp[ant] = z1
                fitphs[ant] = z2
        for length in bandpass.keys():
            if length == 150:
                for ant in bandpass[length].keys():
                    g = global_bp[length]*polyfunc(freq,fitamp[ant])*np.exp(1j*polyfunc(freq,fitphs[ant]))
                    r = np.mean(gains[p][ant][1:53],axis=0)/g - 1
                    for ii in range(0,384):
                        if ii%16 == 0: r[ii] = r[ii+1]
                        if ii%16 ==15: r[ii] = r[ii-1]
                    tau = np.fft.fftfreq(384,(fqs[-1]-fqs[0])/384)
                    delay = np.fft.fft(r,n=384)
                    inds = np.where(abs(np.abs(tau)-1)<0.3)[0]
                    ind = np.where(np.abs(delay)==np.max(np.abs(delay[inds])))[0]
                    for ii in range(delay.size):
                        if not ii in ind: delay[ii] = 0
                    reflect = np.fft.ifft(delay,n=384) + 1
                    gains[p][ant] = np.resize(g*reflect,SH)
            else:
                for ant in bandpass[length].keys():
                    g = global_bp[length]*polyfunc(freq,fitamp[ant])*np.exp(1j*polyfunc(freq,fitphs[ant]))
                    gains[p][ant] = np.resize(g,SH)
        return gains

def poly_bandpass_fit(gains,amp_order=9, phs_order=1,instru='mwa'):
    for p in gains.keys():
        for a in gains[p].keys():
            SH = gains[p][a].shape
            g = copy.copy(gains[p][a])
#            g = np.mean(gains[p][a],axis=0)
            fqs = np.arange(g.size)
            fuse = []
            for ff in range(g.size):
                if instru=='mwa' and ff%16 in [0,15]: continue
                fuse.append(ff)
            fuse = np.array(fuse)
            z1 = np.polyfit(fuse,np.abs(g)[fuse],amp_order)
            z2 = np.polyfit(fuse,np.unwrap(np.angle(g)[fuse]),phs_order)
            gains[p][a] = polyfunc(fqs,z1)*np.exp(1j*polyfunc(fqs,z2))
            gains[p][a] = np.resize(gains[p][a],SH)
    return gains


def ampproj(v2,model_dict,realpos,reds,tave=False):
    amppar = {}
    for p in v2.keys():
        s1,s2 = 0,0
        mdata = model_dict[p]['data']
        mflag = model_dict[p]['flag']
        for bl in v2[p].keys():
            i,j = bl
            ri,rj = realpos[i],realpos[j]
            dr = np.array([ri['top_x']-rj['top_x'],ri['top_y']-rj['top_y'],ri['top_z']-rj['top_z']])
            if np.linalg.norm(dr) < (50*3e8/180e6): continue
            for r in reds:
                if bl in r or bl[::-1] in r:
                    for rbl in r:
                        try:
                            marr = np.ma.masked_array(mdata[rbl][p],mflag[rbl][p])
                        except(KeyError):
                            marr = np.ma.masked_array(mdata[rbl[::-1]][p],mflag[rbl[::-1]][p])
                        if tave:
                            marr = np.mean(marr,axis=0)
                            marr = marr.reshape(1,-1)
                        s1 += (np.abs(v2[p][bl])*np.abs(marr.data)*np.logical_not(marr.mask))
                        s2 += (np.abs(marr.data)*np.abs(marr.data)*np.logical_not(marr.mask))
        ind = np.where(s2==0)
        s2[ind] = np.inf
        A = s1/s2
        amppar[p[0]] = np.sqrt(A)
    return amppar

def phsproj(omni,fhd,realpos,EastHex,SouthHex,ref_antenna):
    phspar = {}
    ax1,ax2 = [],[]
    for ii in range(EastHex.shape[0]):
        ind_east = np.where(EastHex[ii]>0)[0]
        ind_south = np.where(SouthHex[ii]>0)[0]
        ax1.append(EastHex[ii][ind_east])
        ax1.append(SouthHex[ii][ind_south])
    for jj in range(EastHex.shape[1]):
        ind_east = np.where(EastHex[:,jj]>0)[0]
        ind_south = np.where(SouthHex[:,jj]>0)[0]
        ax2.append(EastHex[:,jj][ind_east])
        ax2.append(SouthHex[:,jj][ind_south])
    for p in omni.keys():
        phspar[p] = {}
        phspar[p]['phix'],phspar[p]['phiy'],phspar[p]['offset_east'],phspar[p]['offset_south'] = [],[],[],[]
        SH = omni[p][ref_antenna].shape
        for tt in range(0,SH[0]):
            slp1 = []
            slp2 = []
            for ff in range(0,SH[1]):
                if ff%16 in [0,15]:
                    slp1.append(0)
                    slp2.append(0)
                    continue
            #***** East-West direction fit *****#
                slope = []
                for inds in ax1:
                    x,tau = [],[]
                    for ii in inds:
                        if not ii in omni[p].keys(): continue
                        x.append(realpos[ii]['top_x'])
                        tau.append(np.angle(fhd[p][ii][ff]/omni[p][ii][tt][ff]))
                    tau = np.unwrap(tau)
                    if tau.size < 3: continue
                    z = np.polyfit(x,tau,1)
                    slope.append(z[0])
                slope = np.array(slope)
                slp1.append(np.median(slope)) # slope could be steep, choosing median would be more likely to avoid phase wrapping
            #***** 60 deg East-South direction fit *****#
                slope = []
                for inds in ax2:
                    x,tau = [],[]
                    for ii in inds:
                        if not ii in omni[p].keys(): continue
                        x.append(realpos[ii]['top_x'])
                        tau.append(np.angle(fhd[p][ii][ff]/omni[p][ii][tt][ff]))
                    tau = np.unwrap(tau)
                    if tau.size < 3: continue
                    z = np.polyfit(x,tau,1)
                    slope.append(z[0])
                slope = np.array(slope)
                slp2.append(np.median(slope))
        #****** calculate offset term ************#
            offset1, offset2 = [],[]
            phix = np.array(slp1)
            phiy = (np.array(slp2) - phix)/np.sqrt(3)
            for a in omni[p].keys():
                dx = realpos[a]['top_x'] - realpos[ref_antenna]['top_x']
                dy = realpos[a]['top_y'] - realpos[ref_antenna]['top_y']
                proj = np.exp(1j*(dx*phix+dy*phiy))
                offset = np.exp(1j*np.angle(fhd[p][a]/omni[p][a][tt]/proj))
                if a < 93: offset1.append(offset)
                else: offset2.append(offset)
            offset1 = np.array(offset1)
            offset2 = np.array(offset2)
            offset1 = np.mean(offset1,axis=0)
            offset2 = np.mean(offset2,axis=0)
            phspar[p]['phix'].append(phix)
            phspar[p]['phiy'].append(phiy)
            phspar[p]['offset_east'].append(offset1)
            phspar[p]['offset_south'].append(offset2)
        phspar[p]['phix'] = np.array(phspar[p]['phix'])
        phspar[p]['phiy'] = np.array(phspar[p]['phiy'])
        phspar[p]['offset_east'] = np.array(phspar[p]['offset_east'])
        phspar[p]['offset_south'] = np.array(phspar[p]['offset_south'])
    return phspar

def linproj(omni,fhd,realpos,maxiter=50,conv=1e-6):
    fuse = []
    for ii in range(0,384):
        if not ii%16 in [0,15]: fuse.append(ii)
    proj = {}
    for p in omni.keys():
        a0 = omni[p].keys()[0]
        SH = omni[p][a0].shape
        proj[p] = {}
        M = np.zeros((3,3))
        n = len(omni[p].keys())
        r = {}
        proj[p]['eta'] = np.zeros(SH)
        proj[p]['phix'] = np.zeros(SH)
        proj[p]['phiy'] = np.zeros(SH)
        proj[p]['offset'] = np.zeros(SH)
        for a in omni[p].keys():
            r[a] = fhd[p][a]/omni[p][a]
            x = realpos[a]['top_x']/100
            y = realpos[a]['top_y']/100
            M += np.array([[x*x,x*y,x],
                           [x*y,y*y,y],
                           [x,  y,  1]])
        invM = np.linalg.inv(M)
        for tt in range(0,SH[0]):
            for ii in range(0,maxiter):
                b = np.zeros((3,SH[1]))
                eta = 0
                for a in omni[p].keys():
                    b += np.array([x*r[a][tt].imag,y*r[a][tt].imag,r[a][tt].imag])
                    eta += (r[a][tt].real-1)
                eta /= n
                phs = invM.dot(b)
                if np.max((eta*eta)[fuse])+np.max((phs*phs)[:,fuse]) < conv:
                    print 'maxiter: ',ii
                    break
                proj[p]['eta'][tt] += eta
                proj[p]['phix'][tt] += phs[0]
                proj[p]['phiy'][tt] += phs[1]
                proj[p]['offset'][tt] += phs[2]
                for a in omni[p].keys():
                    x = realpos[a]['top_x']/100
                    y = realpos[a]['top_y']/100
                    factor = np.exp(eta+1j*(x*phs[0]+y*phs[1]+phs[2]))
                    r[a][tt] /= factor
    return proj

def non_hex_cal(data,g2,model_dict,realpos,ex_ants=[]):
    fqflag = []
    for ii in range(384):
        if ii%16 in [0,15]: fqflag.append(ii)
    g3 = {}
    for p in g2.keys():
        g3[p] = {}
        a = g2[p].keys()[0]
        SH = g2[p][a].shape
        pp = p+p
        mvis = model_dict['data']
        mwgt = model_dict['flag']
        for a1 in range(0,57):
            nur,nui,den = 0,0,0
            if a1 in ex_ants: continue
            for a2 in g2[p].keys():
                sep = np.array([realpos[a2]['top_x']-realpos[a1]['top_x'],
                                realpos[a2]['top_y']-realpos[a1]['top_y'],
                                realpos[a2]['top_z']-realpos[a1]['top_z']])
                if np.linalg.norm(sep) < 50*3e8/np.max(model_dict['freqs']): continue
                bl = (a1,a2)
                try: dv = data[bl][pp]
                except(KeyError): dv = data[bl[::-1]][pp].conj()
                try:
                    dm = mvis[bl][pp]*(g2[p][a2].conj())
                    dw = np.logical_not(mwgt[bl][pp])
                except(KeyError):
                    dm = mvis[bl[::-1]][pp].conj()*g2[p][a2]
                    dw = np.logical_not(mwgt[bl[::-1]][pp])
                nur += np.nansum((dv.real*dm.real+dv.imag*dm.imag)*dw,axis=0)
                nui += np.nansum((dv.imag*dm.real-dv.real*dm.imag)*dw,axis=0)
                den += np.nansum((dm.real*dm.real+dm.imag*dm.imag)*dw,axis=0)
            if np.nansum(den) == 0: continue
            den[fqflag] = 1
            g3[p][a1] = nur/den + 1.j*nui/den
    return g3




