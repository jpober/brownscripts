import numpy as np, omnical, aipy, math
import uvdata.uvdata as uvd
import subprocess, datetime, os
from astropy.io import fits

def writefits(npzfiles, repopath, ex_ants):
    
    p2pol = {'EE': 'x','NN': 'y','EN': 'cross', 'NE': 'cross'}

    fn0 = npzfiles[0].split('.')
    if len(npzfiles) > 1: fn0.remove[fn0[-2]]
    fn0[-1] = 'fits'
    outfn = '.'.join(fn0)
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
    tot = ant + ex_ants
    tot.sort()
    time = data['jds']
    freq = data['freqs']/1e6
    pol = ['EE', 'NN', 'EN', 'NE']
    nt = time.shape[0]
    nf = freq.shape[0]
    na = len(tot)
    nam = []
    for nn in range(0,na):
        nam.append('ant'+str(tot[nn]))      #keep this for now, may change in the future
    datarray = []
    flgarray = []
    for ii in range(0,4):
        dd = []
        fl = []
        for jj in range(0,na):
            try: dd.append(datadict[str(tot[jj])+p2pol[pol[ii]]])
            except(KeyError): dd.append(np.ones((nt,nf)))
            if tot[jj] in ex_ants: fl.append(np.ones((nt,nf),dtype=int))
            else: fl.append(np.zeros((nt,nf),dtype=int))
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


def writetxt(npzfiles, repopath, ex_ants):
    
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
    tot = ant + ex_ants
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
                    outfile.write("ant%d, %d, %f, %s, %.8f, %.8f, %.8f, %d\n"%(aa,aa,df,dp,dt,da.real,da.imag,dfl))
    outfile.close()

def uv_read(filenames, filetype=None, polstr=None,antstr='cross',recast_as_array=True):
    info = {'lsts':[], 'times':[]}
    dat, flg = {},{}
    ginfo = [0,0,0]
    freqarr = []
    #    uvdata=uvd.UVData()
    if type(filenames) == str: filenames = [filenames]
    for filename in filenames:
        uvdata = uvd.UVData()
        if filetype == 'miriad':
            uvdata.read_miriad(filename)
        elif filetype == 'uvfits':
            uvdata.read_uvfits(filename)
        elif filetype == 'fhd':                #in this case filename should be a list of files
            uvdata.read_fhd(filename)
        else:
            raise IOError('invalid filetype, it should be miriad, uvfits, or fhd')
        tt = uvdata.time_array
        Nt = uvdata.Ntimes
        blt = len(tt)
        nbl = uvdata.Nbls
        nfreq = uvdata.Nfreqs

        uvdata.set_lsts_from_time_array()
        info['times'] = uvdata.time_array[::nbl]
        info['lsts'] = uvdata.lst_array[::nbl]
        pol = uvdata.polarization_array
        npol = len(pol)
        data = uvdata.data_array
        flag = uvdata.flag_array
        ant1 = uvdata.ant_1_array
        ant2 = uvdata.ant_2_array
        
#        if not (0 in ant1 or 0 in ant2):          #if the index starts from 1
#            ones = np.ones((len(ant1)))
#            ant1 = ant1 - ones
#            ant2 = ant2 - ones

        freqarr = uvdata.freq_array[0]
        auto = 0
        
        dindex = ant1 - ant2
        if 1 in dindex and -1 in dindex: #if both (i,j) and (j,i) are included, use -1 to flag (j,i) (if j>i)
            for ii in range(0,nbl):
                if ant1[ii] > ant2[ii]:
                    ant1[ii]=-1
                    ant2[ii]=-1
        for ii in range(0, nbl):
            if ant1[ii] == ant2[ii]:
                auto += 1
        nbl -= auto
        
        nant = int((1+math.sqrt(1+8*nbl))/2)

        datcut = []
        flgcut = []

        for jj in range(0,npol):
            datcut.append(data[:,0][:,:,jj].reshape(uvdata.Ntimes,uvdata.Nbls,uvdata.Nfreqs))
            flgcut.append(flag[:,0][:,:,jj].reshape(uvdata.Ntimes,uvdata.Nbls,uvdata.Nfreqs))
    
        for ii in range(0,uvdata.Nbls):
            if ant1[ii] < 0: continue
            if ant1[ii] == ant2[ii] and antstr == 'cross': continue
            bl = (ant1[ii],ant2[ii])
            if not dat.has_key(bl): dat[bl],flg[bl] = {},{}
            for jj in range(0,npol):
                pp = aipy.miriad.pol2str[pol[jj]]
                if not dat[bl].has_key(pp):
                    dat[bl][pp],flg[bl][pp] = [],[]
                dat[bl][pp] = np.complex64(datcut[jj][:,ii])
                flg[bl][pp] = np.array(flgcut[jj][:,ii])
#                dat[bl][pp].append(data[:,0][:,:,jj][ii])
#                flg[bl][pp].append(flag[:,0][:,:,jj][ii])
        ginfo[0] = nant
        ginfo[1] = Nt
        ginfo[2] = nfreq
    return info, dat, flg, ginfo, freqarr

def uv_read_v2(filenames, filetype=None, antstr='cross'):
    info = {'lsts':[], 'times':[]}
    dat, flg = {},{}
    ginfo = [0,0,0]
    freqarr = []
    infodict = {}
    #    uvdata=uvd.UVData()
    if type(filenames) == str: filenames = [filenames]
    for filename in filenames:
        uvdata = uvd.UVData()
        if filetype == 'miriad':
            uvdata.read_miriad(filename)
        elif filetype == 'uvfits':
            uvdata.read_uvfits(filename)
        elif filetype == 'fhd':                #in this case filename should be a list of files
            uvdata.read_fhd(filename)
        else:
            raise IOError('invalid filetype, it should be miriad, uvfits, or fhd')
        tt = uvdata.time_array
        Nt = uvdata.Ntimes
        blt = len(tt)
        nbl = uvdata.Nbls
        nfreq = uvdata.Nfreqs
        
        uvdata.set_lsts_from_time_array()
        info['times'] = uvdata.time_array[::nbl]
        info['lsts'] = uvdata.lst_array[::nbl]
        pol = uvdata.polarization_array
        npol = len(pol)
        data = uvdata.data_array
        flag = uvdata.flag_array
        ant1 = uvdata.ant_1_array
        ant2 = uvdata.ant_2_array
        
        #        if not (0 in ant1 or 0 in ant2):          #if the index starts from 1
        #            ones = np.ones((len(ant1)))
        #            ant1 = ant1 - ones
        #            ant2 = ant2 - ones
        
        freqarr = uvdata.freq_array[0]
        auto = 0
        
        dindex = ant1 - ant2
        if 1 in dindex and -1 in dindex: #if both (i,j) and (j,i) are included, use -1 to flag (j,i) (if j>i)
            for ii in range(0,nbl):
                if ant1[ii] > ant2[ii]:
                    ant1[ii]=-1
                    ant2[ii]=-1
        for ii in range(0, nbl):
            if ant1[ii] == ant2[ii]:
                auto += 1
        nbl -= auto
        
        nant = int((1+math.sqrt(1+8*nbl))/2)
        
        datcut = []
        flgcut = []

        for jj in range(0,npol):
            datcut.append(data[:,0][:,:,jj].reshape(uvdata.Ntimes,uvdata.Nbls,uvdata.Nfreqs))
            flgcut.append(flag[:,0][:,:,jj].reshape(uvdata.Ntimes,uvdata.Nbls,uvdata.Nfreqs))

        for jj in range(0,npol):
            pp = aipy.miriad.pol2str[pol[jj]]
            infodict[pp] = {}
            for ii in range(0,uvdata.Nbls):
                if ant1[ii] < 0: continue
                if ant1[ii] == ant2[ii] and antstr == 'cross': continue
                bl = (ant1[ii],ant2[ii])
                if not dat.has_key(bl): dat[bl],flg[bl] = {},{}
                if not dat[bl].has_key(pp):
                    dat[bl][pp],flg[bl][pp] = [],[]
                dat[bl][pp] = np.complex64(datcut[jj][:,ii])
                flg[bl][pp] = np.array(flgcut[jj][:,ii])
            ginfo[0] = nant
            ginfo[1] = Nt
            ginfo[2] = nfreq
            infodict[pp]['timeinfo'] = info
            infodict[pp]['data'] = dat
            infodict[pp]['flag'] = flg
            infodict[pp]['ginfo'] = ginfo
            infodict[pp]['freqs'] = freqarr
            infodict[pp]['pol'] = pp
    return infodict


