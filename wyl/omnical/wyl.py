import numpy as np, omnical, aipy, math
import uvdata.uv as uvd
import subprocess, datetime

def writetxt(npzfiles, repopath, ex_ants):
    
    p2pol = {'EE': 'x','NN': 'y','EN': 'cross', 'NE': 'cross'}  #check the convension
    
    #create output file
    fn0 = npzfiles[0].split('.')
    fn0[-1] = 'txt'
    outfn = '.'.join(fn0)
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
    
    #npzdict = {}
    for f,filename in enumerate(npzfiles):
        data = np.load(filename)
        ant = []
        for ii, ss in enumerate(data):
            if ss[0].isdigit():
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
                        try: da = data[stkey][tt][ff]
                        except(KeyError): da = 1.0
                        outfile.write("ant%d, %d, %f, %s, %.8f, %f, %f, %d\n"%(aa,aa,df,dp,dt,da.real,da.imag,dfl))
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
#        timeinfo = []
#        lstsinfo = []

#        for ii in range(0,Nt):
#            timeinfo.append(tt[ii*nbl])
#            lstsinfo.append(tt[ii*nbl])   #not sure how to calculate lsts
        info['times'] = uvdata.time_array[::nbl]
        info['lsts'] = uvdata.time_array[::nbl]
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
        
#        dindex = ant1 - ant2
#        if 1 in dindex and -1 in dindex: #if both (i,j) and (j,i) are included, use -1 to flag (j,i) (if j>i)
#            for ii in range(0,blt):
#                if ant1[ii] > ant2[ii]:
#                    ant1[ii]=-1
#                    ant2[ii]=-1
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
            #if ant1[ii] < 0: continue
            if ant1[ii] == ant2[ii] and antstr == 'cross': continue
            bl = (ant1[ii],ant2[ii])
            if not dat.has_key(bl): dat[bl],flg[bl] = {},{}
            for jj in range(0,npol):
                pp = aipy.miriad.pol2str[pol[jj]]
                if not dat[bl].has_key(pp):
                    dat[bl][pp],flg[bl][pp] = [],[]
                dat[bl][pp].append(datcut[jj][:,ii])
                flg[bl][pp].append(flgcut[jj][:,ii])
#                dat[bl][pp].append(data[:,0][:,:,jj][ii])
#                flg[bl][pp].append(flag[:,0][:,:,jj][ii])
        #ginfo = [nant, Nt, nfreq]
        ginfo[0] = nant
        ginfo[1] = Nt
        ginfo[2] = nfreq
    if recast_as_array:
        for ii in dat.keys():
            for jj in dat[ii].keys():
                dat[ii][jj] = np.complex64(dat[ii][jj])
                flg[ii][jj] = np.array(flg[ii][jj])
        info['lsts'] = np.array(info['lsts'])
        info['times'] = np.array(info['times'])

    return info, dat, flg, ginfo, freqarr

