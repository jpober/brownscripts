import numpy as np, omnical, aipy, math
import uvdata.uv as uvd
import subprocess, datetime

def writetxt(npzfiles, repopath):
    
    p2pol = {'EE': 'x','NN': 'y','EN': 'cross', 'NE': 'cross'}  #check the convension
    
    #create output file
    fn0 = npzfiles[0].split('.')
    fn0[-1] = 'txt'
    outfn = '.'.join(fn0)
    outfile = open(outfn,'w')
    githash = subprocess.check_output(['git','rev-parse','HEAD'], cwd=repopath)
    today = datetime.date.today().strftime("Date: %d, %b %Y")
    outfile.write("# %s\n"%today)
    outfile.write("# Program of origin: https://github.com/wenyang-li/capo.git\n")
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
        time = data['jds']
        freq = data['freqs']/1e6
        pol = ['EE', 'NN', 'EN', 'NE']
        nt = time.shape[0]
        nf = freq.shape[0]
        na = len(ant)
        for tt in range(0, nt):
            for pp in range(0, 4):
                for ff in range(0, nf):
                    for iaa in range(0, na):
                        aa = ant[iaa]
                        dt = time[tt]
                        dp = pol[pp]
                        df = freq[ff]
                        stkey = str(aa) + p2pol[pol[pp]]
                        try: da = data[stkey][tt][ff]
                        except: da = 1.0
                        outfile.write("ant%d, %d, %f, %s, %.8f, %f, %f, 0\n"%(aa,aa,df,dp,dt,da.real,da.imag))
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
        elif filetype == 'fhd':
            uvdata.read_fhd(filename)
        else:
            raise IOError('invalid filetype, it should be miriad, uvfits, or fhd')
        tt = uvdata.time_array.value
        Nt = uvdata.Ntimes.value
        blt = len(tt)
        nbl = uvdata.Nbls.value
        nfreq = uvdata.Nfreqs.value
        timeinfo = []
        lstsinfo = []
        
        for ii in range(0,Nt):
            timeinfo.append(tt[ii*nbl])
            lstsinfo.append(tt[ii*nbl])   #not sure how to calculate lsts
        info['times'] = timeinfo
        info['lsts'] = lstsinfo
        pol = uvdata.polarization_array.value
        npol = len(pol)
        data = uvdata.data_array.value
        flag = uvdata.flag_array.value
        ant1 = uvdata.ant_1_array.value
        ant2 = uvdata.ant_2_array.value
        
        if filetype == 'fhd':          #for fhd, the index starts from 1
            ones = np.ones((len(ant1)))
            ant1 = ant1 - ones
            ant2 = ant2 - ones
        
        freqarr = uvdata.freq_array.value[0]
        auto = 0
        
        for ii in range(0, nbl):
            if ant1[ii] == ant2[ii]:
                auto += 1
        nbl -= auto
        
        nant = int((1+math.sqrt(1+8*nbl))/2)
        
        for ii in range(0,blt):
            if ant1[ii] == ant2[ii] and antstr == 'cross': continue
            bl = (ant1[ii],ant2[ii])
            if not dat.has_key(bl): dat[bl],flg[bl] = {},{}
            for jj in range(0,npol):
                pp = aipy.miriad.pol2str[pol[jj]]
                if not dat[bl].has_key(pp):
                    dat[bl][pp],flg[bl][pp] = [],[]
                data00,flag00 = [],[]
                for nn in range(0,len(data[ii][0])):
                    data00.append(data[ii][0][nn][jj])
                    flag00.append(flag[ii][0][nn][jj])
                dat[bl][pp].append(data00)
                flg[bl][pp].append(flag00)
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

