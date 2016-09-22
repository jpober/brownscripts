import numpy as np
import subprocess, datetime, os
from astropy.io import fits

def writefits(npzfiles, repopath, ex_ants=[], name_dict={}):
    ### This function writes the solution from output npz files from omni_run to a fits file.   ### 
    ### npzfiles can be a list of npz files with solutions for different polarizations but for  ###.
    ### the same obs id. repopath is for writing the program of origin and git hash, e.g., if   ###
    ### the solution comes from capo, then repopath=/path/to/capo/. ex_ants is used to indicate ###
    ### which antennas are flagged. name_dict is for writing antenna names to fits.             ### 
     
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
        try: nam.append(str(name_dict[tot[nn]]))
        except(KeyError): nam.append('ant'+str(tot[nn]))      #keep this for now, may change in the future
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
    
    
def read_fits(filename, pols):
    ### This function reads in the solution from fits file, whitch returns a dictionary of polarization, ###
    ### each polarization is a dictionary of antenna indexes, which has a value as an numpy array with   ###
    ### the shape (Ntimes, Nfreqs)                                                                       ###
    g0 = {}
    poldict = {'EE': 'xx', 'NN': 'yy', 'EN': 'xy', 'NE': 'yx'}
    hdu = fits.open(filename)
    Ntimes = hdu[0].header['NTIMES']
    Nfreqs = hdu[0].header['NFREQS']
    Npols = hdu[0].header['NPOLS']
    Nants = hdu[0].header['NANTS']
    ant_index = hdu[1].data['ANT INDEX'][0:Nants]
    pol_list = hdu[1].data['POL'][0:Nfreqs*Nants*Npols].reshape(Npols,Nants*Nfreqs)[:,0]
    data_list = hdu[1].data['GAIN'].reshape((Ntimes,Npols,Nfreqs,Nants)).swapaxes(0,1).swapaxes(2,3).swapaxes(1,2) #Npols,Nants,Ntimes,Nfreqs
    for ii in range(0,Npols):
        polarization = poldict[pol_list[ii]]
        if not polarization in pols: continue
        g0[polarization[0]] = {}
        for jj in range(0,Nants):
            g0[polarization[0]][ant_index[jj]]=data_list[ii][jj]
    return g0

