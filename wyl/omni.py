import numpy as np, omnical, aipy, math
import uvdata.uv as uvd
import capo.red as red
import numpy.linalg as la
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)
    import scipy.sparse as sps
    
POL_TYPES = 'xylrab'
#XXX this can't support restarts or changing # pols between runs
POLNUM = {} # factor to multiply ant index for internal ordering, 
NUMPOL = {}

def add_pol(p):
    global NUMPOL
    assert(p in POL_TYPES)
    POLNUM[p] = len(POLNUM)
    NUMPOL = dict(zip(POLNUM.values(), POLNUM.keys()))
    
class Antpol:
    def __init__(self, *args):
        try:
            ant,pol,nant = args
            if not POLNUM.has_key(pol): add_pol(pol)
            self.val, self.nant = POLNUM[pol] * nant + ant, nant
        except(ValueError): self.val, self.nant = args
    def antpol(self): return self.val % self.nant, NUMPOL[self.val / self.nant]
    def ant(self): return self.antpol()[0]
    def pol(self): return self.antpol()[1]
    def __int__(self): return self.val
    def __hash__(self): return self.ant()
    def __str__(self): return ''.join(map(str, self.antpol()))
    def __eq__(self, v): return self.ant() == v
    def __repr__(self): return str(self)
        
## XXX filter_reds w/ pol support should probably be in omnical
def filter_reds(reds, bls=None, ex_bls=None, ants=None, ex_ants=None, ubls=None, ex_ubls=None, crosspols=None, ex_crosspols=None):
    '''Filter redundancies to include/exclude the specified bls, antennas, and unique bl groups and polarizations.
    Assumes reds indices are Antpol objects.'''
    def pol(bl): return bl[0].pol() + bl[1].pol()
    if crosspols: reds = [r for r in reds if pol(r[0]) in crosspols]
    if ex_crosspols: reds = [r for r in reds if not pol(r[0]) in ex_crosspols]
    return omnical.arrayinfo.filter_reds(reds, bls=bls, ex_bls=ex_bls, ants=ants, ex_ants=ex_ants, ubls=ubls, ex_ubls=ex_ubls)

class RedundantInfo(omnical.info.RedundantInfo):
    def __init__(self, nant, filename=None):
        omnical.info.RedundantInfo.__init__(self, filename=filename)
        self.nant = nant
    def bl_order(self):
        '''Return (i,j) baseline tuples in the order that they should appear in data.  Antenna indicies
        are in real-world order (as opposed to the internal ordering used in subsetant).'''
        return [(Antpol(self.subsetant[i],self.nant),Antpol(self.subsetant[j],self.nant)) for (i,j) in self.bl2d]
    def order_data(self, dd):
        '''Create a data array ordered for use in _omnical.redcal.  'dd' is
        a dict whose keys are (i,j) antenna tuples; antennas i,j should be ordered to reflect
        the conjugation convention of the provided data.  'dd' values are 2D arrays
        of (time,freq) data.'''
        d = []
        for i,j in self.bl_order():
            bl = (i.ant(),j.ant())
            pol = i.pol() + j.pol()
            try: d.append(dd[bl][pol])
            except(KeyError): d.append(dd[bl[::-1]][pol[::-1]].conj())
        return np.array(d).transpose((1,2,0))

class FirstCalRedundantInfo(omnical.info.FirstCalRedundantInfo):
    def __init__(self, nant):
        omnical.info.FirstCalRedundantInfo.__init__(self)
        self.nant = nant
        print 'Loading FirstCalRedundantInfo class' 

def compute_reds(nant, pols, *args, **kwargs):
    _reds = omnical.arrayinfo.compute_reds(*args, **kwargs)
    reds = []
    for pi in pols:
        for pj in pols:
            reds += [[(Antpol(i,pi,nant),Antpol(j,pj,nant)) for i,j in gp] for gp in _reds]
    return reds


def aa_to_info(aa, pols=['x'], **kwargs):
    '''Use aa.ant_layout to generate redundances based on ideal placement.
        The remaining arguments are passed to omnical.arrayinfo.filter_reds()'''
    layout = aa.ant_layout
    nant = len(aa)
    antpos = -np.ones((nant*len(pols),3)) # -1 to flag unused antennas
    xs,ys = np.indices(layout.shape)
    for ant,x,y in zip(layout.flatten(), xs.flatten(), ys.flatten()):
        for z,pol in enumerate(pols):
            z = 2**z # exponential ensures diff xpols aren't redundant w/ each other
            i = Antpol(ant,pol,len(aa)) # creates index in POLNUM/NUMPOL for pol
            antpos[i,0],antpos[i,1],antpos[i,2] = x,y,z
    reds = compute_reds(nant, pols, antpos[:nant],tol=.1) # only first nant b/c compute_reds treats pol redundancy separately
    # XXX haven't enforced xy = yx yet.  need to conjoin red groups for that
    ex_ants = [Antpol(i,nant).ant() for i in range(antpos.shape[0]) if antpos[i,0] < 0]
    kwargs['ex_ants'] = kwargs.get('ex_ants',[]) + ex_ants
    reds = filter_reds(reds, **kwargs)
    info = RedundantInfo(nant)
    info.init_from_reds(reds,antpos)
    return info

#generate info from real positions
####################################################################################################
def aa_pos_to_info(aa, pols=['x'], **kwargs):
    '''Use aa.ant_layout to generate redundances based on real placement.
        The remaining arguments are passed to omnical.arrayinfo.filter_reds()'''
    nant = len(aa)
    antpos = -np.ones((nant*len(pols),3)) # -1 to flag unused antennas
    for ant in xrange(nant):
        bl = aa.get_baseline(0,ant,src='z')
        x,y = bl[0], bl[1]#w is currently not included
        for z,pol in enumerate(pols):
            z = 2**z # exponential ensures diff xpols aren't redundant w/ each other
            i = Antpol(ant,pol,len(aa)) # creates index in POLNUM/NUMPOL for pol
            antpos[i,0],antpos[i,1],antpos[i,2] = x,y,z
    reds = compute_reds(nant, pols, antpos[:nant],tol=2) # only first nant b/c compute_reds treats pol redundancy separately
    # XXX haven't enforced xy = yx yet.  need to conjoin red groups for that
    ex_ants = [Antpol(i,nant).ant() for i in range(antpos.shape[0]) if antpos[i,0] < 0]
    kwargs['ex_ants'] = kwargs.get('ex_ants',[]) + ex_ants
    reds = filter_reds(reds, **kwargs)
    info = RedundantInfo(nant)
    info.init_from_reds(reds,antpos)
    return info
####################################################################################################


def redcal(data, info, xtalk=None, gains=None, vis=None,removedegen=False, uselogcal=True, maxiter=50, conv=1e-3, stepsize=.3, computeUBLFit=True, trust_period=1):
    #add layer to support new gains format
    if gains:
        _gains = {}
        for pol in gains:
            for i in gains[pol]:
                ai = Antpol(i,pol,info.nant)
                _gains[int(ai)] = gains[pol][i].conj()
    else: _gains = gains
    if vis:
        _vis = {}
        for pol in vis:
            for i,j in vis[pol]:
                ai,aj = Antpol(i,pol[0],info.nant), Antpol(j,pol[1],info.nant)
                _vis[(int(ai),int(aj))] = vis[pol][(i,j)]
    else: _vis = vis
    meta, gains, vis = omnical.calib.redcal(data, info, xtalk=xtalk, gains=_gains, vis=_vis, removedegen=removedegen, uselogcal=uselogcal, maxiter=maxiter, conv=conv, stepsize=stepsize, computeUBLFit=computeUBLFit, trust_period=trust_period)    
    # rewrap to new format
    def mk_ap(a): return Antpol(a, info.nant)
    for i,j in meta['res'].keys():
        api,apj = mk_ap(i),mk_ap(j)
        pol = api.pol() + apj.pol()
        bl = (api.ant(), apj.ant())
        if not meta['res'].has_key(pol): meta['res'][pol] = {}
        meta['res'][pol][bl] = meta['res'].pop((i,j))
    #XXX make chisq a nested dict, with individual antpol keys?
    for k in [k for k in meta.keys() if k.startswith('chisq')]:
        try:
            ant = int(k.split('chisq')[1])
            meta['chisq'+str(mk_ap(ant))] = meta.pop(k)
        except(ValueError): pass
    for i in gains.keys():
        ap = mk_ap(i)
        if not gains.has_key(ap.pol()): gains[ap.pol()] = {}
        gains[ap.pol()][ap.ant()] = gains.pop(i).conj()
    for i,j in vis.keys():
        api,apj = mk_ap(i),mk_ap(j)
        pol = api.pol() + apj.pol()
        bl = (api.ant(), apj.ant())
        if not vis.has_key(pol): vis[pol] = {}
        vis[pol][bl] = vis.pop((i,j))
    return meta, gains, vis

def compute_xtalk(res, wgts):
    '''Estimate xtalk as time-average of omnical residuals.'''
    xtalk = {}
    for pol in res.keys():
        xtalk[pol] = {}
        for key in res[pol]: 
            r,w = np.where(wgts[pol][key] > 0, res[pol][key], 0), wgts[pol][key].sum(axis=0)
            w = np.where(w == 0, 1, w)
        xtalk[pol][key] = (r.sum(axis=0) / w).astype(res[pol][key].dtype) # avg over time
    return xtalk

def to_npz(filename, meta, gains, vismdl, xtalk):
    '''Write results from omnical.calib.redcal (meta,gains,vismdl) to npz file.
    Each of these is assumed to be a dict keyed by pol, and then by bl/ant/keyword'''
    d = {}
    metakeys = ['jds','lsts','freqs','history']#,chisq]
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

def from_npz(filename, verbose=False):
    '''Reconstitute results from to_npz, returns meta, gains, vismdl, xtalk, each
    keyed first by polarization, and then by bl/ant/keyword.'''
    if type(filename) is str: filename = [filename]
    meta, gains, vismdl, xtalk = {}, {}, {}, {}
    def parse_key(k):
        bl,pol = k.split()
        bl = tuple(map(int,bl[1:-1].split(',')))
        return pol,bl
    for f in filename:
        if verbose: print 'Reading', f
        npz = np.load(f)
        for k in [f for f in npz.files if f.startswith('(')]:
            pol,bl = parse_key(k)
            if not xtalk.has_key(pol): xtalk[pol] = {}
            xtalk[pol][bl] = xtalk[pol].get(bl,[]) + [np.copy(npz[k])]
        for k in [f for f in npz.files if f.startswith('<')]:
            pol,bl = parse_key(k)
            if not vismdl.has_key(pol): vismdl[pol] = {}
            vismdl[pol][bl] = vismdl[pol].get(bl,[]) + [np.copy(npz[k])]
        for k in [f for f in npz.files if f[0].isdigit()]:
            pol,ant = k[-1:],int(k[:-1])
            if not gains.has_key(pol): gains[pol] = {}
            gains[pol][ant] = gains[pol].get(bl,[]) + [np.copy(npz[k])]
        kws = ['chi','hist','j','l','f']
        for kw in kws:
            for k in [f for f in npz.files if f.startswith(kw)]:
                meta[k] = meta.get(k,[]) + [np.copy(npz[k])]
    for pol in xtalk:
        for bl in xtalk[pol]: xtalk[pol][bl] = np.concatenate(xtalk[pol][bl])
    for pol in vismdl:
        for bl in vismdl[pol]: vismdl[pol][bl] = np.concatenate(vismdl[pol][bl])
    for pol in gains:
        for bl in gains[pol]: gains[pol][bl] = np.concatenate(gains[pol][bl])
    for k in meta:
        try: meta[k] = np.concatenate(meta[k])
        except(ValueError): pass
    return meta, gains, vismdl, xtalk

class FirstCal(object):
    def __init__(self, data, fqs, info):
        self.data = data
        self.fqs = fqs
        self.info = info
    def data_to_delays(self):
        '''data = dictionary of visibilities. 
           info = FirstCalRedundantInfo class
           Returns a dictionary with keys baseline pairs and values delays.'''
        self.blpair2delay = {}
        dd = self.info.order_data(self.data)
        for (bl1,bl2) in self.info.bl_pairs:
            d1 = dd[:,:,self.info.bl_index(bl1)]
            d2 = dd[:,:,self.info.bl_index(bl2)]
            delay = red.redundant_bl_cal_simple(d1,d2,self.fqs)
            self.blpair2delay[(bl1,bl2)] = delay    
        return self.blpair2delay
    def get_N(self,nblpairs):
        return np.identity(nblpairs) 
    def get_M(self):
        M = np.zeros((len(self.info.bl_pairs),1))
        blpair2delay = self.data_to_delays()
        for pair in blpair2delay:
            M[self.info.blpair_index(pair)] = blpair2delay[pair]
        return M
    def run(self):
        #make measurement matrix 
        self.M = self.get_M()
        #make noise matrix
        N = self.get_N(len(self.info.bl_pairs)) 
        self._N = np.linalg.inv(N)
        #get coefficients matrix,A
        self.A = self.info.A
        #solve for delays
        invert = np.dot(self.A.T,np.dot(self._N,self.A))
        dontinvert = np.dot(self.A.T,np.dot(self._N,self.M))
        self.xhat = np.dot(np.linalg.pinv(invert), dontinvert)
        #turn solutions into dictionary
        return dict(zip(self.info.subsetant,self.xhat))
    def get_solved_delay(self):
        solved_delays = []
        for pair in self.info.bl_pairs:
            ant_indexes = self.info.blpair2antind(pair)
            dlys = self.xhat[ant_indexes]
            solved_delays.append(dlys[0]-dlys[1]-dlys[2]+dlys[3])
        self.solved_delays = np.array(solved_delays)


def get_phase(fqs,tau):
    return np.exp(-2j*np.pi*fqs*tau)


######################################################################################################
#write npz solutions to a txt file

def writetxt(npzfiles):
    
    p2pol = {'EE': 'x','NN': 'y','EN': 'cross', 'NE': 'cross'}  #check the convension
    
    #create output file
    fn0 = npzfiles[0].split('.')
    fn0[-1] = 'txt'
    outfn = '.'.join(fn0)
    outfile = open(outfn,'w')
    outfile.write("# Program of origin: Omnical\n")
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
                        string = 'ant'+str(aa)+', '+str(aa)+', '+str(df)+', '+dp+', '+str(dt)+', '+str(da.real)+', '+str(da.imag)+', 0\n'
                        outfile.write(string)
    outfile.close()
######################################################################################################

#####################################################################################################
def uv_read(filenames, filetype=None, polstr=None,antstr=None,recast_as_array=True):
    info = {'lsts':[], 'times':[]}
    dat, flg = {},{}
    ginfo = [0,0,0]
    #    uvdata=uvd.UVData()
    if type(filenames) == str: filenames = [filenames]
    for filename in filenames:
        uvdata = uvd.UVData()
        if filetype == 'miriad':
            uvdata.read_miriad(filename)
        elif filetype == 'uvfits':
            uvdata.read_uvfits(filename)
        elif filetype == 'fhd':
            fnames = filename.split(',')
            uvdata.read_fhd(fnames)
        else:
            raise IOError('invalid filetype, it should be miriad, uvfits, or fhd')
        #uvdata.read(filename, filetype)
        tt = uvdata.time_array.value
        Nt = uvdata.Ntimes.value
        blt = len(tt)
        nbl = uvdata.Nbls.value
        nfreq = uvdata.Nfreqs.value
        
        for ii in range(0,Nt):
            info['times'].append(tt[ii*nbl])
            info['lsts'].append(tt[ii*nbl])   #not sure how to calculate lsts
        pol = uvdata.polarization_array.value
        npol = len(pol)
        data = uvdata.data_array.value
        flag = uvdata.flag_array.value
        ant1 = uvdata.ant_1_array.value
        ant2 = uvdata.ant_2_array.value
        freqarr = uvdata.freq_array.value[0]
        
        nant = int((1+math.sqrt(1+8*nbl))/2)
        
        #ginfo=[nant, Nt, nfreq]
        ginfo[0] = nant
        ginfo[1] = Nt
        ginfo[2] = nfreq
        
        for ii in range(0,blt):
            bl = (ant1[ii],ant2[ii])
            if antstr == 'cross':
                if ant1[ii] == ant2[ii]: continue
            if not dat.has_key(bl): dat[bl],flg[bl] = {},{}
            for jj in range(0,npol):
                pp = aipy.miriad.pol2str[pol[jj]]
                if polstr != None:
                    if pp != polstr: continue
                if not dat[bl].has_key(pp):
                    dat[bl][pp],flg[bl][pp] = [],[]
                data00,flag00 = [],[]
                for nn in range(0,len(data[ii][0])):
                    data00.append(data[ii][0][nn][jj])
                    flag00.append(flag[ii][0][nn][jj])
                dat[bl][pp].append(data00)
                flg[bl][pp].append(flag00)
#        if filetype == 'fhd': break
    if recast_as_array:
        for ii in dat.keys():
            for jj in dat[ii].keys():
                dat[ii][jj] = np.array(dat[ii][jj])
                flg[ii][jj] = np.array(flg[ii][jj])
        info['lsts'] = np.array(info['lsts'])
        info['times'] = np.array(info['times'])
    return info, dat, flg, ginfo, freqarr
#####################################################################################################



