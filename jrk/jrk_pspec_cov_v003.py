#! /usr/bin/env python

import aipy as a, numpy as n, capo
from matplotlib import pylab as p

import aipy as a, numpy as n, pylab as p, capo, capo.frf_conv as fringe

import glob, optparse, sys, random
from scipy.linalg import fractional_matrix_power
import capo.zsa as zsa

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True, cal=True)
o.add_option('-b', '--nboot', type='int', default=20,
    help='Number of bootstraps.  Default is 20.')
o.add_option('--plot', action='store_true',
    help='Generate plots')
o.add_option('--window', dest='window', default='blackman-harris',
    help='Windowing function to use in delay transform.  Default is blackman-harris.  Options are: ' + ', '.join(a.dsp.WINDOW_FUNC.keys()))
o.add_option('--sep', default='sep0,1', action='store',
    help='Which separation directory to use for signal loss data.')
o.add_option('--loss', action='store', 
    help='Uses signal loss mode to measure the signal loss. Uses default data in my path. Give it the path to the simulated signal data.')
o.add_option('--level', type='float', default=-1.0,
    help='Scalar by which to multiply the default signal level for simulation runs.')
o.add_option('--rmbls', dest='rmbls',type='string', 
    help='List of baselines (ex:1_4,2_33) to remove from the power spectrum analysis.')
o.add_option('--de', dest='dirtyeven', action='store',
    help='Dirty even LST')
o.add_option('--do', dest='dirtyodd', action='store',
    help='Dirty odd LST')
o.add_option('--re', dest='residualeven', action='store',
    help='Residual even LST')
o.add_option('--ro', dest='residualodd', action='store',
    help='Residual odd LST')

o.add_option('--noise_only',action='store_true',
    help='Instead of injecting noise, Replace data with noise')
o.add_option('--output', type='string', default='',
    help='Output directory for pspec_boot files (default "")')

opts,args = o.parse_args(sys.argv[1:])

#Basic parameters
random.seed(0)
POL = opts.pol #'xx'
LST_STATS = True
DELAY = False
NGPS = 5 #number of groups to break the random sampled bls into
INJECT_SIG = 0.
SAMPLE_WITH_REPLACEMENT = True
PLOT = opts.plot

#Remove baselines if specified
try:
    rmbls = []
    rmbls_list = opts.rmbls.split(',')
    for bl in rmbls_list:
        i,j = bl.split('_')
        rmbls.append(a.miriad.ij2bl(int(i),int(j)))   
    print 'Removing baselines:',rmbls
    #rmbls = map(int, opts.rmbls.split(','))
except:
    rmbls = []


### FUNCTIONS ###

def frf(shape,loc=0,scale=1): #FRF NOISE
    shape = shape[1]*2,shape[0] #(2*times,freqs)
    dij = noise(shape,loc=loc,scale=scale)
    wij = n.ones(shape,dtype=bool) #XXX flags are all true (times,freqs)
    _d,_w,_,_ = fringe.apply_frf(aa,dij,wij,ij[0],ij[1],pol=POL,bins=bins,firs=fir)
    _d = n.transpose(_d)
    _d = _d[:,shape[0]/4:shape[0]/2+shape[0]/4]
    return _d

def get_data(filenames, antstr, polstr, rmbls, verbose=False):
    # XXX could have this only pull channels of interest to save memory
    lsts, dat, flg = [], {}, {}
    if type(filenames) == 'str': filenames = [filenames]
    for filename in filenames:
        if verbose: print '   Reading', filename
        uv = a.miriad.UV(filename)
        a.scripting.uv_selector(uv, antstr, polstr)
        for (crd,t,(i,j)),d,f in uv.all(raw=True):
            bl = a.miriad.ij2bl(i,j)
            #print bl,uv['pol']
            if bl in rmbls: continue
            lst = uv['lst']
            if len(lsts) == 0 or lst != lsts[-1]: lsts.append(lst)
            if not dat.has_key(bl): dat[bl],flg[bl] = [],[]
            dat[bl].append(d)
            flg[bl].append(f)
    order_lst = n.argsort(lsts) #sometimes data is not in LST order!
    lsts = n.array(lsts)[order_lst]
    for bl in dat:
        dat[bl] = n.array(dat[bl])[order_lst]
        flg[bl] = n.array(flg[bl])[order_lst]
    return lsts, dat, flg
    
def cov(m):
    '''Because numpy.cov is stupid and casts as float.'''
    #return n.cov(m)
    X = n.array(m, ndmin=2, dtype=n.complex)
    X -= X.mean(axis=1)[(slice(None),n.newaxis)]
    N = X.shape[1]
    fact = float(N - 1) #normalization
    return (n.dot(X, X.T.conj()) / fact).squeeze()

def get_Q(mode, n_k): #encodes the fourier transform from freq to delay
    if not DELAY:
        _m = n.zeros((n_k,), dtype=n.complex)
        _m[mode] = 1. #delta function at specific delay mode
        m = n.fft.ifft(n.fft.ifftshift(_m*a.dsp.gen_window(nchan, WINDOW)))# * a.dsp.gen_window(nchan, WINDOW) #FFT it to go to freq
        Q = n.einsum('i,j', m, m.conj()) #dot it with its conjugate
        return Q
    else:
        # XXX need to have this depend on window
        Q = n.zeros_like(C)
        Q[mode,mode] = 1
        return Q


SEP = opts.sep #XXX used only for signal loss paths?

#Read even&odd data
#dsets = {
#    'even': [x for x in opts.dirtyeven if 'even' in x],
#    'odd' : [x for x in opts.dirtyodd if 'odd' in x]
#}
dsets = {
    'even': opts.dirtyeven.split(' '),
    'odd' : opts.dirtyodd.split(' ')
}
dsets2 = {
    'even': opts.residualeven.split(' '),
    'odd' : opts.residualodd.split(' ')
}
print 'PRINTING OUT THE DSET VALUES'
print dsets['even'][0]

#Get uv file info
WINDOW = opts.window
uv = a.miriad.UV(dsets['even'][0])#dsets.values()[0][0])
freqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
sdf = uv['sdf']
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])

inttime = uv['inttime']

afreqs = freqs.take(chans)
nchan = chans.size
fq = n.average(afreqs)
z = capo.pspec.f2z(fq)

aa = a.cal.get_aa(opts.cal, afreqs)
bls,conj = capo.red.group_redundant_bls(aa.ant_layout)
sep2ij, blconj, bl2sep = capo.zsa.grid2ij(aa.ant_layout)
jy2T = capo.pspec.jy2T(afreqs)
window = a.dsp.gen_window(nchan, WINDOW)
if not WINDOW == 'none': window.shape=(nchan,1)

# XXX NEED TO FIGURE OUT BW NORMALIZATION
B = sdf * afreqs.size * capo.pfb.NOISE_EQUIV_BW[WINDOW] #proper normalization
etas = n.fft.fftshift(capo.pspec.f2eta(afreqs)) #create etas (fourier dual to frequency)
kpl = etas * capo.pspec.dk_deta(z) 
scalar = 1

print 'Freq:',fq
print 'z:', z
print 'B:', B
print 'scalar:', scalar
sys.stdout.flush()

#Acquire data
antstr = 'cross'
lsts,data,flgs = {},{},{}
lsts2,data2,flgs2 = {},{},{}
days = dsets.keys()
for k in days:
    lsts[k],data[k],flgs[k] = get_data(dsets[k], antstr=antstr, polstr=POL, rmbls=rmbls, verbose=True)
    lsts2[k],data2[k],flgs2[k] = get_data(dsets2[k], antstr=antstr, polstr=POL, rmbls=rmbls, verbose=True)
    lsts
    #data has keys 'even' and 'odd'
    #inside that are baseline keys
    #inside that has shape (#lsts, #freqs)

#Get some statistics
cnt,var = n.ones_like(lsts.values()[0]), n.ones_like(lsts.values()[0])

#Align data in LST (even/odd data might have a different number of LSTs)
if True:
    lstmax = max([lsts[k][0] for k in days]) #the larger of the  initial lsts
    for k in days:
        for i in xrange(len(lsts[k])):
            #allow for small numerical differences (which shouldn't exist!)
            if lsts[k][i] >= lstmax - .001: break
        lsts[k] = lsts[k][i:]
        for bl in data[k]:
            data[k][bl],flgs[k][bl] = data[k][bl][i:],flgs[k][bl][i:]
            data2[k][bl],flgs2[k][bl] = data2[k][bl][i:],flgs2[k][bl][i:]
    j = min([len(lsts[k]) for k in days]) 
    for k in days:
        lsts[k] = lsts[k][:j]
        lsts2[k] = lsts2[k][:j]
        for bl in data[k]:
            data[k][bl],flgs[k][bl] = n.array(data[k][bl][:j]),n.array(flgs[k][bl][:j])
            data2[k][bl],flgs2[k][bl] = n.array(data2[k][bl][:j]),n.array(flgs2[k][bl][:j])
else:
    for k in days:
        for bl in data[k]:
            data[k][bl], flgs[k][bl] = n.array(data[k][bl][:]), n.array(flgs[k][bl][:])
            data2[k][bl], flgs2[k][bl] = n.array(data2[k][bl][:]), n.array(flgs2[k][bl][:])

lsts = lsts.values()[0] #same set of LST values for both even/odd data
daykey = data.keys()[0]
blkey = data[daykey].keys()[0]
ij = a.miriad.bl2ij(blkey)

#Prep FRF Stuff
bins = fringe.gen_frbins(inttime)
frp, bins = fringe.aa_to_fr_profile(aa, ij, len(afreqs)/2, bins=bins)

#Extract frequency range of data
x = {}
x2 = {}
#if opts.noise_only: NOISE = frf((len(chans),len(lsts)),loc=0,scale=1) #same noise for all bl
for k in days:
    x = {}
    x2 = {}
    f = {}
    f2 = {}
    for k in days:
        x[k] = {}
        x2[k] = {}
        f[k] = {}
        f2[k] = {}
        for bl in data[k]:
            d = data[k][bl][:,chans]# * jy2T
            flg = flgs[k][bl][:,chans]
            d2 = data2[k][bl][:,chans]
            flg2 = flgs2[k][bl][:,chans]
            if conj[bl]: 
                d = n.conj(d)
                d2 = n.conj(d2) #conjugate if necessary
            d = n.transpose(d) #now (freqs,times)
            d2 = n.transpose(d2)
            x[k][bl] = d 
            f[k][bl] = n.transpose(flg, [1,0])
            x2[k][bl] = d2
            f2[k][bl] = n.transpose(flg2, [1,0])
bls_master = x.values()[0].keys()
nbls = len(bls_master)
print 'Baselines:', nbls
##PSPEC Stuff
Q = [get_Q(i,nchan) for i in xrange(nchan)]
I,_I,_Ix = {},{},{}
C,_C,_Cx = {},{},{}
I2,_I2,_Ix2 = {},{},{}
C2,_C2,_Cx2 = {},{},{}
for k in days:
    I[k],_I[k],_Ix[k] = {},{},{}
    C[k],_C[k],_Cx[k] = {},{},{}
    I2[k],_I2[k],_Ix2[k] = {},{},{}
    C2[k],_C2[k],_Cx2[k] = {},{},{}
    for bl in bls_master:
        C[k][bl] = cov(x[k][bl])
        I[k][bl] = n.identity(C[k][bl].shape[0])
        U,S,V = n.linalg.svd(C[k][bl].conj()) #singular value decomposition
        _C[k][bl] = n.einsum('ij,j,jk', V.T, 1./S, U.T)
        _I[k][bl] = n.identity(_C[k][bl].shape[0])
        _Cx[k][bl] = n.dot(_C[k][bl], x[k][bl])
        _Ix[k][bl] = x[k][bl].copy()
        ###### x looks like it containts the DATA!!!!!!
        C2[k][bl] = cov(x2[k][bl])
        I2[k][bl] = n.identity(C2[k][bl].shape[0])
        U2,S2,V2 = n.linalg.svd(C2[k][bl].conj()) #singular value decomposition
        _C2[k][bl] = n.einsum('ij,j,jk', V2.T, 1./S2, U2.T)
        _I2[k][bl] = n.identity(_C2[k][bl].shape[0])
        _Cx2[k][bl] = n.dot(_C2[k][bl], x2[k][bl])
        _Ix2[k][bl] = x2[k][bl].copy()


#Make boots        
for boot in xrange(opts.nboot):
    print '%d / %d' % (boot+1,opts.nboot)
    bls = bls_master[:]
    if True: #shuffle and group baselines for bootstrapping
        if not SAMPLE_WITH_REPLACEMENT:
            random.shuffle(bls)
            #bls = bls[:-5] # XXX
        else: #sample with replacement (XXX could be biased)
            bls = [random.choice(bls) for bl in bls]
        gps = [bls[i::NGPS] for i in range(NGPS)]
        gps = [[random.choice(gp) for bl in gp] for gp in gps]
    else: #assign each baseline its own group
        gps = [bls[i::NGPS] for i in range(NGPS)]
    bls = [bl for gp in gps for bl in gp]
    #print '\n'.join([','.join(['%d_%d'%a.miriad.bl2ij(bl) for bl in gp]) for gp in gps])    
    _Iz,_Isum,_IsumQ = {},{},{}
    _Cz,_Csum,_CsumQ = {},{},{}
    _Iz2,_Isum2,_IsumQ2 = {},{},{}
    _Cz2,_Csum2,_CsumQ2 = {},{},{}
    print "   Getting C"
    for k in days:
        _Iz[k],_Isum[k],_IsumQ[k] = {},{},{}
        _Cz[k],_Csum[k],_CsumQ[k] = {},{},{}
        _Iz2[k],_Isum2[k],_IsumQ2[k] = {},{},{}
        _Cz2[k],_Csum2[k],_CsumQ2[k] = {},{},{}
        for i,gp in enumerate(gps): #sum things up over the groups
            _Iz[k][i] = sum([_Ix[k][bl] for bl in gp])
            _Cz[k][i] = sum([_Cx[k][bl] for bl in gp])
            _Isum[k][i] = sum([_I[k][bl] for bl in gp])
            _Csum[k][i] = sum([_C[k][bl] for bl in gp])
            _IsumQ[k][i] = {}
            _CsumQ[k][i] = {}
            _Iz2[k][i] = sum([_Ix2[k][bl] for bl in gp])
            _Cz2[k][i] = sum([_Cx2[k][bl] for bl in gp])
            _Isum2[k][i] = sum([_I2[k][bl] for bl in gp])
            _Csum2[k][i] = sum([_C2[k][bl] for bl in gp])
            _IsumQ2[k][i] = {}
            _CsumQ2[k][i] = {}
            if DELAY: #this is much faster
                _Iz[k][i] = n.fft.fftshift(n.fft.ifft(window*_Iz[k][i], axis=0), axes=0)
                _Cz[k][i] = n.fft.fftshift(n.fft.ifft(window*_Cz[k][i], axis=0), axes=0)
                #XXX need to take fft of _Csum, _Isum here
            for ch in xrange(nchan): #XXX this loop makes computation go as nchan^3
                _IsumQ[k][i][ch] = n.dot(_Isum[k][i], Q[ch])
                _CsumQ[k][i][ch] = n.dot(_Csum[k][i], Q[ch]) #C^-1 Q
                _IsumQ2[k][i][ch] = n.dot(_Isum2[k][i], Q[ch])
                _CsumQ2[k][i][ch] = n.dot(_Csum2[k][i], Q[ch])
    print "   Getting F and q"
    FI = n.zeros((nchan,nchan), dtype=n.complex)
    FI2 = n.zeros((nchan,nchan), dtype=n.complex)
    FC = n.zeros((nchan,nchan), dtype=n.complex)
    FC2 = n.zeros((nchan,nchan), dtype=n.complex)
    qI = n.zeros((nchan,_Iz.values()[0].values()[0].shape[1]), dtype=n.complex)
    qC = n.zeros((nchan,_Cz.values()[0].values()[0].shape[1]), dtype=n.complex)
    qI2 = n.zeros((nchan,_Iz.values()[0].values()[0].shape[1]), dtype=n.complex)
    qC2 = n.zeros((nchan,_Cz.values()[0].values()[0].shape[1]), dtype=n.complex)
    Q_Iz = {}
    Q_Cz = {}
    Q_Iz2 = {}
    Q_Cz2 = {}
    for cnt1,k1 in enumerate(days):
        for k2 in days[cnt1:]: #loop over even with even, even with odd, etc.
            if not Q_Iz.has_key(k2): Q_Iz[k2] = {}
            if not Q_Cz.has_key(k2): Q_Cz[k2] = {}
            for bl1 in _Cz[k1]:
                for bl2 in _Cz[k2]:
                    #if k1 == k2 and bl1 == bl2: continue #this results in a significant bias
                    if k1 == k2 or bl1 == bl2: continue 
                    #print k1, k2, bl1, bl2
                    
                    qI += n.conj(_Iz[k1][bl1]) * _Iz[k2][bl2]
                    qC += n.conj(_Cz[k1][bl1]) * _Cz[k2][bl2]
                    qI2 += n.conj(_Iz2[k1][bl1]) * _Iz2[k2][bl2]
                    qC2 += n.conj(_Cz2[k1][bl1]) * _Cz2[k2][bl2]
                    if DELAY: #by taking FFT of CsumQ above, each channel is already i,j separated
                        FI += n.conj(_IsumQ[k1][bl1]) * _IsumQ[k2][bl2]
                        FC += n.conj(_CsumQ[k1][bl1]) * _CsumQ[k2][bl2]
                    else:
                        for i in xrange(nchan):
                            for j in xrange(nchan):
                                FI[i,j] += n.einsum('ij,ji', _IsumQ[k1][bl1][i], _IsumQ[k2][bl2][j])
                                FC[i,j] += n.einsum('ij,ji', _CsumQ[k1][bl1][i], _CsumQ[k2][bl2][j]) #C^-1 Q C^-1 Q
                                FI2[i,j] += n.einsum('ij,ji', _IsumQ2[k1][bl1][i], _IsumQ2[k2][bl2][j])
                                FC2[i,j] += n.einsum('ij,ji', _CsumQ2[k1][bl1][i], _CsumQ2[k2][bl2][j])
    
    print "   Getting M"
    #Cholesky decomposition to get M
    order = n.array([10,11,9,12,8,20,0,13,7,14,6,15,5,16,4,17,3,18,2,19,1])
    iorder = n.argsort(order)
    FC_o = n.take(n.take(FC,order, axis=0), order, axis=1)
    L_o = n.linalg.cholesky(FC_o)
    U,S,V = n.linalg.svd(L_o.conj())
    MC_o = n.dot(n.transpose(V), n.dot(n.diag(1./S), n.transpose(U)))
    MC = n.take(n.take(MC_o,iorder, axis=0), iorder, axis=1)
    MI  = n.identity(nchan, dtype=n.complex128)
    
    FC_o2 = n.take(n.take(FC2,order, axis=0), order, axis=1)
    L_o2 = n.linalg.cholesky(FC_o2)
    U2,S2,V2 = n.linalg.svd(L_o2.conj())
    MC_o2 = n.dot(n.transpose(V2), n.dot(n.diag(1./S2), n.transpose(U2)))
    MC2 = n.take(n.take(MC_o2,iorder, axis=0), iorder, axis=1)
    MI2  = n.identity(nchan, dtype=n.complex128)

    print "   Getting W"
    #print 'Normalizing M/W'
    WI = n.dot(MI, FI)
    norm  = WI.sum(axis=-1); norm.shape += (1,)
    MI /= norm; WI = n.dot(MI, FI)
    WC = n.dot(MC, FC)
    norm  = WC.sum(axis=-1); norm.shape += (1,)
    MC /= norm; WC = n.dot(MC, FC)
    
    WI2 = n.dot(MI2, FI2)
    norm2  = WI2.sum(axis=-1); norm2.shape += (1,)
    MI2 /= norm2; WI2 = n.dot(MI2, FI2)
    WC2 = n.dot(MC2, FC2)
    norm2  = WC2.sum(axis=-1); norm2.shape += (1,)
    MC2 /= norm2; WC = n.dot(MC2, FC2)

    print '   Generating ps'
    pC = n.dot(MC, qC)# * scalar
    pI = n.dot(MI, qI)# * scalar 
    pC2 = n.dot(MC2, qC2)# * scalar
    pI2 = n.dot(MI2, qI2)# * scalar

    print 'pC ~ ', n.median(pC)
    print 'pI ~ ', n.median(pI)
    print 'pC2 ~ ', n.median(pC2)
    print 'pI2 ~ ', n.median(pI2)
    print 'ratio',n.abs(n.median(pI2)/n.median(pI))
    pI = n.nan_to_num(pI2/pI)#pI2*n.conj(pI)/(n.conj(pI)*pI))
    pC = n.nan_to_num(pC2/pC)#pC2*n.conj(pC)/(n.conj(pC)*pC))
 
    #if pI == 0:
    #    pI = 0+1j

    if len(opts.output) > 0: outpath = opts.output+'/pspec_boot%04d.npz' % boot
    else: outpath = 'pspec_boot%04d.npz' % boot
    print '   Writing '+outpath
    n.savez(outpath, kpl=kpl, scalar=scalar, times=n.array(lsts),
        pk_vs_t=pC, err_vs_t=1./cnt, temp_noise_var=var, nocov_vs_t=pI,
        afreqs=afreqs,chans=chans,cmd=' '.join(sys.argv))


