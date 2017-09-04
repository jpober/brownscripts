import numpy as np
import sys, glob, mp2cal, optparse
from astropy.io import fits
delays = {
'0,5,10,15,1,6,11,16,2,7,12,17,3,8,13,18':-5,
'0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15':-4,
'0,3,6,9,0,3,6,9,0,3,6,9,0,3,6,9':-3,
'0,2,4,6,0,2,4,6,0,2,4,6,0,2,4,6':-2,
'0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3':-1,
'0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0':0,
'3,2,1,0,3,2,1,0,3,2,1,0,3,2,1,0':1,
'6,4,2,0,6,4,2,0,6,4,2,0,6,4,2,0':2,
'9,6,3,0,9,6,3,0,9,6,3,0,9,6,3,0':3,
'12,8,4,0,13,9,5,1,14,10,6,2,15,11,7,3':4,
'15,10,5,0,16,11,6,1,17,12,7,2,18,13,8,3':5,
}

o = optparse.OptionParser()
o.set_usage('average.py [options]')
o.set_description(__doc__)
o.add_option('--fhd_path', dest='fhd_path', default='/users/wl42/data/wl42/Nov2016EoR0/fhd_sol/', help='path to fhd cal solutions')
o.add_option('--dat_path', dest='dat_path', default='/users/wl42/data/wl42/Nov2016EoR0/', help='path to fhd cal solutions')
o.add_option('--ap',dest='ap', default=False, action='store_true', help='average in amplitude and phase, otherwise in real and imag parts, Default=False')
opts,args = o.parse_args(sys.argv[1:])
#p = sys.argv[1]
pols = ['xx','yy']
meta, vismdl, xtalk = {},{},{}
fuse = []
for ii in range(384):
    if not ii%16 in [0,15]: fuse.append(ii)
fn=glob.glob(opts.fhd_path+'*.xx.fhd.npy')
fn.sort()
obsids = []
for f in fn:
    obsids.append(f.split('/')[-1].split('.')[0])
g = {}
#    nfiles = {}
for obs in obsids:
    metafits = opts.dat_path+obs+'.metafits'
    hdu = fits.open(metafits)
    day = int(obs)/86400
    suffix = str(day)+'_'+str(delays[hdu[0].header['DELAYS']])
    if not g.has_key(suffix): g[suffix]={'x':[],'y':[]}
    solx = np.load(opts.fhd_path+obs+'.xx.fhd.npy')
    soly = np.load(opts.fhd_path+obs+'.yy.fhd.npy')
    g[suffix]['x'].append(solx)
    g[suffix]['y'].append(soly)

for suffix in g.keys():
    g[suffix]['x'] = np.array(g[suffix]['x'])
    g[suffix]['y'] = np.array(g[suffix]['y'])
    if opts.ap:
        ampx = np.abs(g[suffix]['x'])
        ampy = np.abs(g[suffix]['y'])
        phsx = np.zeros(g[suffix]['x'].shape)
        phsy = np.zeros(g[suffix]['y'].shape)
        phsx[:,:,fuse] = np.unwrap(np.angle(g[suffix]['x'][:,:,fuse]))
        phsy[:,:,fuse] = np.unwrap(np.angle(g[suffix]['y'][:,:,fuse]))
        ampx = np.nanmean(ampx,axis=0)
        ampy = np.nanmean(ampy,axis=0)
        phsx = np.nanmean(phsx,axis=0)
        phsy = np.nanmean(phsy,axis=0)
        g[suffix]['x'] = ampx*np.exp(1j*phsx)
        g[suffix]['y'] = ampy*np.exp(1j*phsy)
    else:
        g[suffix]['x'] = np.nanmean(g[suffix]['x'],axis=0)
        g[suffix]['y'] = np.nanmean(g[suffix]['y'],axis=0)
    gx = {'x': {}}
    gy = {'y': {}}
    for a in range(128):
        gx['x'][a] = g[suffix]['x'][a]
        gy['y'][a] = g[suffix]['y'][a]
    outfnx = 'fhdave_'+str(suffix)+'.xx.npz'
    outfny = 'fhdave_'+str(suffix)+'.yy.npz'
    mp2cal.wyl.save_gains_omni(outfnx, meta, gx, vismdl, xtalk)
    mp2cal.wyl.save_gains_omni(outfny, meta, gy, vismdl, xtalk)

