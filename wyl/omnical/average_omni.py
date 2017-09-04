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
exec('from PhaseII_cal import *')
o = optparse.OptionParser()
o.set_usage('average.py [options]')
o.set_description(__doc__)
o.add_option('--scale', dest='scale', default=False, action='store_true', help='scale the gains before average,Default=False')
o.add_option('--ap',dest='ap', default=False, action='store_true', help='average in amplitude and phase, otherwise in real and imag parts, Default=False')
opts,args = o.parse_args(sys.argv[1:])
#p = sys.argv[1]
pols = ['xx','yy']
meta, vismdl, xtalk = {},{},{}
fuse = []
for ii in range(384):
    if not ii%16 in [0,15]: fuse.append(ii)
for p in pols:
    fn=glob.glob('./*'+p+'.omni.npz')
    g = {}
    fid = {}
#    nfiles = {}
    for f in fn:
        gains = mp2cal.wyl.quick_load_gains(f)
        if opts.scale: gains = mp2cal.wyl.scale_gains(gains)
        obs = f.split('/')[-1].split('.')[0]
        metafits = '../'+obs+'.metafits'
        hdu = fits.open(metafits)
        day = int(obs)/86400
        suffix = str(day)+'_'+str(delays[hdu[0].header['DELAYS']])
        if not g.has_key(suffix):
            g[suffix]={p[0]:{}}
            fid[suffix]={p[0]:gains}
        else:
#            gains = mp2cal.wyl.degen_project_OF(gains,fid[suffix][p[0]],antpos,EastHex,SouthHex)
            gains = mp2cal.wyl.degen_project_simple(gains,fid[suffix][p[0]],antpos)
#        if not nfiles.has_key(suffix): nfiles[suffix]=0
#        nfiles[suffix]+=1
        for a in gains[p[0]].keys():
            if np.isnan(np.mean(gains[p[0]][a])): continue
            if not g[suffix][p[0]].has_key(a): g[suffix][p[0]][a] = []
            g[suffix][p[0]][a].append(gains[p[0]][a])
        del gains
#    print nfiles
    for suf in g.keys():
        for a in g[suf][p[0]].keys():
            g[suf][p[0]][a] = np.array(g[suf][p[0]][a])
            mask = np.zeros(g[suf][p[0]][a].shape,dtype=bool)
            ind = np.where(g[suf][p[0]][a]==0)
            mask[ind] = True
            if opts.ap:
                amp = np.ma.masked_array(np.abs(g[suf][p[0]][a]),mask,fill_value=0.0)
                gphs = np.angle(g[suf][p[0]][a])
                gphs[:,fuse] = np.unwrap(gphs[:,fuse])
                phs = np.ma.masked_array(gphs,mask,fill_value=0.0)
                g[suf][p[0]][a] = (np.mean(amp,axis=0)).data*np.exp(1j*(np.mean(phs,axis=0)).data)
            else:
                mg = np.ma.masked_array(g[suf][p[0]][a],mask,fill_value=0.0)
                g[suf][p[0]][a] = (np.mean(mg,axis=0)).data
        outfn = 'omniave_'+str(suf)+'.'+p+'.npz'
        mp2cal.wyl.save_gains_omni(outfn, meta, g[suf], vismdl, xtalk)

