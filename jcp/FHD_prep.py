#! /usr/bin/env python

import aipy as a
import numpy as np
import sys,os,optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o,cal=True,src=True)
o.add_option('--onephs',action='store_true',help='Set the phase to a single pointing for the entire dataset.')
o.add_option('--uvfits',action='store_true',help='Perform the miriad task to convert to a uvfits file.')
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal,uv['sdf'],uv['sfreq'],uv['nchan'])

del(uv)

curtime = None
for filename in args:
    if filename.endswith('/'): filename = filename[:-1]
    print filename,'-->',filename+'M'
    if os.path.exists(filename+'M'):
        print '\tFile exists... skipping.'
        continue
    
    uvi = a.miriad.UV(filename)
    if opts.onephs:
        (uvw,t,bl),d = uvi.read()
        aa.set_jultime(t+5.*60./a.const.s_per_day)
        RA = str(aa.sidereal_time())
        dec= str(aa.lat)
        print "opts.onesrc is True: setting phase to %s_%s"%(RA,dec)
        opts.src = RA+'_'+dec
        epoch = (aa.epoch-36525.0)/365.2422 + 2000.
    
    if not opts.src is None:
        if not opts.src.startswith('zen'):
            srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
            src = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs).values()[0]
        else: src = 'z'
    else: src = None

    uvi.rewind()
    D = {}
    for (uvw,t,(i,j)),d in uvi.all():
        aa.set_active_pol(a.miriad.pol2str[uvi['pol']])
        if curtime != t:
            curtime = t
            aa.set_jultime(t)
            print t
            if not src is None and not type(src) == str: src.compute(aa)
        if i == j: continue
        
        try:
            _d = d.copy()
            d = aa.phs2src(d,src,i,j)
            d /= aa.passband(i,j)
            uvw = aa.get_baseline(i,j,src=src)
        except(a.phs.PointingError): d *= 0
    
        bl = a.miriad.ij2bl(i,j)
        p = (uvw,t,(i,j))
        if not t in D.keys(): D[t] = {}
        if not bl in D[t].keys(): D[t][bl] = {}
        if not uvi['pol'] in D[t][bl].keys(): D[t][bl][uvi['pol']] = p,d

    antpos = np.array([ant.pos for ant in aa])
    antpos.shape = antpos.shape[0]*antpos.shape[1]

    uvo = a.pol.UV(filename+'M',status='new')
    ra = src.get_params()['ra']
    #uvo.init_from_uv(uvi,override={'antpos':antpos,'obsra':ra,'ra':ra,'epoch':epoch*2000/36525.})
    uvo.init_from_uv(uvi,override={'antpos':antpos,'obsra':ra,'ra':ra,'epoch':epoch})
    uvo.add_var('restfreq','d')
    uvo['restfreq'] = uvi['sfreq']
    for t in D:
        for bl in D[t]:
            for pol in np.sort(D[t][bl].keys()):
                p,d = D[t][bl][pol]
                crd, t, (i,j) = p
                pi,pj = a.miriad.pol2str[pol]
                uvo.write_pol(a.miriad.pol2str[pol])
                uvo.write(p,d)
    uvo._wrhd('history',uvo['history'] + 'FHD_prep:'+' '.join(sys.argv)+'\n')
    del uvo,uvi

if opts.uvfits:
    for filename in args:
        os.system('fits in=%s op=uvout out=%s'%(filename+'M',filename+'M.uvfits'))
