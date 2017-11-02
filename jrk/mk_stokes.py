#! /usr/bin/env python
import aipy as a, numpy as n
import capo as C
import optparse, sys, os

POL_WGTS = {
    'I': {'xx': 1. , 'yy': 1. },
    'Q': {'xx': 1. , 'yy':-1. },
    'U': {'xy': 1. , 'yx': 1. },
    'V': {'xy':-1.j, 'yx': 1.j},
}

LIN2STOKES = {
    'xx':'I',
    'yy':'Q',
    'xy':'U',
    'yx':'V',
}

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True)
o.add_option('--stokes', help='Comma-delimited list of Stokes parameters to form (e.g. I,Q,U,V).')
opts,args = o.parse_args(sys.argv[1:])

stokes = opts.stokes.split(',')
pol = {}
for s in stokes: pol.update(POL_WGTS[s])
pol = ','.join(pol.keys())
uv = a.miriad.UV(args[0])
ants = a.scripting.parse_ants(opts.ant, uv['nants'])
#aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)

for filename in args:
    outfile = filename + 'P'
    print filename, '->', outfile
    if os.path.exists(outfile):
        print '    File exists, skipping.'
        continue
    dsum,dwgt = {}, {}
    for s in stokes: dsum[s], dwgt[s] = {}, {}
    curtime = None
    print '    Converting %s to %s' % (pol, ','.join(stokes))
    uvi = a.miriad.UV(filename)
    #a.scripting.uv_selector(uvi, ants, pol)
    a.scripting.uv_selector(uvi, ants)
    #uvi.select('decimate', opts.decimate, opts.decphs)
    #uvi.Npols = 2
    for (crd,t,(i,j)),d,f in uvi.all(raw=True):
        bl = a.miriad.ij2bl(i,j)
        if t != curtime:
            curtime = t
            #aa.set_jultime(t)
            for s in stokes:
                dsum[s][t], dwgt[s][t] = {}, {}
        for s in stokes:
            try: wgt = POL_WGTS[s][a.miriad.pol2str[uvi['pol']]]
            except(KeyError): continue
            dsum[s][t][bl] = dsum[s][t].get(bl, 0) + n.where(f, 0, wgt*d)
            dwgt[s][t][bl] = dwgt[s][t].get(bl, 0) + n.abs(wgt)*n.logical_not(f).astype(n.int)
    uvi.rewind()
    
    print '    Writing output file'
    uvo = a.miriad.UV(outfile, status='new')
    uvo.init_from_uv(uvi, override={'pol':a.miriad.str2pol['I'],'npol':1})
    def mfunc(uv, p, d, f):
        uvw,t,(i,j) = p
        bl = a.miriad.ij2bl(i,j)
        s = LIN2STOKES[a.miriad.pol2str[uv['pol']]]
        #print s, a.miriad.pol2str[uv['pol']]
        try: _dsum,_dwgt = dsum[s][t].pop(bl), dwgt[s][t].pop(bl)
        except(KeyError): return p, None, None
        uvo['pol'] = a.miriad.str2pol[s]
        #print uvo['pol']
        uvo['npol']=1
        wgt = _dwgt.clip(1,n.Inf)
        f = n.where(_dwgt == 0, 1, 0)
        d = _dsum / wgt
        return (uvw,t,(i,j)), d, f

    uvo.pipe(uvi, mfunc=mfunc, raw=True, append2hist='COMBINE_POL:' + ' '.join(sys.argv) + '\n')
    del(uvi); del(uvo)




