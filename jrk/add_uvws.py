#!/usr/bin/env python
import aipy as a, numpy as np
import optparse, sys, os

o = optparse.OptionParser()
o.set_usage('add_uvws.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o,cal=True)
opts,args = o.parse_args(sys.argv[1:])


for uvfile in args:
    uvofile = uvfile+'U'
    print uvfile,'->',uvofile
    if os.path.exists(uvofile):
        print uvofile, 'exists, skipping.'
        continue

    #aa = a.phs.ArrayLocation(('-30:43:17.5','21:25:41.9'))
    aa = a.cal.get_aa(opts.cal,np.array([0.15]))
    nints = 0
    curtime = None
    def mfunc(uv, preamble, data, flags):
        global curtime
        global nints
        uvw, t, (i,j) = preamble
        uvw = aa.get_baseline(i,j,'z')
        preamble = (uvw, t, (i,j))
        return preamble, data, flags

    uvi = a.miriad.UV(uvfile)
    uvo = a.miriad.UV(uvofile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=mfunc, raw=True, append2hist=' '.join(sys.argv)+'\n')
