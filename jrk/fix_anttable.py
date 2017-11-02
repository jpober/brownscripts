#!/usr/bin/env python
import aipy as a, numpy as np
import optparse, sys, os

o = optparse.OptionParser()
o.set_usage('fix_anttable.py [options] *.uv')
o.set_description(__doc__)
#a.scripting.add_standard_options(o,cal=True)
opts,args = o.parse_args(sys.argv[1:])

for uvfile in args:
    uvofile = uvfile+'T'
    print uvfile,'->',uvofile
    if os.path.exists(uvofile):
        print uvofile, 'exists, skipping.'
        continue
    ant1,ant2 = [],[]
    uvi = a.miriad.UV(uvfile)
    uvi.Npols = 2
    for (uvw,t,(i,j)),d in uvi.all():
        ant1.append(i)
        ant2.append(j)
    ants = list(set(ant1+ant2))
    print "found {n} antennas".format(n=len(ants))
    print "keeping the following antennas in the antenna table"
    print np.sort(ants)
    antpos = uvi['antpos']
    newantpos = []
    for ant in ants:
        newantpos.append(antpos[ant])
    newantpos = np.array(newantpos).T
    
    def mfunc(uv,p,d,f):
        return p,d,f
    #output the new file
    uvi.rewind()
    uvo = a.miriad.UV(uvofile,status='new')
    uvo.init_from_uv(uvi,override={'antpos':newantpos,'npol':2})
    uvo.pipe(uvi,mfunc=mfunc,raw=True,append2hist=' '.join(sys.argv)+'\n')

