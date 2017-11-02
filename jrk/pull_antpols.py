#! /usr/bin/env python

import aipy as a, sys, optparse, os

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True)
opts,args = o.parse_args(sys.argv[1:])

for filename in args:
    outfile = os.path.basename(filename)+'A'
    print filename, '->', outfile
    uvi = a.miriad.UV(filename)
    if os.path.exists(outfile):
        print '    File exists... skipping.'
        continue
   
    if not opts.pol is None: a.scripting.uv_selector(uvi, ants=opts.ant,pol_str=opts.pol)
    else:  a.scripting.uv_selector(uvi, ants=opts.ant)
    uvo = a.miriad.UV(outfile,status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi,append2hist='PULL ANTPOLS:'+' '.join(sys.argv)+'\n')
