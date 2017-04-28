#!/bin/env python

## Uses pyuvdata's select function to pull specific antennas from input files, appending "A" to the end.

from pyuvdata import UVData
import sys, optparse,os
import aipy as a

o = optparse.OptionParser()
o.add_option("-a", "--ants", help="Antennas to select (AIPY options)")

opts,args = o.parse_args(sys.argv[1:])

uvd = UVData()
for filename in args:
    print filename, '->', filename+'A'
    if os.path.exists(filename+'A'):
        print '    File exists... skipping.'
        continue
    uvd.read_miriad(filename)
    
    baselines = ("_" in opts.ants)   ### If baselines are listed specifically, select by baselines. Otherwise, limit to the antennas listed in ants.

    if baselines:
       ants=a.scripting.parse_ants(opts.ants,uvd.Nants_telescope)
       bl_list=[]
       for cnt,(bl,include,pol) in enumerate(ants):
            i,j = uvd.baseline_to_antnums(bl)
            bl_list.append((int(i),int(j)))
       uvd.select(ant_pairs_nums=bl_list)
    else:
       ants = list(set(map(int,opts.ants.split(","))))
       uvd.select(antenna_nums=ants)
    uvd.write_miriad(filename+'A')

#    a.scripting.uv_selector(uvi, ants=opts.ant)
#    if not opts.pol is None: a.scripting.uv_selector(uvi, pol_str=opts.pol)
#    uvo = a.miriad.UV(filename+'A',status='new')
#    uvo.init_from_uv(uvi)
#    uvo.pipe(uvi,append2hist='PULL ANTPOLS:'+' '.join(sys.argv)+'\n')

