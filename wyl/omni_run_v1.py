#! /usr/bin/env python

import omnical, aipy, numpy, capo
import pickle, optparse, os, sys

o = optparse.OptionParser()
o.set_usage('omni_run.py [options] *uvcRRE')
o.set_description(__doc__)
aipy.scripting.add_standard_options(o,cal=True,pol=True)
o.add_option('--calpar',dest='calpar',type='string',
            help='Path and name of POL.p ("xx.p") calpar file.')
o.add_option('--redinfo',dest='redinfo',type='string',default='',
            help='Path and name of .bin redundant info file.')
o.add_option('--omnipath',dest='omnipath',default='',type='string',
            help='Path to save .npz files. Include final / in path.')
o.add_option('--ba',dest='ba',default=None,
            help='Antennas to exclude, separated by commas.')
opts,args = o.parse_args(sys.argv[1:])

#Dictionary of calpar gains and files
pols = opts.pol.split(',')
files = {}
g0 = {} #firstcal gains
for pp,p in enumerate(pols):
    #dictionary of calpars per pol
    g0[p[0]] = {} #indexing by one pol letter instead of two
    if p in opts.calpar: #your supplied calpar matches a pol
        print 'Reading', opts.calpar
        cp = pickle.load(open(opts.calpar,'rb'))
        for i in xrange(cp[p].shape[1]): #loop over antennas
            g0[p[0]][i] = numpy.conj(cp[p][:,i] / numpy.abs(cp[p][:,i]))
    else: #looks for a calpar you haven't stated in the call
        new_cp = opts.calpar.split('.p')[0][:-2]+p+'.p' #XXX assumes calpar naming is *pol.p
        if os.path.exists(new_cp): #if it exists, use it
            print 'Reading', new_cp
            cp = pickle.load(open(new_cp,'rb'))
            for i in xrange(cp[p].shape[1]): #loop over antennas
                g0[p[0]][i] = numpy.conj(cp[p][:,i] / numpy.abs(cp[p][:,i]))
        elif len(list(set(p))) > 1: #if the crosspol first_cal is missing, don't worry
            #print '%s not found, but that is OK'%new_cp
            continue
        else: #if the linpol first_cal is missing, do worry
            raise IOError('Missing first_cal file %s'%new_cp)
        
for filename in args:
    files[filename] = {}
    for p in pols:
        fn = filename.split('.')
        fn[-2] = p
        files[filename][p] = '.'.join(fn)

#Create info
if opts.redinfo != '': #reading redinfo file
    print 'Reading',opts.redinfo
    info = omnical.info.RedundantInfoLegacy()
    print '   Getting reds from redundantinfo'
    info.fromfile(opts.redinfo)
else: #generate reds from calfile
    aa = aipy.cal.get_aa(opts.cal,numpy.array([.15]))
    print 'Getting reds from calfile'
    if opts.ba: #XXX assumes exclusion of the same antennas for every pol
        ex_ants = []
        for a in opts.ba.split(','):
            ex_ants.append(int(a))
        print '   Excluding antennas:',sorted(ex_ants)
    else: ex_ants = []
    info = capo.omni.aa_pos_to_info(aa, pols=list(set(''.join(pols))), ex_ants=ex_ants, crosspols=pols)
reds = info.get_reds()

### Omnical-ing! Loop Through Compressed Files ###
for f,filename in enumerate(args):
    file_group = files[filename] #dictionary with pol indexed files
    print 'Reading:'
    for key in file_group.keys(): print '   '+file_group[key]

    #pol = filename.split('.')[-2] #XXX assumes 1 pol per file
    timeinfo,d,f = capo.arp.get_dict_of_uv_data([file_group[key] for key in file_group.keys()], antstr='cross', polstr=opts.pol)
    t_jd = timeinfo['times']
    t_lst = timeinfo['lsts']
    freqs = numpy.arange(.1,.2,.1/len(d[d.keys()[0]][pols[0]][0]))
    SH = d.values()[0].values()[0].shape #shape of file data (ex: (19,203))
    data,wgts,xtalk = {}, {}, {}
    m2,g2,v2 = {}, {}, {}
    for p in g0.keys():
        for i in g0[p]: g0[p][i] = numpy.resize(g0[p][i],SH) #resize gains like data
    data = d #indexed by bl and then pol (backwards from everything else)
    for p in pols:
        wgts[p] = {} #weights dictionary by pol
        for bl in f: 
            i,j = bl
            wgts[p][(j,i)] = wgts[p][(i,j)] = numpy.logical_not(f[bl][p]).astype(numpy.int)
    print '   Logcal-ing' 
    m1,g1,v1 = capo.omni.redcal(data,info,gains=g0, removedegen=True) #SAK CHANGE REMOVEDEGEN
    print '   Lincal-ing'
    m2,g2,v2 = capo.omni.redcal(data, info, gains=g1, vis=v1, uselogcal=False, removedegen=True)
    xtalk = capo.omni.compute_xtalk(m2['res'], wgts) #xtalk is time-average of residual
    m2['history'] = 'OMNI_RUN: '+''.join(sys.argv) + '\n'
    m2['jds'] = t_jd
    m2['lsts'] = t_lst
    m2['freqs'] = freqs
    
    if len(pols)>1: #zen.jd.npz
        npzname = opts.omnipath+'.'.join(filename.split('/')[-1].split('.')[0:3])+'.npz'
    else: #zen.jd.pol.npz
        npzname = opts.omnipath+'.'.join(filename.split('/')[-1].split('.')[0:4])+'.npz'
    
    print '   Saving %s'%npzname
    capo.omni.to_npz(npzname, m2, g2, v2, xtalk)
    
