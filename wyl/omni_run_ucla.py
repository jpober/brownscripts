#! /usr/bin/env python

import omnical, aipy, numpy, capo
import pickle, optparse, os, sys
import uvdata.uv as uvd, math

#####################################################################################################
def uv_read(filenames, filetype=None, polstr=None,antstr=None,recast_as_array=True):
    info = {'lsts':[], 'times':[]}
    dat, flg = {},{}
    ginfo = [0,0,0]
    #    uvdata=uvd.UVData()
    if type(filenames) == 'str': filenames = [filenames]
    for filename in filenames:
        uvdata = uvd.UVData()
        if filetype == 'miriad':
            uvdata.read_miriad(filename)
        elif filetype == 'uvfits':
            uvdata.read_uvfits(filename)
        elif filetype == 'fhd':
            uvdata.read_fhd(filenames)
        else:
            raise IOError('invalid filetype, it should be miriad, uvfits, or fhd')
        #uvdata.read(filename, filetype)
        tt = uvdata.time_array.value
        Nt = uvdata.Ntimes.value
        blt = len(tt)
        nbl = uvdata.Nbls.value
        nfreq = uvdata.Nfreqs.value
        
        for ii in range(0,Nt):
            info['times'].append(tt[ii*nbl])
            info['lsts'].append(tt[ii*nbl])   #not sure how to calculate lsts
        pol = uvdata.polarization_array.value
        npol = len(pol)
        data = uvdata.data_array.value
        flag = uvdata.flag_array.value
        ant1 = uvdata.ant_1_array.value
        ant2 = uvdata.ant_2_array.value
        freqarr = uvdata.freq_array.value[0]
        
        nant = int((1+math.sqrt(1+8*nbl))/2)

#ginfo=[nant, Nt, nfreq]
        ginfo[0] = nant
        ginfo[1] = Nt
        ginfo[2] = nfreq
        
        for ii in range(0,blt):
            bl = (ant1[ii],ant2[ii])
            if antstr == 'cross':
                if ant1[ii] == ant2[ii]: continue
            if not dat.has_key(bl): dat[bl],flg[bl] = {},{}
            for jj in range(0,npol):
                pp = aipy.miriad.pol2str[pol[jj]]
                if polstr != None:
                    if pp != polstr: continue
                if not dat[bl].has_key(pp):
                    dat[bl][pp],flg[bl][pp] = [],[]
                data00,flag00 = [],[]
                for nn in range(0,len(data[ii][jj])):
                    data00.append(data[ii][jj][nn][0])
                    flag00.append(flag[ii][jj][nn][0])
                dat[bl][pp].append(data00)
                flg[bl][pp].append(flag00)
        if filetype == 'fhd': break
    if recast_as_array:
        for ii in dat.keys():
            for jj in dat[ii].keys():
                dat[ii][jj] = numpy.array(dat[ii][jj])
                flg[ii][jj] = numpy.array(flg[ii][jj])
        info['lsts'] = numpy.array(info['lsts'])
        info['times'] = numpy.array(info['times'])
    return info, dat, flg, ginfo, freqarr
#####################################################################################################

o = optparse.OptionParser()
o.set_usage('omni_run.py [options] *uvcRRE')
o.set_description(__doc__)
aipy.scripting.add_standard_options(o,cal=True,pol=True)
o.add_option('--calpar',dest='calpar',type='string',default=None,
            help='Path and name of POL.p ("xx.p") calpar file.')
o.add_option('--redinfo',dest='redinfo',type='string',default='',
            help='Path and name of .bin redundant info file.')
o.add_option('--omnipath',dest='omnipath',default='',type='string',
            help='Path to save .npz files. Include final / in path.')
o.add_option('--ba',dest='ba',default=None,
            help='Antennas to exclude, separated by commas.')
o.add_option('--ftype', dest='ftype', default='', type='string',
            help='Type of the input file, .uvfits, or miriad, or fhd, to read fhd, use ',' to separate different save files')
o.add_option('--iftxt', dest='iftxt', default=False, action='store_true',
            help='A switch to write the npz info to a ucla txt file or not')
opts,args = o.parse_args(sys.argv[1:])

#Dictionary of calpar gains and files
pols = opts.pol.split(',')
files = {}
fhdfiles = []
#files=[]
g0 = {} #firstcal gains
if not opts.calpar == None: #create g0 if txt file is provided
    fname = opts.calpar
    if fname.endswith('.txt'):
        f = open(fname,'r')
        tkn = []
        g = {}
        for line in f:
            temp = line.split(',')[:7]
            if temp[0].startswith('#'): continue
            temp2 = []
            for ii, s in enumerate(temp):
                if ii == 0: continue
                elif s.strip() == 'EE': s = 'xx' #need to check the convension
                elif s.strip() == 'NN': s = 'yy'
                temp2.append(s)
            if temp2[2].strip() == 'EN' or temp2[2].strip() == 'NE': continue
            temp3=[temp2[2], int(temp2[0]), float(temp2[3]), float(temp2[1]), float(temp2[4]), float(temp2[5])]  #temp3=[pol,ant,jds,freq,real,imag]
            tkn.append(temp3)
        for p,pp in enumerate(tkn):
            if not g.has_key(pp[0][0]):
                g[pp[0][0]] = {}
                g0[pp[0][0]] = {}
            if not g[pp[0][0]].has_key(pp[1]):
                g[pp[0][0]][pp[1]] = {}
                g0[pp[0][0]][pp[1]] = []
            if not g[pp[0][0]][pp[1]].has_key(pp[2]):
                g[pp[0][0]][pp[1]][pp[2]] = {}
            if not g[pp[0][0]][pp[1]][pp[2]].has_key(pp[3]):
                gg = complex(pp[4],pp[5])
                g[pp[0][0]][pp[1]][pp[2]][pp[3]] = gg.conjugate()/abs(gg)#g{pol:{ant:{jds:{freq:gain}}}}
        for i1, pol in pols:
            for i2, ant in enumerate(g[pol]):
                for i3, jds in enumerate(g[pol][ant]):
                    ff = []
                    for i4, freq in enumerate(g[pol][ant][jds]):
                        ff.append(g[pol][ant][jds][freq])
                    g0[pol][ant].append(ff)
                g0[pol][ant] = numpy.array(g0[pol][ant])      #g0={pol:{ant:{array(jds,freq)}}}
    else:
        raise IOError('invalid txtfile')
#if not provided, will initiate g0 with units in the reading file part

for filename in args:
    if opts.ftype == 'uvfits' or opts.ftype == 'miriad':
        files[filename] = {}
        for p in pols:
            fn = filename.split('.')
            fn[-2] = p
            files[filename][p] = '.'.join(fn)
    elif opts.ftype == 'fhd':
        fhdfiles.append(filename)
    else:
        raise IOError('invalid filetype, it should be miriad, uvfits, or fhd')

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

    if opts.ftype == 'uvfits' or opts.ftype == 'miriad':
        file_group = files[filename] #dictionary with pol indexed files
        print 'Reading:'
        for key in file_group.keys(): print '   '+file_group[key]

    #pol = filename.split('.')[-2] #XXX assumes 1 pol per file
        timeinfo,d,f,ginfo,freqs = uv_read([file_group[key] for key in file_group.keys()], filetype=opts.ftype, polstr=opts.pol, antstr='cross')

    elif opts.ftype == 'fhd':
        print 'Reading:'
        print filename
        timeinfo,d,f,ginfo,freqs = uv_read(filename.split(','), filetype=opts.ftype, polstr=opts.pol, antstr='cross')
    
    #if txt file is not provided, g0 is initiated here, with all of them to be 1.0
    if opts.calpar == None:
        for p in pols:
            if not g0.has_key(p[0]): g0[p[0]] = {}
            for iant in range(0, ginfo[0]):
                g0[p[0]][iant] = numpy.ones((ginfo[1],ginfo[2]))

    t_jd = timeinfo['times']
    t_lst = timeinfo['lsts']

    data,wgts,xtalk = {}, {}, {}
    m2,g2,v2 = {}, {}, {}
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

    if opts.ftype == 'uvfits' or opts.ftype == 'miriad':
        if len(pols)>1: #zen.jd.npz
            npzname = opts.omnipath+'.'.join(filename.split('/')[-1].split('.')[0:3])+'.npz'
        else: #zen.jd.pol.npz
            npzname = opts.omnipath+'.'.join(filename.split('/')[-1].split('.')[0:4])+'.npz'

    elif opts.ftype == 'fhd':
        if len(pols)>1: #obsid.npz
            npzname = opts.omnipath+filename.split('/')[-1].split('_')[0]+'.npz'
        else: #obsid.pol.npz
            npzname = opts.omnipath+filename.split('/')[-1].split('_')[0]+pols[0]+'.npz'

    print '   Saving %s'%npzname
    capo.omni.to_npz(npzname, m2, g2, v2, xtalk)


    if opts.iftxt: #if True, write npz gains to txt files
        print 'writing to txt:'
        capo.omni.writetxt([npzname])
        print 'Saving txt file'


    
