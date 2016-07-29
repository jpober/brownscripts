#! /usr/bin/env python

import omnical, aipy, numpy, capo
import optparse, os, sys, glob
from astropy.io import fits
import pickle

o = optparse.OptionParser()
o.set_usage('omni_run.py [options] *uvcRRE')
o.set_description(__doc__)
aipy.scripting.add_standard_options(o,cal=True,pol=True)
o.add_option('--calpar',dest='calpar',type='string',default=None,
            help='Path and name of calpar file (txt or npz).')
o.add_option('--redinfo',dest='redinfo',type='string',default='',
            help='Path and name of .bin redundant info file.')
o.add_option('--omnipath',dest='omnipath',default='',type='string',
            help='Path to save .npz files. Include final / in path.')
o.add_option('--ba',dest='ba',default=None,
            help='Antennas to exclude, separated by commas.')
o.add_option('--ftype', dest='ftype', default='', type='string',
            help='Type of the input file, .uvfits, or miriad, or fhd, to read fhd, simply type in the path/obsid')
o.add_option('--iftxt', dest='iftxt', default=False, action='store_true',
            help='A switch to write the npz info to a ucla format txt file or not')
o.add_option('--iffits', dest='iffits', default=False, action='store_true',
            help='A switch to write the npz info to a ucla format fits file or not')
opts,args = o.parse_args(sys.argv[1:])

#Dictionary of calpar gains and files
pols = opts.pol.split(',')
files = {}
#files=[]
g0 = {} #firstcal gains
if opts.calpar != None: #create g0 if txt file is provided
    fname = opts.calpar
    print '   Reading: ', fname
    if fname.endswith('.txt'):
        f = open(fname,'r')
        Ntimes = []
        Nfreqs = []
        for line in f:
            temp = line.split(',')[:7]
            if temp[0].startswith('#'): continue
            temp2 = []
            for ii, s in enumerate(temp):
                if ii == 0: continue
                elif s.strip() == 'EE': s = 'xx' #need to check the convension
                elif s.strip() == 'NN': s = 'yy'
                elif s.strip() == 'EN': s = 'xy'
                elif s.strip() == 'NE': s = 'yx'
                temp2.append(s)
            if not temp2[2].strip() in pols: continue
            temp3 = [temp2[2], int(temp2[0]), float(temp2[3]), float(temp2[1]), float(temp2[4]), float(temp2[5])]  #temp3=[pol,ant,jds,freq,real,imag]
            if not temp3[2] in Ntimes: Ntimes.append(temp3[2])
            if not temp3[3] in Nfreqs: Nfreqs.append(temp3[3])
            if not g0.has_key(temp3[0][0]):
                g0[temp3[0][0]] = {}
            if not g0[temp3[0][0]].has_key(temp3[1]):
                g0[temp3[0][0]][temp3[1]] = []
            gg = complex(temp3[4],temp3[5])
            g0[temp3[0][0]][temp3[1]].append(gg.conjugate()/abs(gg))
        for pp in g0.keys():
            for ant in g0[pp].keys():
                g0[pp][ant] = numpy.array(g0[pp][ant])
                g0[pp][ant] = g0[pp][ant].reshape(len(Ntimes),len(Nfreqs))
    elif fname.endswith('.npz'):
        for pp,p in enumerate(pols):
            g0[p[0]] = {}
            if p in fname:
                print '   Reading: ', fname
                cp = numpy.load(fname)
                for i in cp.keys():
                    if i.isdigit():
                        g0[p[0]][int(i)] = cp[i] / numpy.abs(cp[i])
            else:
                new_cp = fname.split('.npz')[0][:-2]+p+'.npz'
                print '   Reading: ', new_cp
                cp = numpy.load(new_cp)
                for i in cp.keys():
                    if i.isdigit():
                        g0[p[0]][int(i)] = cp[i] / numpy.abs(cp[i])
    elif fname.endswith('.fits'):
        poldict = {'EE': 'xx', 'NN': 'yy', 'EN': 'xy', 'NE': 'yx'}
        hdu = fits.open(fname)
        Ntimes = hdu[0].header['NTIMES']
        Nfreqs = hdu[0].header['NFREQS']
        Npols = hdu[0].header['NPOLS']
        Nants = hdu[0].header['NANTS']
        ant_index = hdu[1].data['ANT INDEX'][0:Nants]
        pol_list = hdu[1].data['POL'][0:Nfreqs*Nants*Npols].reshape(Npols,Nants*Nfreqs)[:,0]
        data_list = hdu[1].data['GAIN'].reshape((Ntimes,Npols,Nfreqs,Nants)).swapaxes(0,1).swapaxes(2,3).swapaxes(1,2) #Npols,Nants,Ntimes,Nfreqs
        for ii in range(0,Npols):
            polarization = poldict[pol_list[ii]]
            if not polarization in pols: continue
            g0[polarization[0]] = {}
            for jj in range(0,Nants):
                g0[polarization[0]][ant_index[jj]]=numpy.conj(data_list[ii][jj]/numpy.abs(data_list[ii][jj]))
    else:
        raise IOError('invalid calpar file')

#if not provided, will initiate g0 with units in the reading file part

for filename in args:
    files[filename] = {}
    for p in pols:
        if opts.ftype == 'uvfits' or opts.ftype == 'miriad':
            fn = filename.split('.')
            fn[-2] = p
            files[filename][p] = '.'.join(fn)
        elif opts.ftype == 'fhd':
            obs = filename + '*'
            filelist = glob.glob(obs)
            if p == 'xx':
                try: filelist.remove(filename + '_vis_YY.sav')
                except: pass
            elif p == 'yy':
                try: filelist.remove(filename + '_vis_XX.sav')
                except: pass
            else:
                raise IOError('do not support cross pol')
            files[filename][p] = filelist
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


    file_group = files[filename] #dictionary with pol indexed files
    print '   Reading:'
    for key in file_group.keys(): print '   ', file_group[key]

    #pol = filename.split('.')[-2] #XXX assumes 1 pol per file
    timeinfo,d,f,ginfo,freqs = capo.wyl.uv_read([file_group[key] for key in file_group.keys()], filetype=opts.ftype, polstr=opts.pol, antstr='cross')
    
    #if txt file or first cal is not provided, g0 is initiated here, with all of them to be 1.0
    if opts.calpar == None:
        for p in pols:
            if not g0.has_key(p[0]): g0[p[0]] = {}
            for iant in range(0, ginfo[0]):
                g0[p[0]][iant] = numpy.ones((ginfo[1],ginfo[2]))
    elif opts.calpar.endswith('.npz'):
        SH = (ginfo[1],ginfo[2])
        for p in g0.keys():
            for i in g0[p]: g0[p][i] = numpy.resize(g0[p][i],SH)

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
            npzname = opts.omnipath+filename.split('/')[-1]+'.npz'
        else: #obsid.pol.npz
            npzname = opts.omnipath+filename.split('/')[-1]+'.'+pols[0]+'.npz'

    print '   Saving %s'%npzname
    capo.omni.to_npz(npzname, m2, g2, v2, xtalk)


    if opts.iftxt: #if True, write npz gains to txt files
        scrpath = os.path.abspath(sys.argv[0])
        pathlist = os.path.split(scrpath)[0].split('/')
        repopath = '/'.join(pathlist[0:-1])+'/'
        print '   Writing to txt:'
        capo.wyl.writetxt([npzname], repopath, ex_ants)
        print '   Finish'

    if opts.iffits: #if True, write npz gains to fits files
        scrpath = os.path.abspath(sys.argv[0])
        pathlist = os.path.split(scrpath)[0].split('/')
        repopath = '/'.join(pathlist[0:-1])+'/'
        print '   Writing to fits:'
        capo.wyl.writefits([npzname], repopath, ex_ants)
        print '   Finish'



