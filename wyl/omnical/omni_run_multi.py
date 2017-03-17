#! /usr/bin/env python
import numpy as np
import omnical, aipy, capo
import optparse, os, sys, glob
from astropy.io import fits
import pickle
from multiprocessing import Pool
from scipy.io.idl import readsav
#from IPython import embed

o = optparse.OptionParser()
o.set_usage('omni_run_multi.py [options] *uvcRRE/obsid') #only takes 1 obsid
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
o.add_option('--flength',dest='flength',default=None,
             help='a threshold for baseline lengths to use, in meters')
o.add_option('--ftype', dest='ftype', default='', type='string',
            help='Type of the input file, .uvfits, or miriad, or fhd, to read fhd, simply type in the path/obsid')
o.add_option('--tave', dest='tave', default=False, action='store_true',
             help='choose to average data over time before calibration or not')
o.add_option('--gave', dest='gave', default=False, action='store_true',
             help='choose to average solution over time after calibration or not')
o.add_option('--iftxt', dest='iftxt', default=False, action='store_true',
            help='Toggle: write the npz info to a ucla format txt file or not')
o.add_option('--iffits', dest='iffits', default=False, action='store_true',
            help='Toggle: write the npz info to a ucla format fits file or not')
o.add_option('--removedegen',dest='removedegen',default=False,action='store_true',
             help='Toggle: turn removedegen on')
o.add_option('--initauto',dest='initauto',default=False,action='store_true',
             help='Toggle: use auto_corr as initial guess for gains')
o.add_option('--instru', dest='instru', default='mwa', type='string',
             help='instrument type. Default=mwa')
o.add_option('--projdegen', dest='projdegen', default=False, action='store_true',
             help='Toggle: project degeneracy to raw fhd solutions')
o.add_option('--fitdegen', dest='fitdegen', default=False, action='store_true',
             help='Toggle: project degeneracy to fitted fhd solutions')
o.add_option('--divauto', dest='divauto', default=False, action='store_true',
             help='Toggle: use auto corr to weight visibilities before cal')
o.add_option('--fhdpath', dest='fhdpath', default='/users/wl42/data/wl42/FHD_out/fhd_PhaseII_EoR0/', type='string',
             help='path to fhd solutions for projecting degen parameters. Default=/path/to/calibration/')
o.add_option('--metafits', dest='metafits', default='/users/wl42/data/wl42/EoR0_PhaseII/', type='string',
             help='path to metafits files')
o.add_option('--ex_dipole', dest='ex_dipole', default=False, action='store_true',
             help='Toggle: exclude tiles which have dead dipoles')
opts,args = o.parse_args(sys.argv[1:])

#Dictionary of calpar gains and files
pols = opts.pol.split(',')
files = {}
#files=[]
g0 = {} #firstcal gains
g_scale = {'x': 1.0, 'y': 1.0}
if opts.calpar != None: #create g0 if txt file is provided
    fname = opts.calpar
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
                g0[pp][ant] = np.array(g0[pp][ant])
                g0[pp][ant] = g0[pp][ant].reshape(len(Ntimes),len(Nfreqs))
    elif fname.endswith('.npz'):
        for pp,p in enumerate(pols):
            g0[p[0]] = {}   #obs(or jds).pol.fc.npz
            fpname = fname.split('.')
            fpname[-3] = p
            fpname = '.'.join(fpname)
            print '   Reading: ', fpname
            cp = np.load(fpname)
            for i in cp.keys():
                if i[0].isdigit():
                    g0[p[0]][int(i[:-1])] = cp[i] / np.abs(cp[i])
    elif fname.endswith('.fits'):
        g0 = capo.omni.fc_gains_from_fits(opts.calpar)
        for key1 in g0:
            for key2 in g0[key1]:
                g0[key1][key2] /= np.abs(g0[key1][key2])
    elif fname.endswith('.sav'):
        cal = readsav(opts.calpar,python_dict=True)
        fqfl = np.zeros((128,384),dtype=bool)
        for ff in range(384):
            if ff%16==0 or ff%16==15: fqfl[:,ff]=True
        g = cal['cal']['GAIN'][0]
        g0['x'] = {}
        g0['y'] = {}
        gx = np.ma.masked_array(g[0],fqfl,fill_value=1.0)
        gy = np.ma.masked_array(g[1],fqfl,fill_value=1.0)
        gnan = np.where(np.isnan(np.mean(gx,axis=1)))[0]
        g_scale['x'] = np.nanmean(np.abs(gx[56:-1]))
        g_scale['y'] = np.nanmean(np.abs(gy[56:-1]))
        for nn in gnan:
            gx[nn] = g_scale['x']
            gy[nn] = g_scale['y']
        for ii in range(0,cal['cal']['N_TILE'][0]):
            g0['x'][ii] = gx[ii].filled()
            g0['y'][ii] = gy[ii].filled()
    else:
        raise IOError('invalid calpar file')

if opts.projdegen or opts.fitdegen:
    fhd_cal = readsav(opts.fhdpath+'calibration/'+args[0]+'_cal.sav',python_dict=True)
    gfhd = {'x':{},'y':{}}
    if opts.fitdegen:
        for a in range(fhd_cal['cal']['N_TILE'][0]):
            gfhd['x'][a] = fhd_cal['cal']['GAIN'][0][0][a]
            gfhd['y'][a] = fhd_cal['cal']['GAIN'][0][1][a]
    else:
        for a in range(fhd_cal['cal']['N_TILE'][0]):
            gfhd['x'][a] = fhd_cal['cal']['GAIN'][0][0][a] + fhd_cal['cal']['GAIN_RESIDUAL'][0][0][a]
            gfhd['y'][a] = fhd_cal['cal']['GAIN'][0][1][a] + fhd_cal['cal']['GAIN_RESIDUAL'][0][1][a]

for filename in args:
    files[filename] = {}
    if opts.ftype == 'miriad':
        for p in pols:
            fn = filename.split('.')
            fn[-2] = p
            files[filename][p] = '.'.join(fn)
    elif opts.ftype == 'uvfits':
        files[filename][filename] = filename + '.uvfits'
    elif opts.ftype == 'fhd':
        obs = filename + '*'
        filelist = glob.glob(obs)
        if len(pols) == 1:
            p = pols[0]
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
            files[filename][filename] = filelist
    else:
        raise IOError('invalid filetype, it should be miriad, uvfits, or fhd')

exec('from %s import *'% opts.cal) # Including antpos, realpos, EastHex, SouthHex

if opts.instru == 'mwa':
    print "   Loading model"
    model_files = glob.glob(opts.fhdpath+'vis_data/'+args[0]+'*') + glob.glob(opts.fhdpath+'metadata/'+args[0]+'*')
    model_dict = capo.wyl.uv_read_omni([model_files],filetype='fhd', antstr='cross', p_list=pols, use_model=True)
#################################################################################################
def calibration(infodict):#dict=[filename, g0, timeinfo, d, f, ginfo, freqs, pol, auto_corr]
    filename = infodict['filename']
    g0 = infodict['g0']
    p = infodict['pol']
    d = infodict['data']
    f = infodict['flag']
    ginfo = infodict['ginfo']
    freqs = infodict['freqs']
    timeinfo = infodict['timeinfo']
    ex_ants = infodict['ex_ants']
    auto = infodict['auto_corr']
    print 'Getting reds from calfile'
    print 'generating info:'
    filter_length = None
    if not opts.flength == None: filter_length = float(opts.flength)
    info = capo.omni.pos_to_info(antpos, pols=list(set(''.join([p]))), filter_length=filter_length, ex_ants=ex_ants, crosspols=[p])

    ### Omnical-ing! Loop Through Compressed Files ###

    print '   Calibrating ' + p + ': ' + filename
    
    #if txt file or first cal is not provided, g0 is initiated here, with all of them to be 1.0
    if opts.calpar == None:
        if not g0.has_key(p[0]): g0[p[0]] = {}
        for iant in range(0, ginfo[0]):
            g0[p[0]][iant] = np.ones((1,ginfo[2]))
            if opts.initauto: g0[p[0]][iant] *= auto[iant]
            if opts.tave: g0[p[0]][iant] = np.mean(g0[p[0]][iant],axis=0)
    elif opts.calpar.endswith('.sav'):
        for key in g0[p[0]].keys():
            g0_temp = g0[p[0]][key]
            if opts.tave: g0[p[0]][key] = np.resize(g0_temp,(1,ginfo[2]))
            else: g0[p[0]][key] = np.resize(g0_temp,(ginfo[1],ginfo[2]))
    elif opts.calpar.endswith('.npz'):
        for key in g0[p[0]].keys():
            g0_temp = g0[p[0]][key]
            if opts.tave: g0[p[0]][key] = np.resize(g0_temp,(1,ginfo[2]))
            else: g0[p[0]][key] = np.resize(g0_temp,(ginfo[1],ginfo[2]))
            if opts.initauto: g0[p[0]][key] *= auto[key]

    t_jd = timeinfo['times']
    t_lst = timeinfo['lsts']

    data,wgts,xtalk = {}, {}, {}
    m2,g2,v2 = {}, {}, {}
    if opts.divauto:
        for bl in d.keys():
            i,j = bl
            data[bl] = {}
            data[bl][p] = d[bl][p]/(auto[i]*auto[j])
    else: data = d #indexed by bl and then pol (backwards from everything else)

    wgts[p] = {} #weights dictionary by pol
    for bl in f:
        i,j = bl
        wgts[p][(j,i)] = wgts[p][(i,j)] = np.logical_not(f[bl][p]).astype(np.int)
    print '   Logcal-ing' 
    m1,g1,v1 = capo.omni.redcal(data,info,gains=g0, removedegen=opts.removedegen) #SAK CHANGE REMOVEDEGEN
    print '   Lincal-ing'
    m2,g2,v2 = capo.omni.redcal(data, info, gains=g1, vis=v1, uselogcal=False, removedegen=opts.removedegen)
    if opts.divauto:
        for a in g2[p[0]].keys():
            g2[p[0]][a] *= auto[a]
    if opts.tave:
        for a in g2[p[0]].keys():
            g2[p[0]][a] = np.resize(g2[p[0]][a],(ginfo[1],ginfo[2]))
        for bl in v2[p].keys():
            v2[p][bl] = np.resize(v2[p][bl],(ginfo[1],ginfo[2]))
    if opts.gave:
        for a in g2[p[0]].keys():
            gmean = np.mean(g2[p[0]][a],axis=0)
            g2[p[0]][a] = np.resize(gmean,(ginfo[1],ginfo[2]))
    xtalk = capo.omni.compute_xtalk(m2['res'], wgts) #xtalk is time-average of residual
    ############# correct the center of each coarse band and coarse band edge if instrument is mwa ######################
#    if opts.instru == 'mwa':
#        for a in g2[p[0]].keys():
#            for ff in range(0,384):
#                if ff%16==8:
#                    g2[p[0]][a][:,ff] = (g2[p[0]][a][:,ff-1]+g2[p[0]][a][:,ff+1])/2
#                if ff%16 in [0,15]:
#                    g2[p[0]][a][:,ff] = 0
    ############# To project out degeneracy parameters ####################
    if opts.projdegen or opts.fitdegen:
        print '   Projecting degeneracy'
        for a in g2[p[0]].keys():
            if g2[p[0]][a].ndim == 2 : g2[p[0]][a] = np.mean(g2[p[0]][a][1:53],axis=0)
        ref = g2[p[0]].keys()[0] # pick a reference tile to reduce the effect of phase wrapping
        ref_exp = np.exp(1j*np.angle(g2[p[0]][ref]/gfhd[p[0]][ref]))
        for a in g2[p[0]].keys(): g2[p[0]][a] /= ref_exp
        print '   projecting amplitude'
        amppar = capo.wyl.ampproj(g2,gfhd)
        print '   projecting phase'
        phspar = capo.wyl.phsproj(g2,gfhd,realpos,EastHex,SouthHex,ref)
        for a in g2[p[0]].keys():
            dx = realpos[a]['top_x']-realpos[ref]['top_x']
            dy = realpos[a]['top_y']-realpos[ref]['top_y']
            proj = amppar[p[0]]*np.exp(1j*(dx*phspar[p[0]]['phix']+dy*phspar[p[0]]['phiy']))
            if a < 92: proj *= phspar[p[0]]['offset_east']
            else: proj *= phspar[p[0]]['offset_south']
            g2[p[0]][a] *= proj
        print '   linear projecting'
        lp = capo.wyl.linproj(g2,gfhd,realpos)
        for a in g2[p[0]].keys():
            dx = realpos[a]['top_x']/100
            dy = realpos[a]['top_y']/100
            proj = np.exp(lp[p[0]]['eta']+1j*(dx*lp[p[0]]['phix']+dy*lp[p[0]]['phiy']+lp[p[0]]['offset']))
            g2[p[0]][a] *= proj
            g2[p[0]][a] = np.resize(g2[p[0]][a],(ginfo[1],ginfo[2]))
            for ff in range(384):
                if ff%16 in [0,15]:
                    g2[p[0]][a][:,ff] = 0    #clean nans
    ###########################################################################################
    m2['history'] = 'OMNI_RUN: '+''.join(sys.argv) + '\n'
    m2['jds'] = t_jd
    m2['lsts'] = t_lst
    m2['freqs'] = freqs
    if opts.instru == 'mwa':
        print '   start non-hex tiles calibration using fhd model'
        g3 = capo.wyl.non_hex_cal(d,g2,model_dict[p],realpos,ex_ants=ex_ants)
        for a in g3[p[0]].keys():
            if not g2[p[0]].has_key(a): g2[p[0]][a] = g3[p[0]][a]
    if opts.ftype == 'miriad':
        npzname = opts.omnipath+'.'.join(filename.split('/')[-1].split('.')[0:4])+'.npz'
    else:
        npzname = opts.omnipath+filename.split('/')[-1]+'.'+ p +'.npz'

    print '   Saving %s'%npzname
    capo.omni.to_npz(npzname, m2, g2, v2, xtalk)
    return npzname
#######################################################################################################

for f,filename in enumerate(args):

    npzlist = []
    infodict = {}
    filegroup = files[filename]
    info_dict = []
    print "  Reading data: " + filename
    if opts.ftype == 'miriad':
        for p in pols:
            dict0 = capo.wyl.uv_read_omni([filegroup[p]], filetype = 'miriad', antstr='cross', p_list=[p], tave=opts.tave)
            infodict[p] = dict0[p]
            infodict[p]['filename'] = filegroup[p]
            infodict['name_dict'] = dict0['name_dict']
    else:
        infodict = capo.wyl.uv_read_omni([filegroup[key] for key in filegroup.keys()], filetype=opts.ftype, antstr='cross', p_list=pols, tave=opts.tave)
        for p in pols:
            infodict[p]['filename'] = filename
    print "  Finish reading."
    for p in pols:
        if opts.calpar == None:
            infodict[p]['g0'] = {}
        else:
            infodict[p]['g0'] = {}
            infodict[p]['g0'][p[0]] = g0[p[0]]
        if opts.ba:
            for a in opts.ba.split(','):
                if not int(a) in infodict[p]['ex_ants']:
                    infodict[p]['ex_ants'].append(int(a))
        if opts.ex_dipole:
            metafits_path = opts.metafits + args[0] + '.metafits'
            if os.path.exists(metafits_path):
                print '    Finding dead dipoles in metafits'
                hdu = fits.open(metafits_path)
                inds = np.where(hdu[1].data['Delays']==32)[0]
                dead_dipole = np.unique(hdu[1].data['Antenna'][inds])
                for dip in dead_dipole:
                    if not dip in infodict[p]['ex_ants']:
                        infodict[p]['ex_ants'].append(dip)
            else: print '    Warning: Metafits not found. Cannot get the information of dead dipoles'
        ex_ants = sorted(infodict[p]['ex_ants'])
        print '   Excluding antennas:', ex_ants

        info_dict.append(infodict[p])
    print "  Start Parallelism:"
    par = Pool(2)
    npzlist = par.map(calibration, info_dict)
    par.close()
    name_dict = infodict['name_dict']

    if opts.iftxt: #if True, write npz gains to txt files
        scrpath = os.path.abspath(sys.argv[0])
        pathlist = os.path.split(scrpath)[0].split('/')
        repopath = '/'.join(pathlist[0:-1])+'/'
        print '   Writing to txt:'
        capo.wyl.writetxt(npzlist, repopath, ex_ants=ex_ants, name_dict=name_dict)
        print '   Finish'

    if opts.iffits: #if True, write npz gains to fits files
        scrpath = os.path.abspath(sys.argv[0])
        pathlist = os.path.split(scrpath)[0].split('/')
        repopath = '/'.join(pathlist[0:-1])+'/'
        print '   Writing to fits:'
        capo.wyl.writefits(npzlist, repopath, ex_ants=ex_ants, name_dict=name_dict)
        print '   Finish'



