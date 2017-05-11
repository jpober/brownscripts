#! /usr/bin/env python
import numpy as np
import omnical, aipy, capo
import optparse, os, sys, glob
from astropy.io import fits
import pickle, copy
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
o.add_option('--ex_bls',dest='ex_bls',default=None,
             help='baselines to exclude, separated by commas. eg:0_1,2_3')
o.add_option('--ex_ubls',dest='ex_ubls',default=None,
             help='exclude unique type of baselines,separated by comma, eg: 0_1,1_2')
o.add_option('--ftype', dest='ftype', default='', type='string',
            help='Type of the input file, .uvfits, or miriad, or fhd, to read fhd, simply type in the path/obsid')
o.add_option('--tave', dest='tave', default=False, action='store_true',
             help='choose to average data over time before calibration or not')
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
o.add_option('--cal_all', dest='cal_all', default='', type='string',
             help='options: model(use sky model to cal) ,copy(directly copy from FHD), or empty')
o.add_option('--projdegen', dest='projdegen', default=False, action='store_true',
             help='Toggle: project degeneracy to raw fhd solutions')
o.add_option('--fitdegen', dest='fitdegen', default=False, action='store_true',
             help='Toggle: project degeneracy to fitted fhd solutions')
o.add_option('--wgt_cal', dest='wgt_cal', default=False, action='store_true',
             help='Toggle: weight each gain by its amplitude before cal')
o.add_option('--divauto', dest='divauto', default=False, action='store_true',
             help='Toggle: use auto corr to weight visibilities before cal, need wgt_cal on')
o.add_option('--smooth', dest='smooth', default=False, action='store_true',
             help='Toggle: smooth data before cal by removing any signal beyond horizon, need divauto on')
o.add_option('--fhdpath', dest='fhdpath', default='/users/wl42/data/wl42/FHD_out/fhd_PhaseII_Longrun_EoR0/', type='string',
             help='path to fhd dir for projecting degen parameters, or fhd output visibilities if ftype is fhd.')
o.add_option('--metafits', dest='metafits', default='/users/wl42/data/wl42/Nov2016EoR0/', type='string',
             help='path to metafits files')
o.add_option('--ex_dipole', dest='ex_dipole', default=False, action='store_true',
             help='Toggle: exclude tiles which have dead dipoles')
o.add_option('--min_size', dest='min_size', default=40, type='int',
             help='minimun size of redundant groups to use to do diagnostic')
o.add_option('--sigma_tol', dest='sigma_tol', default=2.0, type='float',
             help='The tolerance of excluding bad vis data in diagnostic')
o.add_option('--snr', dest='snr', default=0.5, type='float',
             help='The tolerance of SNR for excluding bad red gp in diagnostic')
opts,args = o.parse_args(sys.argv[1:])

#Dictionary of calpar gains and files
pols = opts.pol.split(',')
files = {}
#files=[]
ex_bls = []
if opts.ex_bls:
    ex_bl_list = opts.ex_bls.split(',')
    for bl in ex_bl_list:
        _bl = bl.split('_')
        ex_bls.append((int(_bl[0]),int(_bl[1])))
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

if opts.projdegen or opts.fitdegen or opts.cal_all == 'model' or opts.cal_all == 'copy':
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
        filelist = glob.glob(opts.fhdpath+'/vis_data/'+filename+'*')+glob.glob(opts.fhdpath+'/metadata/'+filename+'*')
        if len(pols) == 1:
            p = pols[0]
            if p == 'xx':
                try: filelist.remove(opts.fhdpath+'/vis_data/'+filename + '_vis_YY.sav')
                except: pass
            elif p == 'yy':
                try: filelist.remove(opts.fhdpath+'/vis_data/'+filename + '_vis_XX.sav')
                except: pass
            else:
                raise IOError('do not support cross pol')
            files[filename][p] = filelist
        else:
            files[filename][filename] = filelist
    else:
        raise IOError('invalid filetype, it should be miriad, uvfits, or fhd')

exec('from %s import *'% opts.cal) # Including antpos, realpos, EastHex, SouthHex

if opts.cal_all == 'model':
    print "   Loading model"
    model_files = glob.glob(opts.fhdpath+'vis_data/'+args[0]+'*') + glob.glob(opts.fhdpath+'metadata/'+args[0]+'*')
    model_dict = capo.wyl.uv_read_omni([model_files],filetype='fhd', antstr='cross', p_list=pols, use_model=True)
#################################################################################################

ex_ubls = []
for ubls in opts.ex_ubls.split(','):
    i,j = ubls.split('_')
    ex_ubls.append((int(i),int(j)))

def diagnostic(infodict):
    min_size_bl_gp = int(opts.min_size)
    exclude_bls = []
    g0 = infodict['g0']
    p = infodict['pol']
    d = infodict['data']
    f = infodict['flag']
    ginfo = infodict['ginfo']
    freqs = infodict['freqs']
    ex_ants = infodict['ex_ants']
    info = capo.omni.pos_to_info(antpos, pols=list(set(''.join([p]))), ex_ants=ex_ants, ex_ubls=ex_ubls, crosspols=[p])
    reds = info.get_reds()
    data = {}
    for bl in d.keys():
        i,j = bl
        if not (i in info.subsetant and j in info.subsetant): continue
        m = np.ma.masked_array(d[bl][p],mask=f[bl][p])
        m = np.mean(m,axis=0)
        data[bl] = {p: np.complex64(m.data.reshape(1,-1))}
    #if txt file or first cal is not provided, g0 is initiated here, with all of them to be 1.0
    if opts.calpar == None:
        if not g0.has_key(p[0]): g0[p[0]] = {}
        for iant in range(0, ginfo[0]):
            g0[p[0]][iant] = np.ones((1,ginfo[2]))
    elif opts.calpar.endswith('.sav'):
        for key in g0[p[0]].keys():
            g0_temp = g0[p[0]][key]
            g0[p[0]][key] = np.resize(g0_temp,(1,ginfo[2]))
    elif opts.calpar.endswith('.npz'):
        for key in g0[p[0]].keys():
            g0_temp = g0[p[0]][key]
            g0[p[0]][key] = np.resize(g0_temp,(1,ginfo[2]))
            if opts.initauto: g0[p[0]][key] *= auto[key]
    m1,g1,v1 = capo.omni.redcal(data,info,gains=g0, removedegen=False)
    m2,g2,v2 = capo.omni.redcal(data, info, gains=g1, vis=v1, uselogcal=False, removedegen=False)
    for r in reds:
        if len(r) < min_size_bl_gp: continue
        stack_data = []
        stack_bl = []
        for bl in r:
            stack_bl.append(bl)
            try: stack_data.append(data[bl][p][0]/(g2[p[0]][bl[0]][0]*g2[p[0]][bl[1]][0].conj()))
            except(KeyError): stack_data.append(data[bl[::-1]][p][0]/(g2[p[0]][bl[1]][0]*g2[p[0]][bl[0]][0].conj()))
        stack_data = np.array(stack_data)
        stack_bl = np.array(stack_bl)
        vis_std = np.nanstd(stack_data,axis=0)
        vis_ave = np.nanmean(stack_data,axis=0)
        nonzero = np.where(vis_std>0)
        if np.mean(np.abs(vis_ave)[nonzero]/vis_std[nonzero]) < opts.snr: exclude_bls = exclude_bls + r
        else:
            n_sigmas = np.nanmean(np.abs(stack_data-vis_ave)/vis_std,axis=1)
            ind = np.where(n_sigmas > float(opts.sigma_tol))
            for ii in ind[0]: exclude_bls.append(tuple(stack_bl[ii]))
    return [exclude_bls, g2]



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
    mask_arr = infodict['mask']
    amp_wgt = infodict['amp_wgt'][p[0]]
    for bl in ex_bls:
        if not bl in infodict['ex_bls']: infodict['ex_bls'].append(bl)
    print 'Getting reds from calfile'
    print 'generating info:'
    filter_length = None
    info = capo.omni.pos_to_info(antpos, pols=list(set(''.join([p]))), ex_ants=ex_ants, ex_ubls=ex_ubls, ex_bls=infodict['ex_bls'], crosspols=[p])

    reds = info.get_reds()
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

    # organize data for redundant cal
    for bl in d.keys():
        i,j = bl
        if not (i in info.subsetant and j in info.subsetant): continue
        if opts.tave:
            m = np.ma.masked_array(d[bl][p],mask=f[bl][p])
            m = np.mean(m,axis=0)
            data[bl] = {p: np.complex64(m.data.reshape(1,-1))}
        else: data[bl] = {p: copy.copy(d[bl][p])}
    if opts.wgt_cal:
        # if weight antennas by auto corr:
        if opts.divauto:
            for bl in data.keys():
                i,j = bl
                data[bl][p] /= (auto[i]*auto[j])
                if opts.smooth:
                    tfq = np.fft.fftfreq(freqs.size,(freqs[1]-freqs[0]))
                    fftdata = np.fft.fft(data[bl][p],axis=1)
                    ri,rj = realpos[i],realpos[j]
                    rij = np.array([ri['top_x']-rj['top_x'],ri['top_y']-rj['top_y'],ri['top_z']-rj['top_z']])
                    inds = np.where(np.abs(tfq)>(np.linalg.norm(rij)/3e8+50e-9))
                    fftdata[:,inds]=0
                    data[bl][p] = np.complex64(np.fft.ifft(fftdata,axis=1))
        else:
            for bl in data.keys():
                i,j = bl
                data[bl][p] /= (np.abs(amp_wgt[i][0])*np.abs(amp_wgt[j][0]))
     #indexed by bl and then pol (backwards from everything else)

    wgts[p] = {} #weights dictionary by pol
    for bl in f:
        i,j = bl
        wgts[p][(j,i)] = wgts[p][(i,j)] = np.logical_not(f[bl][p]).astype(np.int)
    print '   Logcal-ing'
    m1,g1,v1 = capo.omni.redcal(data,info,gains=g0, removedegen=opts.removedegen) #SAK CHANGE REMOVEDEGEN
    print '   Lincal-ing'
    m2,g2,v2 = capo.omni.redcal(data, info, gains=g1, vis=v1, uselogcal=False, removedegen=opts.removedegen)
    if opts.wgt_cal:
        if opts.divauto:
            for a in g2[p[0]].keys():
                g2[p[0]][a] *= auto[a]
        else:
            for a in g2[p[0]].keys():
                g2[p[0]][a] *= np.abs(amp_wgt[a][0])
    xtalk = capo.omni.compute_xtalk(m2['res'], wgts) #xtalk is time-average of residual

    ############# To project out degeneracy parameters ####################
    if opts.projdegen or opts.fitdegen:
        fuse = []
        for ff in range(384):
            if not ff%16 in [0,15]: fuse.append(ff)
        print '   Projecting degeneracy'
        ref = min(g2[p[0]].keys()) # pick a reference tile to reduce the effect of phase wrapping, it has to be a tile in east hex
        ref_exp = np.exp(1j*np.angle(g2[p[0]][ref][:,fuse]/gfhd[p[0]][ref][fuse]))
        for a in g2[p[0]].keys(): g2[p[0]][a][:,fuse] /= ref_exp
        print '   projecting amplitude'
        amppar = capo.wyl.ampproj(g2,gfhd)
        print '   projecting phase'
        phspar = capo.wyl.phsproj(g2,gfhd,realpos,EastHex,SouthHex,ref)
        degen_proj = {}
        for a in g2[p[0]].keys():
            dx = realpos[a]['top_x']-realpos[ref]['top_x']
            dy = realpos[a]['top_y']-realpos[ref]['top_y']
            proj = amppar[p[0]]*np.exp(1j*(dx*phspar[p[0]]['phix']+dy*phspar[p[0]]['phiy']))
#            if a < 93 and ref > 92: proj *= phspar[p[0]]['offset_east']
            if a > 92: proj *= phspar[p[0]]['offset_south']
            degen_proj[a] = proj
            g2[p[0]][a] *= proj
        for bl in v2[p].keys():
            i,j = bl
            degenij = (degen_proj[j].conj()*degen_proj[i])
            fuse = np.where(degenij!=0)
            fnot = np.where(degenij==0)
            v2[p][bl][fuse] /= degenij[fuse]
            v2[p][bl][fnot] *= 0
    #compute chi-square
    if not opts.tave:
        print '   compute chi-square'
        chisq = 0
        for r in reds:
            for bl in r:
                if v2[p].has_key(bl): yij = v2[p][bl]
            for bl in r:
                try: md = np.ma.masked_array(d[bl][p],mask=f[bl][p])
                except(KeyError): md = np.ma.masked_array(d[bl[::-1]][p].conj(),mask=f[bl[::-1]][p],fill_value=0.0)
                i,j = bl
                chisq += (np.abs(md.data-g2[p[0]][i]*g2[p[0]][j].conj()*yij))**2/(np.var(md,axis=0).data+1e-7)
        DOF = (info.nBaseline - info.nAntenna - info.ublcount.size)
        m2['chisq'] = chisq / float(DOF)
    for a in g2[p[0]].keys():
        if opts.tave:
            g2[p[0]][a] = np.resize(g2[p[0]][a],(ginfo[2]))
            stack_mask = np.sum(np.logical_not(mask_arr),axis=0).astype(bool)
            g2[p[0]][a] *= stack_mask
        else:
            chi = m2['chisq']
            chi_mask = np.zeros(chi.shape,dtype=bool)
            ind = np.where(chi>1.2)
            chi_mask[ind] = True
            or_mask = np.logical_or(chi_mask,mask_arr)
            g_temp = np.ma.masked_array(g2[p[0]][a],or_mask,fill_value=0.0)
            g_temp = np.mean(g_temp,axis=0)
            g2[p[0]][a] = g_temp.data
            if opts.instru == 'mwa':
                for ii in range(384):
                    if ii%16 == 8: g2[p[0]][a][ii] = (g2[p[0]][a][ii+1]+g2[p[0]][a][ii-1])/2
    ###########################################################################################
    m2['history'] = 'OMNI_RUN: '+''.join(sys.argv) + '\n'
    m2['jds'] = t_jd
    m2['lsts'] = t_lst
    m2['freqs'] = freqs
    m2['ex_bls'] = infodict['ex_bls']
    if opts.cal_all == 'model':
        print '   start absolute cal'
        ref = min(g2[p[0]].keys())
        g2 = capo.wyl.absoulte_cal(d,g2,model_dict[p],realpos,ref,ex_ants=ex_ants)
    elif opts.cal_all == 'copy':
        print '   copying non-hex cal solution from FHD run'
        for a in gfhd[p[0]].keys():
            if a in g2[p[0]].keys() or a in ex_ants: continue
            if np.isnan(np.mean(gfhd[p[0]][a])): g2[p[0]][a] = np.zeros(gfhd[p[0]][a].shape,dtype=np.complex)
            else: g2[p[0]][a] = gfhd[p[0]][a]
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
            dict0 = capo.wyl.uv_read_omni([filegroup[p]], filetype = 'miriad', antstr='cross', p_list=[p])
            infodict[p] = dict0[p]
            infodict[p]['filename'] = filegroup[p]
            infodict['name_dict'] = dict0['name_dict']
    else:
        infodict = capo.wyl.uv_read_omni([filegroup[key] for key in filegroup.keys()], filetype=opts.ftype, antstr='cross', p_list=pols)
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
        if opts.projdegen or opts.fitdegen or opts.cal_all == 'model' or opts.cal_all == 'copy':
            for a in gfhd[p[0]].keys():
                if np.isnan(np.mean(gfhd[p[0]][a])):
                    if not a in infodict[p]['ex_ants']: infodict[p]['ex_ants'].append(a)
        if opts.ex_dipole:
            metafits_path = opts.metafits + filename + '.metafits'
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
    print "  Start Diagnostic:"
    par1 = Pool(2)
    list_exclude_bls_g2 = par1.map(diagnostic, info_dict)
    par1.close()
    for ii in range(len(info_dict)):
        print '   excluded baselines:', list_exclude_bls_g2[ii][0]
        info_dict[ii]['ex_bls'] = list_exclude_bls_g2[ii][0]
        info_dict[ii]['amp_wgt'] = list_exclude_bls_g2[ii][1]
    print "  Start Calibration:"
    par2 = Pool(2)
    npzlist = par2.map(calibration, info_dict)
    par2.close()
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



