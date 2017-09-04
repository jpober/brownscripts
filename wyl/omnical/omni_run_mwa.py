#! /usr/bin/env python
import numpy as np
import heracal, aipy, mp2cal
import optparse, os, sys, glob
from astropy.io import fits
import pickle, copy
from multiprocessing import Pool
from scipy.io.idl import readsav
import pyuvdata.uvdata as uvd
#from IPython import embed

o = optparse.OptionParser()
o.set_usage('omni_run_multi.py [options] obsid') #only takes 1 obsid
o.set_description(__doc__)
aipy.scripting.add_standard_options(o,cal=True,pol=True)
o.add_option('--ubls', dest='ubls', default='', help='Unique baselines to use, separated by commas (ex: 1_4,64_49).')
o.add_option('--ex_ubls', dest='ex_ubls', default='', help='Unique baselines to exclude, separated by commas (ex: 1_4,64_49).')
o.add_option('--bls', dest='bls', default='', help='Baselines to use, separated by commas (ex: 1_4,64_49).')
o.add_option('--ex_bls', dest='ex_bls', default='', help='Baselines to exclude, separated by commas (ex: 1_4,64_49).')
o.add_option('--ants', dest='ants', default='', help='Antennas to use, separated by commas (ex: 1,4,64,49).')
o.add_option('--ex_ants', dest='ex_ants', default='', help='Antennas to exclude, separated by commas (ex: 1,4,64,49).')
o.add_option('--ftype', dest='ftype', default='', type='string', help='Type of the input file, uvfits or fhd save files')
o.add_option('--omnipath',dest='omnipath',default='',type='string', help='Path to load firstcal files and to save solution files. Include final / in path.')
o.add_option('--fhdpath', dest='fhdpath', default='/users/wl42/data/wl42/FHD_out/fhd_PhaseII_Longrun_EoR0/', type='string', help='path to fhd dir for projecting degen parameters, or fhd output visibilities if ftype is fhd.')
o.add_option('--metafits', dest='metafits', default='/users/wl42/data/wl42/Nov2016EoR0/', type='string', help='path to metafits files')
o.add_option('--tave', dest='tave', default=False, action='store_true', help='choose to average data over time before calibration or not')
o.add_option('--projdegen', dest='projdegen', default=False, action='store_true', help='Toggle: Project degen to FHD solutions')
o.add_option('--ex_dipole', dest='ex_dipole', default=False, action='store_true', help='Toggle: exclude tiles which have dead dipoles')
o.add_option('--wgt_cal', dest='wgt_cal', default=False, action='store_true', help='Toggle: weight each gain by auto corr before cal')
opts,args = o.parse_args(sys.argv[1:])

#*****************************************************************************
def bl_parse(bl_args):
    bl_list = []
    for bl in bl_args.split(','):
        i,j = bl.split('_')
        bl_list.append((int(i),int(j)))
    return bl_list

def ant_parse(ant_args):
    ant_list = []
    for a in ant_args.split(','): ant_list.append(int(a))
    return ant_list

#************************************************************************************************
exec('from %s import *'% opts.cal) # Including antpos, realpos, EastHex, SouthHex
pols = opts.pol.split(',')
ubls, ex_ubls, bls, ex_bls, ants, ex_ants = None, [], None, [], None, []
if not opts.ubls == '':
    ubls = bl_parse(opts.ubls)
    print '     ubls: ', ubls
if not opts.ex_ubls == '':
    ex_ubls = bl_parse(opts.ex_ubls)
    print '     ex_ubls: ', ex_ubls
if not opts.bls == '':
    bls = bl_parse(opts.bls)
    print '     bls: ', bls
if not opts.ex_bls == '':
    ex_bls = bl_parse(opts.ex_bls)
    print '     ex_bls: ', ex_bls
if not opts.ants == '':
    ants = ant_parse(opts.ants)
    print '   ants: ', ants
if not opts.ex_ants == '': ex_ants = ant_parse(opts.ex_ants)

#********************************** load and wrap data ******************************************
if not len(args) == 1: raise IOError('Do not support multiple files.')
obsid = args[0]
uv = uvd.UVData()
if opts.ftype == 'uvfits':
    uv.read_uvfits(obsid+'.uvfits',run_check=False,run_check_acceptability=False)
elif opts.ftype == 'fhd':
    uv.read_fhd(glob.glob(opts.fhdpath+'/vis_data/'+obsid+'*')+glob.glob(opts.fhdpath+'/metadata/'+obsid+'*'),use_model=False,run_check=False,run_check_acceptability=False)
else: IOError('invalid filetype, it should be uvfits or fhd')
data_wrap = mp2cal.wyl.uv_wrap_omni(uv,pols=pols)
t_jd = uv.time_array[::uv.Nbls]
t_lst = uv.lst_array[::uv.Nbls]
freqs = uv.freq_array[0]
SH = (uv.Ntimes, uv.Nfreqs)

#********************************** load fhd ***************************************************
gfhd = mp2cal.wyl.load_gains_fhd(opts.fhdpath+'calibration/'+obsid+'_cal.sav')

#*********************************** ex_ants *****************************************************
ex_ants_find = mp2cal.wyl.find_ex_ant(uv)
del uv
for a in ex_ants_find:
    if not a in ex_ants: ex_ants.append(a)
for a in gfhd['x'].keys():
    if np.isnan(np.mean(gfhd['x'][a])) or np.isnan(np.mean(gfhd['y'][a])):
        if not a in ex_ants: ex_ants.append(a)
if opts.ex_dipole:
    metafits_path = opts.metafits + filename + '.metafits'
    if os.path.exists(metafits_path):
        print '    Finding dead dipoles in metafits'
        hdu = fits.open(metafits_path)
        inds = np.where(hdu[1].data['Delays']==32)[0]
        dead_dipole = np.unique(hdu[1].data['Antenna'][inds])
        for dip in dead_dipole:
            if not dip in ex_ants: ex_ants.append(dip)
    else: print '    Warning: Metafits not found. Cannot get the information of dead dipoles'
print '     ex_ants: ', ex_ants

#################################################################################################
data_list = []
for pp in pols: data_list.append(data_wrap[pp])

def omnirun(data_wrap):
    pp = data_wrap['pol']
    p = pp[0]
    data = data_wrap['data']
    flag = data_wrap['flag']
    auto = data_wrap['auto']
    mask_arr = data_wrap['mask']
    flag_bls = []
    for bl in flag.keys():
        wgt_data = np.logical_not(flag[bl][pp])
        wgt_data = np.sum(wgt_data,axis=0)
        ind = np.where(wgt_data==0)
        if ind[0].size > 48: flag_bls.append(bl)
    print 'flagged baselines: ', flag_bls
    omnisol = opts.omnipath + obsid + '.' + pp + '.omni.npz'
    info = mp2cal.wyl.pos_to_info(antpos,pols=[p],fcal=False,ubls=ubls,ex_ubls=ex_ubls,bls=bls,ex_bls=ex_bls+flag_bls,ants=ants,ex_ants=ex_ants)
    reds = info.get_reds()
    redbls = [bl for red in reds for bl in red]

    #*********************** organize data *************************************
    dat,wgts,xtalk = {}, {}, {}
    for bl in data.keys():
        i,j = bl
        if not (i in info.subsetant and j in info.subsetant): continue
        if bl in ex_bls+flag_bls: continue
        if opts.tave:
            m = np.ma.masked_array(data[bl][pp],mask=flag[bl][pp])
            m = np.mean(m,axis=0)
            dat[bl] = {pp: np.complex64(m.data.reshape(1,-1))}
        else: dat[bl] = {pp: np.copy(data[bl][pp])}
        if opts.wgt_cal: dat[bl][pp] /= (auto[i]*auto[j])

    #*********************** generate g0 ***************************************
    if opts.ftype == 'fhd':
        print '     setting g0 as units'
        g0 = {p:{}}
        for a in info.subsetant: g0[p][a] = np.ones((freqs.size),dtype=np.complex64)
    else:
        fcfile = opts.omnipath + obsid + '.' + pp + '.fc.npz'
        if os.path.exists(fcfile):
            print '     loading firstcal file: ', fcfile
            g0 = mp2cal.wyl.load_gains_fc(fcfile)
        else:
            print '     firstcal not found, start rough cal'
            info_rough = mp2cal.wyl.pos_to_info(antpos,pols=[p],ubls=[(57,61),(57,62)],ex_ants=ex_ants)
            g0 = mp2cal.wyl.rough_cal(dat,info_rough,pol=pp)

    #*********************** Calibrate ******************************************
    wgts[pp] = {} #weights dictionary by pol
    for bl in flag:
        i,j = bl
        wgts[pp][(j,i)] = wgts[pp][(i,j)] = np.logical_not(flag[bl][pp]).astype(np.int)
    print '   Run omnical'
#    m2,g2,v2 = mp2cal.wyl.run_omnical(dat,info,gains0=g0, maxiter=500, conv=1e-9)
    m2,g2,v2 = heracal.omni.run_omnical(dat,info,gains0=g0, maxiter=500, conv=1e-12)
    if opts.wgt_cal:
        for a in g2[p].keys(): g2[p][a] *= auto[a]
    xtalk = heracal.omni.compute_xtalk(m2['res'], wgts) #xtalk is time-average of residual

    #************************ Average cal solutions ************************************
    if not opts.tave:
        print '   compute chi-square'
        chisq = 0
        for r in reds:
            for bl in r:
                if v2[pp].has_key(bl): yij = v2[pp][bl]
            for bl in r:
                try: md = np.ma.masked_array(data[bl][pp],mask=flag[bl][pp])
                except(KeyError): md = np.ma.masked_array(data[bl[::-1]][pp].conj(),mask=flag[bl[::-1]][pp],fill_value=0.0)
                i,j = bl
                chisq += (np.abs(md.data-g2[p][i]*g2[p][j].conj()*yij))**2/(np.var(md,axis=0).data+1e-7)
        DOF = (info.nBaseline - info.nAntenna - info.ublcount.size)
        m2['chisq'] = chisq / float(DOF)
        chi = m2['chisq']
        chi_mask = np.zeros(chi.shape,dtype=bool)
        ind = np.where(chi>1.2)
        chi_mask[ind] = True
        or_mask = np.logical_or(chi_mask,mask_arr)
    for a in g2[p].keys():
#        ind = np.where(np.isnan(g2[p][a]))
#        g2[p][a][ind] = 0
        if opts.tave:
            g2[p][a] = np.resize(g2[p][a],(SH[1]))
            stack_mask = np.sum(np.logical_not(mask_arr),axis=0).astype(bool)
            g2[p[0]][a] *= stack_mask
        else:
            g_temp = np.ma.masked_array(g2[p][a],or_mask,fill_value=0.0)
            g_temp = np.mean(g_temp,axis=0)
            g2[p[0]][a] = g_temp.data
#            for ii in range(384):
#                if ii%16 == 8: g2[p[0]][a][ii] = (g2[p[0]][a][ii+1]+g2[p[0]][a][ii-1])/2

    #*********************** project degeneracy *********************************
    if opts.projdegen:
        print '   Projecting degeneracy'
        if opts.ftype == 'fhd':
            g2 = mp2cal.wyl.degen_project_FO(g2,antpos,v2)
        elif opts.ftype == 'uvfits':
            g2 = mp2cal.wyl.degen_project_OF(g2,gfhd,antpos,EastHex,SouthHex,v2)
    else: g2 = mp2cal.wyl.remove_degen_hex(g2, antpos)

    #************************* metadata parameters ***************************************
    m2['history'] = 'OMNI_RUN: '+' '.join(sys.argv) + '\n'
    m2['jds'] = t_jd
    m2['lsts'] = t_lst
    m2['freqs'] = freqs

    #************************** Saving cal ************************************************
    print '     saving %s' % omnisol
    mp2cal.wyl.save_gains_omni(omnisol,m2,g2,v2,xtalk)

par = Pool(2)
npzlist = par.map(omnirun,data_list)
par.close()
