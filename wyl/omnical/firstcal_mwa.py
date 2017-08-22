# edited version of firstcal in heracal
#! /usr/bin/env python
import heracal, mp2cal
import pylab as p, aipy as a
import sys,optparse,glob,os
import numpy as np
from multiprocessing import Pool
import pyuvdata.uvdata as uvd

o = optparse.OptionParser()
o.set_usage('firstcal_mwa.py [options] obsid')
a.scripting.add_standard_options(o,cal=True,pol=True)
o.add_option('--ubls', dest='ubls', default='57_61,57_62', help='Unique baselines to use, separated by commas (ex: 1_4,64_49).')
o.add_option('--ex_ubls', dest='ex_ubls', default='', help='Unique baselines to exclude, separated by commas (ex: 1_4,64_49).')
o.add_option('--bls', dest='bls', default='', help='Baselines to use, separated by commas (ex: 1_4,64_49).')
o.add_option('--ex_bls', dest='ex_bls', default='', help='Baselines to exclude, separated by commas (ex: 1_4,64_49).')
o.add_option('--ants', dest='ants', default='', help='Antennas to use, separated by commas (ex: 1,4,64,49).')
o.add_option('--ex_ants', dest='ex_ants', default='', help='Antennas to exclude, separated by commas (ex: 1,4,64,49).')
o.add_option('--outpath', default='/users/wl42/data/wl42/Nov2016EoR0/omni_sol/',help='Output path of solutions.')
o.add_option('--verbose', action='store_true', default=False, help='Turn on verbose.')
o.add_option('--ftype', dest='ftype', default='', type='string',
             help='Type of the input file, uvfits or fhd')
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
ex_ants_find = mp2cal.wyl.find_ex_ant(uv)
for a in ex_ants_find:
    if not a in ex_ants: ex_ants.append(a)
print '     ex_ants: ', ex_ants
reds = mp2cal.wyl.cal_reds_from_pos(antpos)
redbls = [bl for red in reds for bl in red]
print 'Number of redundant baselines:',len(redbls)
wrap_list = mp2cal.wyl.uv_wrap_fc(uv,redbls,pols=pols)
fqs = uv.freq_array[0]/1e9
del uv

#************************************************************************************************
def firstcal(data_wrap):
    pp = data_wrap['pol']
    p = pp[0]
    outname = opts.outpath + obsid + '.' + pp + '.fc.npz'
    #if os.path.exists(outname): raise IOError("File {0} already exists".format(outname))
    datpack = data_wrap['data']
    wgtpack = data_wrap['flag']
    flag_bls = []
    for bl in wgtpack.keys():
        wgt_data = np.logical_not(wgtpack[bl][pp])
        wgt_data = np.sum(wgt_data,axis=0)
        ind = np.where(wgt_data==0)
        if ind[0].size > 48: flag_bls.append(bl)
    print 'flagged baselines: ', flag_bls
    info = mp2cal.wyl.pos_to_info(antpos,pols=[p],fcal=True,ubls=ubls,ex_ubls=ex_ubls,bls=bls,ex_bls=ex_bls+flag_bls,ants=ants,ex_ants=ex_ants)
    wgtpack = {k : { qp : np.logical_not(wgtpack[k][qp]) for qp in wgtpack[k]} for k in wgtpack}
    fc = heracal.firstcal.FirstCal(datpack,wgtpack,fqs,info)
    print "     running firstcal"
    sols = fc.run(finetune=True,verbose=False,average=True,window='none')
    print('     Saving {0}'.format(outname))
    mp2cal.wyl.save_gains_fc(sols,fqs*1e9,outname)

#*****************************************************************************
par = Pool(2)
npzlist = par.map(firstcal, wrap_list)
par.close()
