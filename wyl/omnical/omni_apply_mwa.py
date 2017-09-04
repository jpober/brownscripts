#! /usr/bin/env python
# Do not support miriad

import numpy as np
import aipy, mp2cal
import pickle, optparse, os, sys, glob
import pyuvdata.uvdata as uvd
from astropy.io import fits
from scipy.io.idl import readsav

### Options ###
o = optparse.OptionParser()
o.set_usage('omni_apply_fhd.py [options] obsid(do not include .uvfits) or zen.jds.pol.uv')
o.set_description(__doc__)
aipy.scripting.add_standard_options(o,pol=True,cal=True)
o.add_option('--xtalk',dest='xtalk',default=False,action='store_true',
            help='Toggle: apply xtalk solutions to data. Default=False')
o.add_option('--bpfit',dest='bpfit',default=False,action='store_true',
             help='Toggle: do a global bp fit to sols. Default=False')
o.add_option('--polyfit',dest='polyfit',default=False,action='store_true',
             help='Toggle: do a polyfit to sols over the band. Default=False')
o.add_option('--omnipath',dest='omnipath',default='',type='string', help='Path to load omnical solution files. Include final / in path.')

o.add_option('--outtype', dest='outtype', default='uvfits', type='string',
             help='Type of the output file, .uvfits, or miriad, or fhd')
o.add_option('--intype', dest='intype', default=None, type='string',
             help='Type of the input file, .uvfits or fhd')
o.add_option('--metafits', dest='metafits', default='/users/wl42/data/wl42/Nov2016EoR0/', type='string',
             help='path to metafits files')
o.add_option('--fhdpath', dest='fhdpath', default='/users/wl42/data/wl42/FHD_out/fhd_PhaseII_Longrun_EoR0/', type='string',
             help='path to fhd dir for fhd output visibilities if ftype is fhd.')
o.add_option('--appfhd',dest='appfhd',default=False,action='store_true',
             help='Toggle: apply FHD solutions to non-hex tiles. Default=False')
o.add_option('--ave',dest='ave',default=False,action='store_true',
             help='Toggle: apply averaged calibration solution. Default=False')
o.add_option('--nodegen', dest='nodegen', default=False, action='store_true',
             help='Toggle: Do not project degeneracy')
opts,args = o.parse_args(sys.argv[1:])

delays = {
'0,5,10,15,1,6,11,16,2,7,12,17,3,8,13,18':-5,
'0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15':-4,
'0,3,6,9,0,3,6,9,0,3,6,9,0,3,6,9':-3,
'0,2,4,6,0,2,4,6,0,2,4,6,0,2,4,6':-2,
'0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3':-1,
'0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0':0,
'3,2,1,0,3,2,1,0,3,2,1,0,3,2,1,0':1,
'6,4,2,0,6,4,2,0,6,4,2,0,6,4,2,0':2,
'9,6,3,0,9,6,3,0,9,6,3,0,9,6,3,0':3,
'12,8,4,0,13,9,5,1,14,10,6,2,15,11,7,3':4,
'15,10,5,0,16,11,6,1,17,12,7,2,18,13,8,3':5,
}


#File Dictionary
pols = opts.pol.split(',')
files = {}

if not len(args) == 1: raise IOError('Do not support multiple files.')
obsid = args[0]

#start processing
    #create an out put filename
if opts.outtype == 'uvfits':
    if opts.intype == 'fhd': suffix = 'FO'
    else: suffix = 'O'
    if opts.ave: suffix = suffix + 'A'
    if opts.bpfit:
        suffix = suffix + 'B'
    if opts.polyfit:
        suffix = suffix + 'P'
    newfile = obsid + '_' + suffix + '.uvfits'
if os.path.exists(newfile): raise IOError('   %s exists.  Skipping...' % newfile)

    #read in the file
print '  Reading', obsid
uvi = uvd.UVData()
if opts.intype == 'fhd':
    uvi.read_fhd(glob.glob(opts.fhdpath+'/vis_data/'+obsid+'*')+glob.glob(opts.fhdpath+'/metadata/'+obsid+'*'),use_model=False,run_check=False,run_check_acceptability=False)
elif opts.intype == 'uvfits':
    uvi.read_uvfits(obsid+'.uvfits',run_check=False,run_check_acceptability=False)
Nblts = uvi.Nblts
Nfreqs = uvi.Nfreqs
Nbls = uvi.Nbls
pollist = uvi.polarization_array
freqs = uvi.freq_array[0]

gfhd = mp2cal.wyl.load_gains_fhd(opts.fhdpath+'calibration/'+obsid+'_cal.sav')

    #find npz for each pol, then apply
for ip,p in enumerate(pols):
    pid = np.where(pollist == aipy.miriad.str2pol[p])[0][0]
    omnifile_ave = ''
    if opts.ave:
        day = int(obsid)/86400
        hdu = fits.open(opts.metafits+obsid+'.metafits')
        pointing = delays[hdu[0].header['DELAYS']]
        omnifile_ave = opts.omnipath + 'omniave' + '_' + str(day) + '_' + str(pointing) + '.' + p + '.npz'
    omnifile = opts.omnipath + obsid.split('/')[-1]+'.'+p+'.omni.npz'
    print '  Reading and applying:', omnifile, omnifile_ave
    if opts.ave:
        if opts.xtalk:
            _,gains,_,_ = mp2cal.wyl.load_gains_omni(omnifile_ave)
            _,_,_,xtalk = mp2cal.wyl.load_gains_omni(omnifile)
        else:
            xtalk = {}
            gains = mp2cal.wyl.quick_load_gains(omnifile_ave)
    else:
        if opts.xtalk:
            _,gains,_,xtalk = mp2cal.wyl.load_gains_omni(omnifile) #loads npz outputs from omni_run
        else:
            xtalk = {}
            gains = mp2cal.wyl.quick_load_gains(omnifile)
#********************** if project degeneracy ***********************************
    if not opts.nodegen:
        print "   Project Degeneracy"
        exec('from %s import *'% opts.cal)
        ref = min(gains[p[0]].keys())
        if opts.intype == 'fhd':
            gains = mp2cal.wyl.degen_project_FO(gains,antpos)
        elif opts.intype == 'uvfits':
            gains = mp2cal.wyl.degen_project_OF(gains,gfhd,antpos,EastHex,SouthHex)
#********************** if choose to make sols smooth ***************************
    if opts.bpfit:
        print '   bandpass fitting'
        exec('from %s import tile_info'% opts.cal)
        auto = {}
        bi = 128*uvi.ant_1_array+uvi.ant_2_array
        for a in gains[p[0]].keys():
            auto_inds = np.where(bi==128*a+a)
            auto_corr = uvi.data_array[auto_inds][:,0][:,:,pid].real
            auto[a] = np.mean(np.sqrt(auto_corr)[1:-2],axis=0)
            auto[a] /= np.mean(auto[a])
        gains = mp2cal.wyl.mwa_bandpass_fit(gains,auto,tile_info)
    if opts.polyfit:
        print '   polyfitting'
        gains = mp2cal.wyl.poly_bandpass_fit(gains)
    ex_ants = []
#*********************************************************************************************
    if opts.appfhd:
        for a in gfhd[p[0]].keys():
            if np.isnan(np.mean(gfhd[p[0]][a])):
                ex_ants.append(a)
                continue
            if a > 56: continue
            gains[p[0]][a] = gfhd[p[0]][a]
    for ii in range(0,Nblts):
        a1 = uvi.ant_1_array[ii]
        a2 = uvi.ant_2_array[ii]
        if a1 in ex_ants or a2 in ex_ants:
            uvi.flag_array[:,0][:,:,pid][ii] = True
            continue
        p1,p2 = p
        if opts.xtalk:
            try: uvi.data_array[:,0][:,:,pid][ii] -= xtalk[p][(a1,a2)]
            except(KeyError):
                try: uvi.data_array[:,0][:,:,pid][ii] -= xtalk[p][(a2,a1)].conj()
                except(KeyError): pass
        try:
            fuse = np.where(gains[p1][a1]!=0)
            fnot = np.where(gains[p1][a1]==0)
            uvi.data_array[:,0][:,:,pid][ii][fuse] /= gains[p1][a1][fuse]
            uvi.data_array[:,0][:,:,pid][ii][fnot] *= 0
            uvi.flag_array[:,0][:,:,pid][ii][fnot] = True
        except(KeyError): pass
        try:
            fuse = np.where(gains[p2][a2]!=0)
            fnot = np.where(gains[p2][a2]==0)
            uvi.data_array[:,0][:,:,pid][ii][fuse] /= gains[p2][a2][fuse].conj()
            uvi.data_array[:,0][:,:,pid][ii][fnot] *= 0
            uvi.flag_array[:,0][:,:,pid][ii][fnot] = True
        except(KeyError): pass

    #write file
#uvi.history = ''
if opts.outtype == 'uvfits':
    print 'writing:' + newfile
    uvi.write_uvfits(newfile,spoof_nonessential=True)
    print 'saving ' + newfile


