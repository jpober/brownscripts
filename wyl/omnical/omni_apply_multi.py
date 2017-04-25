#! /usr/bin/env python
# Do not support miriad

import numpy as np
import omnical, aipy, capo
import pickle, optparse, os, sys, glob
import pyuvdata.uvdata as uvd
from astropy.io import fits

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
o.add_option('--omnipath',dest='omnipath',default='%s.npz',type='string',
            help='Format string (e.g. "path/%s.npz", where you actually type the "%s") which converts the input file name to the omnical npz path/file.')
o.add_option('--npz',dest='npz',default=None,type='string',
             help='specify npz file names for gain solutions, (format: path/name, without .pol.npz), otherwise read gain solutions from npz according to obsid and omnipath')
o.add_option('--outtype', dest='outtype', default='uvfits', type='string',
             help='Type of the output file, .uvfits, or miriad, or fhd')
o.add_option('--intype', dest='intype', default=None, type='string',
             help='Type of the input file, .uvfits or fhd')
o.add_option('--instru', dest='instru', default='mwa', type='string',
             help='instrument type. Default=mwa')
o.add_option('--flag_bls',dest='flag_bls',default=False,action='store_true',
             help='Toggle: Flag baselines which are excluded by omnical. Default=False')
o.add_option('--metafits', dest='metafits', default='/users/wl42/data/wl42/Nov2016EoR0/', type='string',
             help='path to metafits files')
o.add_option('--fhdpath', dest='fhdpath', default='/users/wl42/data/wl42/FHD_out/fhd_PhaseII_Longrun_EoR0/', type='string',
             help='path to fhd dir for fhd output visibilities if ftype is fhd.')
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

#create a dictionary of file lists
for filename in args:
    if opts.intype == 'fhd':
        files[filename] = []
        filelist = glob.glob(opts.fhdpath+'/vis_data/'+filename+'*')+glob.glob(opts.fhdpath+'/metadata/'+filename+'*')
        files[filename] = filelist
    elif opts.intype == 'uvfits':
        files[filename] = filename + '.uvfits'
    elif opts.intype == 'miriad':
        files[filename] = {}
        for p in pols:
            fn = filename.split('.')
            fn[-2] = p
            files[filename][p] = '.'.join(fn)
    else:
        raise IOError('invalid filetype, it should be miriad, uvfits, or fhd')

#start processing
for f,filename in enumerate(args):
    
    #create an out put filename
    if opts.outtype == 'uvfits':
        if opts.intype == 'fhd': suffix = 'FO'
        else: suffix = 'O'
        if opts.bpfit:
            suffix = suffix + 'B'
        if opts.polyfit:
            suffix = suffix + 'P'
        newfile = filename + '_' + suffix + '.uvfits'
    if os.path.exists(newfile):
        print '    %s exists.  Skipping...' % newfile
        continue

    #read in the file
    print '  Reading', files[filename]
    uvi = uvd.UVData()
    if opts.intype == 'fhd':
        uvi.read_fhd(files[filename])
    elif opts.intype == 'uvfits':
        uvi.read_uvfits(files[filename])
    Nblts = uvi.Nblts
    Nfreqs = uvi.Nfreqs
    Nbls = uvi.Nbls
    pollist = uvi.polarization_array
    freqs = uvi.freq_array[0]

    #find npz for each pol, then apply
    for ip,p in enumerate(pols):
        pid = np.where(pollist == aipy.miriad.str2pol[p])[0][0]
        omnifile_ave = ''
        if not opts.npz == None:
            if opts.instru == 'mwa':
                day = int(filename)/86400
                hdu = fits.open(opts.metafits+filename+'.metafits')
                pointing = delays[hdu[0].header['DELAYS']]
                omnifile_ave = opts.npz + '_' + str(day) + '_' + str(pointing) + '.' + p + '.npz'
            else: omnifile_ave = opts.npz + '.' + p + '.npz'
        omnifile = opts.omnipath % (filename.split('/')[-1]+'.'+p)
        print '  Reading and applying:', omnifile, omnifile_ave
        if not opts.npz == None:
            _,gains,_,_ = capo.omni.from_npz(omnifile_ave)
            meta,_,_,xtalk = capo.omni.from_npz(omnifile)
        else:
            meta,gains,_,xtalk = capo.omni.from_npz(omnifile) #loads npz outputs from omni_run
        if opts.flag_bls:
            ex_bls = []
            for ii in range(meta['ex_bls'].shape[0]):
                ex_bls.append(tuple(meta['ex_bls'][ii]))
            print '   bls to flag:', ex_bls
#********************** if choose to make sols smooth ***************************
        if opts.bpfit and opts.instru == 'mwa':
            print '   bandpass fitting'
            exec('from %s import tile_info'% opts.cal)
            auto = {}
            bi = 128*uvi.ant_1_array+uvi.ant_2_array
            for a in gains[p[0]].keys():
                auto_inds = np.where(bi==128*a+a)
                auto_corr = uvi.data_array[auto_inds][:,0][:,:,pid].real
                auto[a] = np.mean(np.sqrt(auto_corr)[1:-2],axis=0)
                auto[a] /= np.mean(auto[a])
            gains = capo.wyl.mwa_bandpass_fit(gains,auto,tile_info)
        if opts.polyfit:
            print '   polyfitting'
            gains = capo.wyl.poly_bandpass_fit(gains)
#*********************************************************************************************
        for ii in range(0,Nblts):
            a1 = uvi.ant_1_array[ii]
            a2 = uvi.ant_2_array[ii]
            if (a1,a2) in ex_bls or (a2,a1) in ex_bls:
                if opts.flag_bls: uvi.flag_array[:,0][:,:,pid][ii] = True
            p1,p2 = p
#            ti = ii/Nbls
                #for jj in range(0,Nfreqs):
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


