#! /usr/bin/env python
# Do not support miriad

import omnical, aipy, numpy, capo
import pickle, optparse, os, sys, glob
import uvdata.uvdata as uvd

### Options ###
o = optparse.OptionParser()
o.set_usage('omni_apply_fhd.py [options] obsid (do not include .uvfits)')
o.set_description(__doc__)
aipy.scripting.add_standard_options(o,pol=True)
o.add_option('--xtalk',dest='xtalk',default=False,action='store_true',
            help='Toggle: apply xtalk solutions to data. Default=False')
o.add_option('--omnipath',dest='omnipath',default='%s.npz',type='string',
            help='Format string (e.g. "path/%s.npz", where you actually type the "%s") which converts the input file name to the omnical npz path/file.')
o.add_option('--outtype', dest='outtype', default='uvfits', type='string',
             help='Type of the output file, .uvfits, or miriad, or fhd')
o.add_option('--intype', dest='intype', default=None, type='string',
             help='Type of the input file, .uvfits or fhd')
opts,args = o.parse_args(sys.argv[1:])


#File Dictionary
pols = opts.pol.split(',')
files = {}

#create a dictionary of file lists
for filename in args:
    if opts.intype == 'fhd':
        files[filename] = []
        obs = filename + '*'
        filelist = glob.glob(obs)
        files[filename] = filelist
    elif opts.intype == 'uvfits':
        files[filename] = filename + '.uvfits'
    else:
        raise IOError('invalid filetype, it should be uvfits, or fhd')

#start processing
for f,filename in enumerate(args):
    
    #create an out put filename
    if opts.outtype == 'uvfits':
        newfile = filename + '_O.uvfits'
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
    pollist = list(uvi.polarization_array)

    #find npz for each pol, then apply
    for ip,p in enumerate(pols):
        omnifile = opts.omnipath % (filename.split('/')[-1]+'.'+p)
        print '  Reading and applying:', omnifile
        _,gains,_,xtalk = capo.omni.from_npz(omnifile) #loads npz outputs from omni_run
        pid = pollist.index(aipy.miriad.str2pol[p])
        for ii in range(0,Nblts):
            a1 = uvi.ant_1_array[ii]
            a2 = uvi.ant_2_array[ii]
            p1,p2 = p
            ti = ii/Nbls
                #for jj in range(0,Nfreqs):
            if opts.xtalk:
                try: uvi.data_array[:,0][:,:,pid][ii] -= xtalk[p][(a1,a2)]
                except(KeyError):
                    try: uvi.data_array[:,0][:,:,pid][ii] -= xtalk[p][(a2,a1)].conj()
                    except(KeyError): pass
            try: uvi.data_array[:,0][:,:,pid][ii] /= gains[p1][a1][ti]
            except(KeyError): pass
            try: uvi.data_array[:,0][:,:,pid][ii] /= gains[p2][a2][ti].conj()
            except(KeyError): pass

    #write file
#uvi.history = ''
    if opts.outtype == 'uvfits':
        print 'writing:' + newfile
        uvi.write_uvfits(newfile,spoof_nonessential=True)
        print 'saving ' + newfile


