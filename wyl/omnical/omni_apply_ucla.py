#! /usr/bin/env python

import omnical, aipy, numpy, capo
import pickle, optparse, os, sys, glob
import uvdata.uvdata as uvd

### Options ###
o = optparse.OptionParser()
o.set_usage('omni_apply.py [options] *uvcRRE')
o.set_description(__doc__)
aipy.scripting.add_standard_options(o,pol=True)
o.add_option('--xtalk',dest='xtalk',default=False,action='store_true',
            help='Toggle: apply xtalk solutions to data. Default=False')
o.add_option('--omnipath',dest='omnipath',default='%s.npz',type='string',
            help='Format string (e.g. "path/%s.npz", where you actually type the "%s") which converts the input file name to the omnical npz path/file.')
o.add_option('--intype', dest='intype', default='', type='string',
             help='Type of the input file, .uvfits, or miriad, or fhd, to read fhd, simply type in the path/obsid')
o.add_option('--outtype', dest='outtype', default='uvfits', type='string',
             help='Type of the output file, .uvfits, or miriad, or fhd')
opts,args = o.parse_args(sys.argv[1:])


#File Dictionary
pols = opts.pol.split(',')
files = {}
for filename in args:
    files[filename] = {}
    if opts.intype == 'miriad':
        for p in pols:
            fn = filename.split('.')
            fn[-2] = p
            files[filename][p] = '.'.join(fn)
    elif opts.intype == 'uvfits':
        files[filename][filename] = filename
    elif opts.intype == 'fhd':
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

### Read Data and Solutions ###
for f,filename in enumerate(args):
    if opts.intype == 'miriad':
        if len(pols)>1:
            omnifile = opts.omnipath % '.'.join(filename.split('/')[-1].split('.')[0:3])
        else:
            omnifile = opts.omnipath % '.'.join(filename.split('/')[-1].split('.')[0:4])
    else:
        if len(pols)>1:
            omnifile = opts.omnipath % filename.split('/')[-1]
        else:
            omnifile = opts.omnipath % filename.split('/')[-1]+'.'+pols[0]
    print '   Omnical npz:', omnifile
    _,gains,_,xtalk = capo.omni.from_npz(omnifile) #loads npz outputs from omni_run

#    for ip,p in enumerate(pols):
    for key in files[filename].keys():
        print 'Reading', files[filename][key]
        if opts.outtype == 'uvfits':
            newfn = files[filename][key].split('.')
            newfn[-1] = 'O.uvfits'
            newfile = '.'.join(newfn)
        if os.path.exists(newfile):
            print '    %s exists.  Skipping...' % newfile
            continue

        uvi = uvd.UVData()
        if opts.intype == 'uvfits':
            uvi.read_uvfits(files[filename][key])
        elif opts.intype == 'miriad':
            uvi.read_miriad(files[filename][key])
        elif opts.intype == 'fhd':
            uvi.read_fhd(files[filename][key])

        Nblts = uvi.Nblts
        Nfreqs = uvi.Nfreqs
        Nbls = uvi.Nbls
        pollist = list(uvi.polarization_array)
        for p in pols:
            pid = pollist.index(aipy.miriad.str2pol[p])

            for ii in range(0,Nblts):
                a1 = uvi.ant_1_array[ii]
                a2 = uvi.ant_2_array[ii]
                p1,p2 = p
                ti = ii/Nbls
                if opts.xtalk:
                    try: uvi.data_array[:,0][:,:,pid][ii] -= xtalk[p][(a1,a2)]
                    except(KeyError):
                        try: uvi.data_array[:,0][:,:,pid][ii] -= xtalk[p][(a2,a1)].conj()
                        except(KeyError): pass
                try: uvi.data_array[:,0][:,:,pid][ii] /= gains[p1][a1][ti]
                except(KeyError): pass
                try: uvi.data_array[:,0][:,:,pid][ii] /= gains[p2][a2][ti].conj()
                except(KeyError): pass
#uvi.history = ''
        if opts.outtype == 'uvfits':
            print 'writing:' + newfile
            uvi.write_uvfits(newfile,spoof_nonessential=True)
            print 'saving ' + newfile
    


