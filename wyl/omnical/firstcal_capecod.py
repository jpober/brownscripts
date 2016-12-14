# edited version of firstcal in capo/omni/
#! /usr/bin/env python
import capo.hex as hx, capo.wyl as wyl, capo.red as red, capo.omni as omni
import pylab as p, aipy as a
import sys,optparse,glob
import numpy as np
from multiprocessing import Pool
from IPython import embed

o = optparse.OptionParser()
o.set_usage('firstcal_capecod.py [options] *zen.jds.pol.uv/obsid')
a.scripting.add_standard_options(o,cal=True,pol=True)
o.add_option('--ubls', default='', help='Unique baselines to use, separated by commas (ex: 1_4,64_49).')
o.add_option('--ex_ants', default='', help='Antennas to exclude, separated by commas (ex: 1,4,64,49).')
o.add_option('--outpath', default=None,help='Output path of solution npz files. Default will be the same directory as the data files.')
o.add_option('--plot', action='store_true', default=False, help='Turn on plotting in firstcal class.')
o.add_option('--verbose', action='store_true', default=False, help='Turn on verbose.')
o.add_option('--ftype', dest='ftype', default='', type='string',
             help='Type of the input file, .uvfits, or miriad, or fhd, to read fhd, simply type in the path/obsid')
opts,args = o.parse_args(sys.argv[1:])
print opts.plot
print opts.verbose

#*****************************************************************************
def flatten_reds(reds):
    freds = []
    for r in reds:
        freds += r
    return freds

def firstcal(datdict):
    datapack = datdict['dat']
    wgtpack = datdict['wgt']
    pp = datdict['pol']
    fqs = datdict['fqs']
    fqs = fqs/1e9 #  in GHz
    dlys = np.fft.fftshift(np.fft.fftfreq(fqs.size, np.diff(fqs)[0]))
    #gets phase solutions per frequency.
    info = omni.pos_to_info(datdict['antpos'], pols=[pp[0]], fcal=True, ubls=datdict['ubls'], ex_ants=datdict['ex_ants'])
    fc = omni.FirstCal(datapack,wgtpack,fqs,info)
    sols = fc.run(finetune=True,verbose=opts.verbose,plot=False,noclean=False,offset=False,average=True,window='none')
    
    #Save solutions
    filename = datdict['filename']
#    if len(args)==1: filename=args[0]
#    else: filename='fcgains.%s.npz'%pp #if averaging a bunch together of files together.
    if not datdict['outpath'] is None:
        outname='%s/%s'%(datdict['outpath'],filename.split('/')[-1])+'.'+pp
    else:
        outname='%s'%filename+'.'+pp
    #    embed()
    omni.save_gains_fc(sols,fqs, pp[0], outname, ubls=ubls, ex_ants=ex_ants)
    return (outname+'.fc.npz')
#*****************************************************************************

### Main ###
exec('from %s import antpos as _antpos'% opts.cal)
pols = opts.pol.split(',')
ex_ants = []
ubls = []
for a in opts.ex_ants.split(','):
    try: ex_ants.append(int(a))
    except: pass
for bl in opts.ubls.split(','):
    try:
        i,j = bl.split('_')
        ubls.append((int(i),int(j)))
    except: pass
print 'Excluding Antennas:',ex_ants
if len(ubls) != None: print 'Using Unique Baselines:',ubls
reds = omni.cal_reds_from_pos(_antpos)
reds = flatten_reds(reds)
#redstest = infotest.get_reds()#for plotting

print 'Number of redundant baselines:',len(reds)
#Read in data here.
bl_string = ','.join(['_'.join(map(str,k)) for k in reds])

file_group = []
if not len(args) == 1: raise IOError('Do not support multiple files.')
for fn in args:
    if opts.ftype == 'fhd':
        file_group.append(glob.glob(fn + '*'))
    elif opts.ftype == 'uvfits':
        file_group.append(fn+'.uvfits')
    elif opts.ftype == 'miriad':
        if len(pols) == 1:
            file_group.append(fn)
        else:
            for pp in pols:
                fnlist = fn.split('.')
                fnlist[-2] = pp
                file_group.append('.'.join(fnlist))
    else:
        raise IOError('invalid filetype, it should be miriad, uvfits, or fhd')
npzlist = []
times, data, flags, ginfo, fqs, exa = wyl.uv_read_fc(file_group, filetype=opts.ftype, bl_str=bl_string, p_list=pols)
if len(exa) > 0:
    for ii in exa:
        if not ii in ex_ants: ex_ants.append(ii)
print 'Excluding Antennas:',ex_ants
arglist = []
if opts.ftype == 'miriad': fn = '.'.join(args[0].split('.')[0:-2])
else: fn = args[0]
for pp in pols:
    dict = {}
    dict['pol'] = pp
    dict['fqs'] = fqs
    datapack,wgtpack = {},{}
    for (i,j) in data.keys():
        datapack[(i,j)] = data[(i,j)][pp]
        wgtpack[(i,j)] = np.logical_not(flags[(i,j)][pp])
    dict['dat'] = datapack
    dict['wgt'] = wgtpack
    dict['outpath'] = opts.outpath
    dict['filename'] = fn
    dict['ubls'] = ubls
    dict['ex_ants'] = ex_ants
    dict['antpos'] = _antpos
    arglist.append(dict)

print "  Start Parallelism:"
par = Pool(2)
npzlist = par.map(firstcal, arglist)
par.close()

#npzlist.append(outname+'.fc.npz')
#omni.fc_gains_to_fits(npzlist,fn)







