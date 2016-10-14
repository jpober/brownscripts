# edited version of firstcal in capo/omni/
#! /usr/bin/env python
import capo.hex as hx, capo.wyl as wyl, capo.red as red, capo.omni as omni
import pylab as p, aipy as a
import sys,optparse,glob
import numpy as np
from IPython import embed

o = optparse.OptionParser()
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

def flatten_reds(reds):
    freds = []
    for r in reds:
        freds += r
    return freds


#hera info assuming a hex of 19 and 128 antennas
#aa = a.cal.get_aa(opts.cal, n.array([.150]))
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
info = omni.pos_to_info(_antpos, fcal=True, ubls=ubls, ex_ants=ex_ants)
reds = flatten_reds(info.get_reds())
#redstest = infotest.get_reds()#for plotting

print 'Number of redundant baselines:',len(reds)
#Read in data here.
ant_string =','.join(map(str,info.subsetant))
bl_string = ','.join(['_'.join(map(str,k)) for k in reds])
file_group = []
for fn in args:
    if opts.ftype == 'fhd':
        file_group.append(glob.glob(fn + '*'))
    elif opts.ftype == 'uvfits' or opts.ftype == 'miriad':
        file_group.append(fn)
    else:
        raise IOError('invalid filetype, it should be miriad, uvfits, or fhd')
npzlist = []
times, data, flags, ginfo, fqs = wyl.uv_read(file_group, filetype=opts.ftype, bl_str=bl_string, p_list=pols)
for pp in pols:
#arp.get_dict_of_uv_data(args, bl_string, opts.pol, verbose=True)
    datapack,wgtpack = {},{}
    for (i,j) in data.keys():
        datapack[(i,j)] = data[(i,j)][pp]
        wgtpack[(i,j)] = np.logical_not(flags[(i,j)][pp])
#    nfreq = datapack[datapack.keys()[0]].shape[1] #XXX less hacky than previous hardcode, but always safe?
    fqs = fqs/1e9 #  in GHz
    dlys = np.fft.fftshift(np.fft.fftfreq(fqs.size, np.diff(fqs)[0]))

#gets phase solutions per frequency.
    fc = omni.FirstCal(datapack,wgtpack,fqs,info)
    sols = fc.run(tune=True,verbose=opts.verbose,offset=True,plot=opts.plot)

#Save solutions
    if len(args)==1: filename=args[0]
    else: filename='fcgains.%s.npz'%pp #if averaging a bunch together of files together.
    if not opts.outpath is None:
        outname='%s/%s'%(opts.outpath,filename.split('/')[-1])
    else:
        outname='%s'%filename
#    embed()
    omni.save_gains_fc(sols,fqs, pp[0], outname, ubls=ubls, ex_ants=ex_ants)
    npzlist.append(outname+'.fc.npz')







