# edited version of firstcal in capo/omni/
#! /usr/bin/env python
import capo,aipy
import numpy as n
import sys,optparse

o = optparse.OptionParser()
aipy.scripting.add_standard_options(o,cal=True,pol=True)
o.add_option('--ubls', default='', help='Unique baselines to use, separated by commas (ex: 1_4,64_49).')
o.add_option('--ex_ants', default='', help='Antennas to exclude, separated by commas (ex: 1,4,64,49).')
o.add_option('--ftype', dest='ftype', default='', type='string',
             help='Type of the input file, .uvfits, or miriad, or fhd, to read fhd, simply type in the path/obsid')
o.add_option('--fcpath',dest='fcpath',default='',type='string',
             help='Path to save .npz files. Include final / in path.')
opts,args = o.parse_args(sys.argv[1:])

def flatten_reds(reds):
    freds = []
    for r in reds:
        freds += r
    return freds

def save_gains(s,f,pol,fn):
    s2 = {}
    for k,i in s.iteritems():
        s2[str(k)] = capo.omni.get_phase(f,i)
    print 'Saving %s_fcgains.%s.npz'%(fn,pol)
    n.savez('%s_fcgains.%s.npz'%(fn,pol),**s2)

def normalize_data(datadict):
    d = {}
    for key in datadict.keys():
        d[key] = datadict[key]/n.where(n.abs(datadict[key]) == 0., 1., n.abs(datadict[key]))
    return d 

files = {}
pols = opts.pol.split(',')
for filename in args:
    files[filename] = {}
    for p in pols:
        if opts.ftype == 'uvfits' or opts.ftype == 'miriad':
            fn = filename.split('.')
            fn[-2] = p
            files[filename][p] = '.'.join(fn)
        elif opts.ftype == 'fhd':
            obs = filename + '*'
            filelist = glob.glob(obs)
            if p == 'xx':
                filelist.remove(filename + '_vis_YY.sav')
            elif p == 'yy':
                filelist.remove(filename + '_vis_XX.sav')
            else:
                raise IOError('do not support cross pol')
            files[filename][p] = filelist
        else:
            raise IOError('invalid filetype, it should be miriad, uvfits, or fhd')

aa = aipy.cal.get_aa(opts.cal, n.array([.150]))
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
info = capo.omni.aa_pos_to_info(aa, pols=list(set(''.join(pols))), fcal=True, ubls=ubls, ex_ants=ex_ants, crosspols=pols)
reds = flatten_reds(info.get_reds())

print 'Number of redundant baselines:',len(reds)

ant_string =','.join(map(str,info.subsetant))
bl_string = ','.join(['_'.join(map(str,k)) for k in reds])

for f,filename in enumerate(args):
    
    file_group = files[filename]
    print 'Reading:'
    for key in file_group.keys(): print file_group[key]
    times,data,flags,ginfo,fqs = capo.wyl.uv_read([file_group[key] for key in file_group.keys()], filetype=opts.ftype, polstr=opts.pol, antstr='cross')
    fqs /= 1e9
    dd = {}
    for (i,j) in data.keys():
        for pp in file_group.keys():
            dd[(i,j)] = data[(i,j)][pp]

    dlys = n.fft.fftshift(n.fft.fftfreq(fqs.size, fqs[1]-fqs[0]))


#gets phase solutions per frequency.
    fc = capo.omni.FirstCal(dd,fqs,info)
    sols = fc.run(tune=True)
    fn = ''
    if opts.ftype == 'miriad' or 'uvfits':
        fn = opts.fcpath + '.'.join(filename.split('.'))[0:-2]
    elif opts.ftype == 'fhd':
        fn = opts.fcpath + filename

    pstr = ''.join(pols)
#Save solutions
    save_gains(sols,fqs, pstr, fn)




