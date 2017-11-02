#! /usr/bin/env python
import numpy as n
import aipy as a
import sys, os, optparse 

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
o.add_option('--plot', action='store_true', 
            help='plot bandpass.')
opts, args = o.parse_args(sys.argv[1:])


uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)
freqs = aa.get_afreqs() * 1e3
print freqs
#This was from 'zaki_passband_9th_order.npz'
bpfile = n.load('zaki_passband_9th_order.npz')
gianni_bandpass_coeffs = bpfile['coeff']

print gianni_bandpass_coeffs
gianni_bandpass = n.polyval(gianni_bandpass_coeffs, freqs)
if opts.plot:
    import pylab as p
    p.plot(freqs, gianni_bandpass)
    p.show()




def mfunc(uv, p, d, f):
    uvw,t,(ij) = p
    pol = a.miriad.pol2str[uv['pol']]
    aa.set_active_pol(pol)
    d = d*gianni_bandpass
    return p, d, f
     

for infile in args:
    outfile = infile+'G'
    print infile, '-->', outfile
    if os.path.exists(outfile):
        print 'File exists, skipping...'
        continue
    
    uvi = a.miriad.UV(infile)
    uvo = a.miriad.UV(outfile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc, raw=True, append2hist='APPLY-GIANNI-BANDPASS: ' + ' '.join(sys.argv) + '\n')


