#! /usr/bin/env python
"""
Add or subtract the data in the second UV file from the data in the 
first UV file.  They'll need to have corresponding integrations/data order.
"""

import aipy as a, numpy as n, sys, optparse, os

o = optparse.OptionParser()
o.set_usage('uv_addsub.py [options] file1.uv file2.uv')
o.set_description(__doc__)
o.add_option('--div', dest='div', action='store_true',
    help='Subtract the data instead of adding.')
opts,args = o.parse_args(sys.argv[1:])

#assert(len(args) == 2)
uv1 = a.miriad.UV(args[0])
uv2 = a.miriad.UV(args[1])

def mfunc(uv, p, d, f):
    p2,d2,f2 = uv2.read(raw=True)
    f = n.logical_or(f,f2)
    if opts.div:
        D = n.fft.fft(d)
        D2 = n.fft.fft(d2)
        D2D = (D*n.conj(D2))/(n.conj(D2)*D2)
        d = n.fft.ifft(D2D)
        #d = d*n.conj(d2)/(n.conj(d2)*d2)
    else: 
        d = d + d2
    return p, n.where(f, 0, d), f

if opts.div: filename = args[0] + 'd'
else: filename = args[0] + 'a'
if os.path.exists(filename):
    print 'File exists: skipping'
    sys.exit(0)
print args[0], '->', filename
uvo = a.miriad.UV(filename, status='new')
uvo.init_from_uv(uv1)
uvo.pipe(uv1, mfunc=mfunc, raw=True, 
    append2hist='UVADDSUB: --div=%s file=%s\n' % (opts.div, args[1]))
