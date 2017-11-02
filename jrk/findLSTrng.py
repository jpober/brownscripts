import aipy as a
import capo.arp as arp
import capo.zsa as zsa
import numpy as n, pylab as p
import sys, os, optparse
#### Currently set for a 8.5 Hr LST window to check for consistency ###
fileList=' '
for arg in sys.argv[1:]:
    uv = a.miriad.UV(arg)
    lstTime=24.0*uv['lst']/(2*n.pi) + 14.3/60.0
    if lstTime >= 0 and lstTime <= 8.8:
        fileList=arg+' '+fileList
    del(uv)
        
print fileList

