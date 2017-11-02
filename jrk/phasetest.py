from scipy.io import readsav
import uvdata
import numpy as n
import optparse, sys, os

o = optparse.OptionParser()
opts,args = o.parse_args(sys.argv[1:])
print args

filelist = args
for filename in filelist:
    base = filename.split('/')[-1]
    phase_me = uvdata.miriad.Miriad()
    phase_me.read_miriad(filelist[0])
    phase_me.phase(time=n.sort(phase_me.time_array)[0])
    phase_me.antenna_positions = n.zeros((63,3))
    phase_me.write_miriad(base+'n')



