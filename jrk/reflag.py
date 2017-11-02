import uvdata
import numpy as n
import optparse, sys, os

o = optparse.OptionParser()
opts,args = o.parse_args(sys.argv[1:])
print args

filelist = args
for file in filelist:
    print file.split(':')[0]
    dirty = uvdata.miriad.Miriad()
    dirty.read_miriad(file.split(':')[0])
    residual = uvdata.miriad.Miriad()
    residual.read_miriad(file.split(':')[1])
    residual.flag_array=dirty.flag_array
    #dirty.flag_array[:,:,:,:]=False
    #dirty.flag_array[:,0,152:203,:]=True
    residual.write_miriad(file.split(':')[1]+'F')
