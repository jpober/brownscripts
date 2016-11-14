import numpy as n
import uvdata
from numpy import fft
import aipy as a
import optparse, sys, os
import pylab as pl
from numpy.random import randn

o = optparse.OptionParser()
o.add_option("-s","--suffix",dest="suffix",default="SP")
opts,args = o.parse_args(sys.argv[1:])
print args
print opts.suffix
t=0
for files in args:
    try:
        if t == 0:
            mir = uvdata.miriad.Miriad()
            mir.read_miriad(files)
        else:
            model = uvdata.miriad.Miriad()
            model.read_miriad(files)
    except:
        print 'Something bad happened.'
        pass
    t+=1

for i in range(0,64):
    for j in range(0,64):
        mir_id = (mir.antnums_to_baseline(i,j)==mir.baseline_array)
        model_id = (model.antnums_to_baseline(i,j)==model.baseline_array)
        try:
            mir.data_array[mir_id,:,:,:] = mir.data_array[mir_id,:,:,:] - model.data_array[model_id,:,:,:]
        except:
            print 'Bad things keeps happening.'
            continue
mir.antenna_positions = n.zeros((63,3))
mir.write_miriad(args[0].split('HP')[0]+opts.suffix) #replace HP tp uv
    
