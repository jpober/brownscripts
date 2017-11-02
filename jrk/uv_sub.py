import numpy as n
import pyuvdata
from numpy import fft
import aipy as a
import optparse, sys, os
import pylab as pl
from numpy.random import randn
from glob import glob

o = optparse.OptionParser()
o.add_option("--dirty", action='store', dest='dirty')
o.add_option("--model", action='store', dest='model')
o.add_option("-s","--suffix",dest="suffix",default="SP")
opts,args = o.parse_args(sys.argv[1:])

dirty = glob(opts.dirty)
dirty = n.sort(dirty)
models = glob(opts.model)
models = n.sort(models)
print models
#try:
for f in range(len(dirty)):
    mir = pyuvdata.miriad.Miriad()
    mir.read_miriad(dirty[f])

    model = pyuvdata.miriad.Miriad()
    model.read_miriad(models[f])

    for i in n.unique(mir.baseline_array):
        mir_id = mir.baseline_array==i#(mir.antnums_to_baseline(i,j)==mir.baseline_array)
        model_id = model.baseline_array==i#(model.antnums_to_baseline(i,j)==model.baseline_array)
        try:
            mir.data_array[mir_id,:,:,:] = mir.data_array[mir_id,:,:,:] - model.data_array[model_id,:,:,:]
        except:
            print 'Bad things keeps happening.'
            continue
#mir.antenna_positions = n.zeros((63,3))
    mir.write_miriad(dirty[f].split('HP')[0]+opts.suffix) #replace HP tp uv
    del(mir)
    del(model)
