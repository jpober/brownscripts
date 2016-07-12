import aipy as a, numpy as n
import sys, os
from scipy.io import readsav
import pylab as plt
from optparse import OptionParser
import matplotlib.pyplot as mp

parser = OptionParser()
parser.add_option("-d", "--dirty", dest="dirty",help="Dirty Visibility from fhd to plot", metavar="FILE",action='store')
parser.add_option("-m", "--model", dest="model",help="Model Visibility from fhd to plot", metavar="FILE",action='store')
parser.add_option("-b", "--baseline", dest="baseline",help="Choose the baseline to plot", metavar="FILE",action='store')
(options,args) = parser.parse_args()
#filename = args
#print options.filename
#savfile = readsav(options.filename)

def makeVisArray(X):
    savfile = readsav(X)
    for m in X.split('_'):
        if m == 'model':
            vis_array = savfile['vis_model_ptr']
            break
    if m != 'model':
        vis_array = savfile['vis_ptr']
    obs = savfile['obs']
    times = obs['baseline_info'][0]['JDATE'][0]
    ants1 = obs['baseline_info'][0]['TILE_A'][0]
    ants2 = obs['baseline_info'][0]['TILE_B'][0]
    ntimes = len(times)
    nbls = obs['NBASELINES'][0]
    time_order = n.argsort(times)
    times = times[time_order]
    _vis_array = n.zeros_like(vis_array)
    for ntime in xrange(ntimes):
        correct_time = time_order[ntime]
        _vis_array[ntime*nbls:(ntime+1)*nbls] = vis_array[correct_time*nbls:(correct_time+1)*nbls]
    vis_array = _vis_array
    return vis_array
#for m in options.filename.split('_'):
#    print m
#    if m == 'model':
#        vis_array = savfile['vis_model_ptr']
#        break
#if m != 'model':
#    vis_array = savfile['vis_ptr']
#obs = savfile['obs']
#times = obs['baseline_info'][0]['JDATE'][0]
#ants1 = obs['baseline_info'][0]['TILE_A'][0]
#ants2 = obs['baseline_info'][0]['TILE_B'][0]
#ntimes = len(times)
#nbls = obs['NBASELINES'][0]

#reorder vis_array
#time_order = n.argsort(times)
#times = times[time_order]
#_vis_array = n.zeros_like(vis_array)
#for ntime in xrange(ntimes):
#    correct_time = time_order[ntime]
#    _vis_array[ntime*nbls:(ntime+1)*nbls] = vis_array[correct_time*nbls:(correct_time+1)*nbls]
#vis_array = _vis_array
#(_,nfreqs) = vis_array.shape
#print _,nfreqs
visDirty = n.absolute(makeVisArray(options.dirty))
visModel = n.absolute(makeVisArray(options.model))*100
#visDirty = visDirty/n.max(visDirty)
#visModel = visModel/n.max(visModel)
print n.max(visDirty),n.max(visModel)
(_,nfreqs) = visDirty.shape
freqtime = n.zeros((14,nfreqs))
#visRes = visDirty-visModel
nbl = 63*62/2
#for i in range(0,13):
#    for j in range(0,13):
#freqtime[:,:] = n.abs(visRes[14*int(options.baseline):14*(int(options.baseline) + 1),:])
bsl = int(options.baseline)
mp.figure()

mp.subplot(3,1,1)
mp.imshow(visDirty[0:200,:])
mp.colorbar()
mp.subplot(3,1,2)
mp.imshow(visModel[0:200,:])
mp.colorbar()
mp.subplot(3,1,3)
mp.imshow(visDirty[0:200,:] - visModel[0:200,:])
mp.colorbar()
mp.savefig('freqtime_%i.png' %bsl)

    
