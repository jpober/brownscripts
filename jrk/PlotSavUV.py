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
parser.add_option("-f", "--flags", dest="flags",help="Choose the baseline to plot",action='store')
(options,args) = parser.parse_args()
#filename = args
#print options.filename
#savfile = readsav(options.filename)

def makeVisArray(X):
    #savfile = []
    #vis_array = []
    savfile = readsav(X)
    for m in X.split('_'):
        if m == 'model':
            vis_array = savfile['vis_model_ptr']
            print 'Model found!'
            break
    if m != 'model':
        vis_array = savfile['vis_ptr']
    try:
        flag = readsav(X.split('_')[0]+'_flags.sav')
        X.split('_')[0]+'_flags.sav'
    except IOError:
        flag = readsav(('_').join(X.split('_')[0:3])+'_flags.sav')
        ('_').join(X.split('_')[0:3])+'_flags.sav'
    print flag.keys()
    flag = flag['flag_arr'][0]
    #vis_array[flag==1] = 0
    #vis_array[flag==1] = 0
    print n.max(vis_array)
    #vis_array[flag==-1]=0
    obs = savfile['obs']
    times = obs['baseline_info'][0]['JDATE'][0]
    #ants1 = obs['baseline_info'][0]['TILE_A'][0]
    #ants2 = obs['baseline_info'][0]['TILE_B'][0]
    #print n.where(vis_array[5000,:]==0)[0].shape
    ntimes = len(times)
    nbls = obs['NBASELINES'][0]
    time_order = n.argsort(times)
    ant1,ant2 = options.baseline.split('_')
    ind1 = (obs[0]['baseline_info']['tile_a'][0]==int(ant1)).astype(int)
    ind2 = (obs[0]['baseline_info']['tile_b'][0]==int(ant2)).astype(int)
    intersect = ((ind1+ind2)/2).astype(bool)
    print intersect.max()
    bsl_array = vis_array[intersect]
    print bsl_array.shape
    bsl_array = bsl_array[time_order]
    return bsl_array

visDirty = makeVisArray(options.dirty)
visModel = makeVisArray(options.model)
(_,nfreqs) = visDirty.shape
nbl = 63*62/2
#bsl = int(options.baseline)

mp.figure()
ampDirty = n.abs(visDirty)#n.sqrt(n.real(visDirty*n.conj(visDirty)))
ampModel = n.abs(visModel)#n.sqrt(n.real(visModel*n.conj(visModel)))
#scaleFactor = n.mean(ampModel)/n.mean(ampDirty)
#print 'Scale Factor: ',scaleFactor
ampRes = ampDirty/ampModel#n.abs(visDirty-visModel)
#print n.linspace(100,200,203)[n.where(ampRes[0,:]>0)[0]]
#n.sqrt(n.real(ampRes*n.conj(ampRes)))
#z = n.where(ampRes<0)
#ampRes[z] = 0
#print n.min(n.log10(ampDirty)[n.log10(ampDirty)>0])
print 'Average (Dirty):(Model):(Residual)  ',n.mean(ampDirty),':',n.mean(ampModel),':',n.mean(ampRes)
#mp.title(options.baseline)
mp.subplot(1,3,1)
mp.title('Dirty')
#mp.imshow(ampDirty,extent=(100,200,0,14),vmin=n.min(ampModel),vmax=n.max(ampModel))
#mp.imshow(ampDirty,extent=(100,200,0,14))
mp.imshow(n.log10(ampDirty),aspect='auto',interpolation='nearest',extent=(100,200,0,14))
mp.colorbar()
mp.subplot(1,3,2)
mp.title('Model')
mp.imshow(n.log10(ampModel),aspect='auto',interpolation='nearest',extent=(100,200,0,14),vmin=n.min(n.log10(ampDirty)[n.log10(ampDirty)>0]),vmax=n.max(n.log10(ampDirty)))
mp.colorbar()
mp.subplot(1,3,3)
mp.title('Residual')
mp.imshow(n.log10(ampRes),aspect='auto',interpolation='nearest',extent=(100,200,0,14),vmin=n.min(n.log10(ampDirty)[n.log10(ampDirty)>0]),vmax=n.max(n.log10(ampDirty)))
mp.colorbar()
#mp.colorbar(orientation='horizontal')
#mp.xlabel('Freq.(Mhz)')
#mp.ylabel('Time Sample')
#mp.colorbar()
#mp.subplot(1,3,3)
#mp.title('Model')
#mp.imshow(ampModel,extent=(100,200,0,14))
#mp.colorbar(orientation='horizontal')
#mp.xlabel('Freq.(Mhz)')
#mp.ylabel('Time Sample')
#mp.colorbar()
#mp.subplot(3,1,3)
#mp.title('Residual')
#mp.imshow(ampRes,extent=(100,200,0,14),vmin=n.min(ampModel),vmax=n.max(ampModel))
#mp.imshow(ampRes,extent=(100,200,0,14))
#mp.colorbar(orientation='horizontal')
#mp.xlabel('Freq.(Mhz)')
#mp.ylabel('Time Sample')
#mp.colorbar()
mp.show()
#mp.savefig('freqtime_%i.png' %bsl)

    
