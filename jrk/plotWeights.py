import aipy as a, numpy as n
import sys, os
from scipy.io import readsav
import pylab as plt
from optparse import OptionParser
import matplotlib.pyplot as mp

parser = OptionParser()
parser.add_option("-w", "--weights", dest="weights",help="UV Weights to FFT and plot.", metavar="FILE",action='store')
(options,args) = parser.parse_args()

def makeWeightsArray(X):
    savfile = readsav(X)
    wuv = savfile['weights_uv']
    wxy = n.fft.fft2(wuv)
    wuv = n.abs(wuv)
    wxy = n.abs(wxy)
    return wuv,wxy

uv,xy = makeWeightsArray(options.weights)
#mp.subplot(1,2,1)
#mp.imshow(n.log(uv))
#mp.colorbar()
#mp.subplot(1,2,2)
mp.imshow(n.fft.fftshift(xy))
mp.colorbar()
mp.savefig('weights.png')

    
