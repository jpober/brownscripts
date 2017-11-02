import numpy as n
import uvdata
from numpy import fft
import aipy as a
import optparse, sys, os
import pylab as pl
from numpy.random import randn

o = optparse.OptionParser()
o.add_option("--dirty", action='store', dest='dirty')
o.add_option("--residual", action='store', dest='residual')
o.add_option("-s","--suffix",dest="suffix",default="r")
opts,args = o.parse_args(sys.argv[1:])
print args
print opts.suffix
print opts.dirty
print opts.residual
#t=0
#for files in args:

#try:
mir = uvdata.miriad.Miriad()
mir.read_miriad(opts.dirty)

res = uvdata.miriad.Miriad()
res.read_miriad(opts.residual)
#    except:
#        print 'Something bad happened.'
#        pass
#    t+=1
window =  a.dsp.gen_window(20,window='blackman-harris').astype(complex)
ratArray = n.zeros_like(mir.data_array,dtype=complex)
for i in range(0,64):
    for j in range(0,64):
        mir_id = (mir.antnums_to_baseline(i,j)==mir.baseline_array)
        res_id = (res.antnums_to_baseline(i,j)==res.baseline_array)
        #try:
        D = n.fft.fftshift(n.fft.fft(window*mir.data_array[mir_id,0,95:115,0],axis=1),axes=1)
        R = n.fft.fftshift(n.fft.fft(window*res.data_array[res_id,0,95:115,0],axis=1),axes=1)
        
        RATIO = R*n.conj(D)/(n.conj(D)*D)
        ratArray[mir_id,0,95:115,0] = RATIO
        ratio = n.fft.ifft(n.fft.ifftshift(RATIO,axes=1),axis=1)
        print ratio.shape
        #mir.data_array[mir_id,0,:,0] = ratio
            
           # mir.data_array[mir_id,:,:,:] = (res.data_array[res_id,:,:,:]*mir.data_array[mir_id,:,:,:])/(n.conj(mir.data_array[mir_id,:,:,:])*mir.data_array[mir_id,:,:,:])
        #except:
        #    print 'Bad things keeps happening.'
        #    continue
print n.mean(n.abs(ratArray[:,0,95:115,0]))
mir.data_array = n.nan_to_num(mir.data_array)
ratArray = n.nan_to_num(ratArray)
mir.flag_array[mir.data_array==0j] = True
mir.antenna_positions = n.zeros((63,3))
bsl = mir.antnums_to_baseline(41,49)==mir.baseline_array
bidx = mir.flag_array==bsl
pl.imshow(n.abs(ratArray[:,0,95:115,0]),aspect='auto',interpolation='none',vmin=0,vmax=1.0)
pl.colorbar()
pl.show()

mir.write_miriad(opts.dirty.split('HP')[0]+opts.suffix) #replace HP tp uv
    
