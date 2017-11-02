import uvdata
import numpy as n
import optparse, sys, os
from scipy.io import readsav
from glob import glob

o = optparse.OptionParser()
#o.add_option("-s", "--subtract", dest="subtract",help="Return residual observation.", action='store_true') 
opts,args = o.parse_args(sys.argv[1:])
print args
#caloc=glob('../calibration/*sav')[0]
#calibrations = readsav(caloc)
#gains = calibrations['cal'][0]['GAIN']
filelist = args
base = args[0].split('_')[0]
dirty = uvdata.fhd.FHD()
dirty.read_fhd(filelist)
#dirty.unphase_to_drift()
dirty.antenna_positions = n.zeros((dirty.Nants_data,3))
#model = uvdata.fhd.FHD()
#model.read_fhd(filelist,use_model=True)
#model.unphase_to_drift()


# Rescale all gains back to original calibration gains
#for p in range(2):
#    for i in range(0,64):
#        for j in range(0,64):
#            bsl=dirty.antnums_to_baseline(i,j)
#            bslMsk = bsl==dirty.baseline_array
#            gainFactor = gains[p][i,:]*n.conj(gains[p][j,:])
#           print gainFactor,bslMsk.sum()
#            dirty.data_array[bslMsk,0,:,p] = dirty.data_array[bslMsk,0,:,p]*gainFactor/2.
#            model.data_array[bslMsk,0,:,p] = model.data_array[bslMsk,0,:,p]*gainFactor/2.

dirty.unphase_to_drift()
#model.unphase_to_drift()
dirty.data_array = n.conj(dirty.data_array)
#model.data_array = n.conj(model.data_array)

# Remove data outside of channels 42(24) and 165, these cause huge issues for the
# Wideband filter
dirty.flag_array[dirty.data_array==0j]=True
#Try and make stokes I...I guess
dirty.flag_array = dirty.flag_array[:,:,:,0]
dirty.flag_array = dirty.flag_array.reshape(n.append(dirty.flag_array.shape,1))
dirty.data_array = dirty.data_array[:,:,:,0]+dirty.data_array[:,:,:,1]
dirty.data_array = dirty.data_array.reshape(n.append(dirty.data_array.shape,1))
dirty.nsample_array = dirty.nsample_array[:,:,:,0]+dirty.nsample_array[:,:,:,1]
dirty.nsample_array = dirty.nsample_array.reshape(n.append(dirty.nsample_array.shape,1))
print dirty.freq_array.shape
dirty.freq_array = dirty.freq_array[:,:]
print dirty.freq_array.shape
dirty.polarization_array = n.array([1])
dirty.Npols=1


#model.data_array = model.data_array[:,:,:,0]+model.data_array[:,:,:,1]
#model.data_array = model.data_array.reshape(n.append(model.data_array.shape,1))
#dirty.freq_array=dirty.freq_array[:,24:165] #channels were 24:165
#dirty.flag_array = dirty.flag_array[:,:,24:165,:]
#dirty.data_array = dirty.data_array[:,:,24:165,:]
#dirty.nsample_array = dirty.nsample_array[:,:,24:165,:]
dirty.Nfreqs = len(dirty.data_array[0,0,:,0])

dirty.write_miriad(base+'HP')

print 'Applying subtractions to data and outputting to miriad...'
#model.data_array=model.data_array[:,:,24:165,:]

#residual = dirty.data_array - model.data_array
#residual[:,0,:,0] = dirty.data_array[:,0,:,0] - model.data_array[:,0,:,0]
#residual[:,0,:,1] = dirty.data_array[:,0,:,1]# - model.data_array[:,0,:,0]
#print 'Dirty: ',n.sum(n.abs(dirty.data_array)),'Residual: ',n.sum(n.abs(residual))
#dirty.data_array = residual
#dirty.write_miriad(base+'SP')




