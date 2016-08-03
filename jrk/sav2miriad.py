import uvdata
import numpy as n
import optparse, sys, os

o = optparse.OptionParser()
o.add_option("-s", "--subtract", dest="subtract",help="Return residual observation.", action='store_true') 
opts,args = o.parse_args(sys.argv[1:])
print args

filelist = args
base = args[0].split('_')[0]
dirty = uvdata.fhd.FHD()
dirty.read_fhd(filelist)
dirty.unphase_to_drift()
dirty.antenna_positions = n.zeros((dirty.Nants_data,3))
#dirty.zenith_ra = n.ones(dirty.time_array.shape)*dirty.phase_center_ra
#dirty.zenith_dec =n.ones(dirty.time_array.shape)*dirty.phase_center_dec
dirty.write_miriad(base+'H')
print opts.subtract
if opts.subtract == True:
    model = uvdata.fhd.FHD()
    model.read_fhd(filelist,use_model=True)
    model.unphase_to_drift()
    print 'Applying subtractions to data and outputting to miriad...'
    residual = dirty.data_array - model.data_array
    print 'Dirty: ',n.mean(dirty.data_array),'Residual: ',n.mean(residual)
    dirty.data_array = residual
    dirty.write_miriad(base+'S')
#else:
#    print os.path.dirname(os.path.realpath(__file__))
#    print 'Saving Dirty to miriad...'
#    dirty.write_miriad(base+'H')



