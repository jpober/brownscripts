import pyuvdata
import numpy as n
import optparse,sys
import matplotlib.pyplot as mp
import os
import ephem

parser = optparse.OptionParser()
parser.add_option("-x", "--x", dest="xx",help="xx pol", metavar="FILE",action='store')
parser.add_option("-y", "--y", dest="yy",help="yy pol", metavar="FILE",action='store')
(options,args) = parser.parse_args()

opts,args = parser.parse_args(sys.argv[1:])
#cwd = os.getcwd()

xx = pyuvdata.miriad.Miriad()
xx.read_miriad(opts.xx)

yy = pyuvdata.miriad.Miriad()
yy.read_miriad(opts.yy)

xx.data_array = n.nan_to_num(xx.data_array)
xx.flag_array[xx.data_array==0j] = True

yy.data_array = n.nan_to_num(yy.data_array)
yy.flag_array[yy.data_array==0j] = True

shp = n.shape(xx.data_array)

newData = n.zeros((shp[0],1,shp[2],2)).astype(complex)
newData[:,0,:,0] = xx.data_array[:,0,:,0]
newData[:,0,:,1] = yy.data_array[:,0,:,0]

newFlag = n.zeros((shp[0],1,shp[2],2)).astype(bool)
newFlag[:,0,:,0] = xx.flag_array[:,0,:,0]
newFlag[:,0,:,1] = yy.flag_array[:,0,:,0]

newNsamp = n.zeros((shp[0],1,shp[2],2))
newNsamp[:,0,:,0] = xx.nsample_array[:,0,:,0]
newNsamp[:,0,:,1] = yy.nsample_array[:,0,:,0]


xx.Npols = 2
xx.data_array = newData
xx.flag_array = newFlag
xx.nsample_array = newNsamp
xx.polarization_array = n.array([-5,-8])
    #os.chdir("/users/jkerriga/data/jkerriga/FHDOutput")
xx.phase(dec=xx.zenith_dec[0],ra=xx.zenith_ra[0],epoch=2000.0)
fname = opts.xx.split('.xx')
fname = fname[0]+fname[1]
    #save as uvfits
print 'Saving...'+opts.xx+'--->'+ fname +'.uvfits'
#xx.write_uvfits(fname +'.uvfits',spoof_nonessential=True,force_phase=False)
xx.write_miriad(fname)
