import uvdata
import numpy as n
from optparse import OptionParser
import matplotlib.pyplot as mp
import os
import ephem

parser = OptionParser()
parser.add_option("-f", "--file", dest="filename",help="Observation to be converted.", metavar="FILE",action='store')
(options,args) = parser.parse_args()
print options.filename.split('/')[-1:]
#import observation
a = uvdata.miriad.Miriad()
a.read_miriad(options.filename)
#print a.phase_center_ra,a.phase_center_epoch
print a.Npols
a.Npols = 2
print 'Npols set to: ',a.Npols
#set phase
#a.altitude = 377.827
#a.phase(dec=a.zenith_dec[0],ra=a.zenith_ra[0],epoch=2257.21020508)
a.data_array = n.nan_to_num(a.data_array)
print a.phase_center_ra,a.zenith_ra
#a.phase(time=a.time_array[0])
os.chdir("/users/jkerriga/data/jkerriga/FHDOutput")
#save as uvfits
print 'Saving...'+options.filename.split('/')[-1:][0]+'--->'+'P'+options.filename.split('/')[-1:][0]+'.uvfits'
a.write_uvfits('P'+options.filename.split('/')[-1:][0]+'.uvfits',spoof_nonessential=True,force_phase=False)
