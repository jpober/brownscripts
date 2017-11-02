import uvdata,numpy as n
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-f", dest="filename",help="File to apply mask to", metavar="FILE",action='store')
(options,args) = parser.parse_args()

file = str(options.filename).split('/')[-1:]
print file[0]

pyuv = uvdata.miriad.Miriad()
pyuv.read_miriad(options.filename)
mask = []
for i in range(0,pyuv.Nblts):
    if n.sum(pyuv.data_array[i,:,:,:]) != 0j:
        mask.append(i)

pyuv.Nbls = len(mask)/pyuv.Ntimes
pyuv.Nblts = len(mask)
#pyuv.Npols = 2
pyuv.ant_1_array = pyuv.ant_1_array[mask]
pyuv.ant_2_array = pyuv.ant_2_array[mask]

antnames = []

for i in range(1,65):
    antnames.append('ANT'+str(i))

pyuv.antenna_names = antnames
pyuv.baseline_array = pyuv.baseline_array[mask]
pyuv.data_array = pyuv.data_array[mask,:,:,0:2]
pyuv.flag_array = pyuv.flag_array[mask,:,:,0:2]
#pyuv.flag_array = n.zeros(pyuv.flag_array.shape,dtype='Bool')
pyuv.lst_array = pyuv.lst_array[mask]
pyuv.nsample_array = pyuv.nsample_array[mask,:,:,0:2]
pyuv.polarization_array = pyuv.polarization_array[0:2]
pyuv.time_array = pyuv.time_array[mask]
pyuv.uvw_array = pyuv.uvw_array[:,mask]
pyuv.data_array = n.nan_to_num(pyuv.data_array)
print 'Writing out uvfits...'
pyuv.write_uvfits('S'+file[0]+'.uvfits',spoof_nonessential=True,force_phase=True)
