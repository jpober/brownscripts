import numpy as n
import uvdata
from glob import glob
import pylab as pl
import optparse,sys

o = optparse.OptionParser()
opts,args = o.parse_args(sys.argv[1:])
#x = glob('/users/jkerriga/data/jkerriga/2DayOutput/fhd_TestCalFHDPRSV/vis_data/*')
#y = glob('/users/jkerriga/data/jkerriga/2DayOutput/fhd_TestCalFHDPRSV/metadata/*')
#z = x+y
#file2=z
#file1='/users/jkerriga/data/jkerriga/2DayOutput/fhd_CalSrc100Poly0/vis_data/Pzen.2456242.30605.uvcRREcACOTUcHPB'
#file2='/users/jkerriga/data/jkerriga/2DayOutput/fhd_CalSrc100Poly0/vis_data/Pzen.2456242.30605.uvcRREcACOTUcSPB'
file1=args[0]
file2=args[1]
#sav = uvdata.fhd.FHD()
#sav.read_fhd(z)
try:
    sav=uvdata.miriad.Miriad()#dirty
    sav.read_miriad(file1,run_sanity_check=False)
except:
    pass
#sav.unphase_to_drift()
try:
    mir = uvdata.miriad.Miriad()#residual
    mir.read_miriad(file2,run_sanity_check=False)
except:
    pass
#mir.unphase_to_drift()
#sav.unphase_to_drift()

bsl=mir.antnums_to_baseline(6,33)

bslSav = sav.antnums_to_baseline(6,33)
print bsl,bslSav
savmsk = sav.baseline_array==bslSav
savMasked = sav.data_array[savmsk,0,:,0]
#savMasked = n.abs(n.nan_to_num(savMasked))

savMaskedy = sav.data_array[savmsk,0,:,0]
#savMaskedy = n.abs(n.nan_to_num(savMaskedy))

#print n.max(n.abs(savMasked))
mirmsk = mir.baseline_array==bsl

mirMasked = mir.data_array[mirmsk,0,:,0]
#mirMasked = n.abs(n.nan_to_num(mirMasked))
print 'Dirty: ',n.mean(n.abs(sav.data_array[:,0,:,0])),n.sum(n.abs(sav.data_array)),' Residual: ',n.mean(n.abs(mir.data_array[:,0,:,0])),n.sum(n.abs(mir.data_array))
#mirMaskedy = mir.data_array[mirmsk,0,:,0]
#mirMaskedy = n.abs(n.nan_to_num(mirMaskedy))
#print n.max(n.abs(mirMasked))
#mirMasked=n.nan_to_num(mirMasked)
#savMasked=n.nan_to_num(savMasked)
div = n.abs(savMasked)/n.abs(mirMasked)
#div = n.abs(n.nan_to_num(div))
div[div>10**1]=0
print 'Mean of ratio: ',n.mean(div)
#divy = savMaskedy/mirMaskedy
#divy = n.abs(n.nan_to_num(divy))
print 'Mean: ',n.mean(div),'Max: ',n.max(div),'Min: ',n.min(div),'Sum: ',n.sum(div)
print div.shape
pl.subplot(311)
pl.imshow(n.log10(n.abs(savMasked)),aspect='auto')
pl.colorbar()
pl.subplot(312)
pl.imshow(n.log10(n.abs(mirMasked)),aspect='auto')
pl.colorbar()
pl.subplot(313)
pl.imshow((div),aspect='auto')
pl.colorbar()
pl.show()
