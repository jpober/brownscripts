#!/bin/env python

from astropy.io import fits
from astropy.time import Time
import aipy as a
import sys, optparse, csv
import numpy as n
import ephem
from uvdata.uv import UVData, UVProperty


## NB -- pyuvdata and uvfits use Hz and Seconds for everything. AIPY assumes GHz and ns.


o = optparse.OptionParser()
o.set_usage('uvfits_zeros_gen [-o <outfile>] <cal file of coords> <sim_settings_file>')

o.add_option('-o', dest='ofile', help='Destination filename')

opts,args = o.parse_args(sys.argv[1:])

if not opts.ofile == None:
    ofile=opts.ofile
else:
    ofile='empty.uvfits'

#uv.uv_property

#Import the sim and instrument settings from the file.
exec('from %s import sim_prms as _sim_params' % args[1].split('.')[0] )
exec('from %s import instr_prms as _instr_params' % args[1].split('.')[0] )

uvd = UVData()

for key, val in _sim_params.iteritems():
    uvp = UVProperty()
    uvp.value = val
    setattr(uvd,key,uvp)
for key, val in _instr_params.iteritems():
    uvp = UVProperty()
    uvp.value = val
    setattr(uvd,key,uvp)
	
# Generate an antenna array from a cal file.
aa = a.cal.get_aa(args[0].split('.')[0], uvd.channel_width.value/1e9, uvd.sfreq.value/1e9, uvd.Nfreqs.value)  #Saved in Hz, AIPY expects GHz
Nants = len(aa)
uvd.Nants.value = Nants
uvd.latitude.value, uvd.longitude.value, uvd.altitude.value= aa.lat, aa.long, aa.elev ### Altitude vs elevation?

polar = uvd.polarization_array.value
if polar == 'stokes':
    setattr(uvd.polarization_array, 'value', n.arange(1,5,1.).tolist())
elif polar == 'circ':
    setattr(uvd.polarization_array, 'value', n.arange(-1,-5,-1.).tolist())
elif polar == 'lin':
    setattr(uvd.polarization_array, 'value', n.arange(-5,-9,-1.).tolist())


#Time -- The rest of this uses julian date. Take param from sim_prms (in uvd) and 
# convert to ephem object. Then use that to get julian date.

if uvd.time_format.value == 'julian':
    tzero = ephem.date(uvd.start_time.value - 2415020.)
    tfin = ephem.date(uvd.end_time.value - 2415020.)

tzero = ephem.julian_date(tzero)
tfin = ephem.julian_date(tfin)

uvd.dateobs.value = tzero
dt = (tfin - tzero)/ float(uvd.Ntimes.value)

try:
   uvd.x_telescope.value, uvd.y_telescope.value, uvd.z_telescope.value = aa.get_xyz_telescope()
except AttributeError:
   pass

#Frequency array
uvd.freq_array.value = n.array([aa.get_afreqs()* 1e9])

#Baseline array:
## Format -- Repeat baseline numbers for each time.
bls = n.array(aa.bl_indices(auto=False))
uvd.baseline_array.value = n.tile(bls, uvd.Ntimes.value)

#number of baselines
nbl = len(bls)

uvd.Nblts.value = nbl*uvd.Ntimes.value

#Time array:
## Format -- Repeat times for each baseline number.
tims = n.arange(uvd.Ntimes.value, dtype=n.float) * dt + tzero
uvd.time_array.value = n.tile(tims, nbl)

uvd.Nbls.value = nbl

#Data array
uvd.data_array.value = n.zeros((nbl * uvd.Ntimes.value, uvd.Nspws.value, uvd.Nfreqs.value,uvd.Npols.value), dtype=n.intc)

#Antennas
uvd.antenna_indices.value = n.arange(Nants)
uvd.antenna_names.value = ["ANT"+str(i) for i in uvd.antenna_indices.value]
uvd.antenna_positions.value = n.array([ant.pos for ant in aa])
uvd.ant_1_array.value = uvd.antenna_indices.value
uvd.ant_2_array.value = uvd.antenna_indices.value


## NOTE -- AIPY's "get_uvw" converts the baseline to units of 1/lambda, and thus requires dotting with the frequency array. 
## uvfits prefers a set of uvw coords for each time, since it's tracking a source.

#uvw_array = n.array([ aa.get_baseline(i,j,src='z') for i,j in [aa.bl2ij(bl) for bl in bls] ]).T
#uvd.uvw_array = n.tile(uvw_array, Nt) /(a.const.len_ns * 100.)  #in meters

### Phase to zenith at the halfway point:
#t = (tfin - tzero)/2.
t = tzero
aa.set_jultime(t)
RA = str(aa.sidereal_time())
dec = str(aa.lat)
src = RA+'_'+dec

#Phase to 0,-27
#src="0:00:00.0_-27:00:00.00"
#src="355:46:28.8_-25:57:51.84"

print src
uvd.phase_center_ra.value= RA
uvd.phase_center_dec.value  = dec
uvd.object_name.value= "zenith"
uvd.phase_center_epoch.value = 2000
uvd.history.value = ''

srclist,cutoff,catalogs = a.scripting.parse_srcs(src, 'helm')
src = a.cal.get_catalog(args[0], srclist, cutoff, catalogs).values()[0]
#src._epoch=36525.
uvw_array = []
curtime=-1
for t in tims:
	if curtime != t:
		curtime = t
		aa.set_jultime(t)
	src.compute(aa)
	for bl in bls:
		i,j = aa.bl2ij(bl)
		uvw_array.append(aa.get_baseline(i,j,src=src))
#		uvw_array.append(aa.get_baseline(i,j,src='z'))
#		if (i,j) == (88,95):
#			print aa.get_baseline(i,j,src=src)/a.const.len_ns * 100.
#	uvw_array.append([ aa.get_baseline(i,j,src=src) for i,j in [aa.bl2ij(bl) for bl in bls] ])

uvd.uvw_array.value = n.array(uvw_array).T * a.const.len_ns / 100.


del uvw_array

#TODO --- Check line below. uvd needs a flag array. True = flagged
uvd.flag_array.value = n.zeros(shape=uvd.data_array.value.shape, dtype=n.intc)
uvd.nsample_array.value = n.ones(shape=uvd.data_array.value.shape, dtype=n.intc)

#for atr in [a for a in dir(uvd) if not a.startswith('__')]:
#    print atr, type(atr)

#sys.exit(0)


uvd.write_uvfits(ofile)
