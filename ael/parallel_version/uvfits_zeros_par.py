#!/bin/env python

from astropy.io import fits
from astropy.time import Time
import aipy as a
import sys, optparse, csv
import numpy as n
import ephem
from uvdata.uv import UVData, UVProperty


## NB -- pyuvdata and uvfits use Hz and Seconds for everything. AIPY assumes GHz and ns.

#Parallel version of uvfits_zeros_gen


o = optparse.OptionParser()
o.set_usage('uvfits_zeros_gen [-o <outfile>] <cal file of coords> <sim_settings_file>')

o.add_option('-o', dest='ofile', help='Destination filename')
o.add_option('-i', dest='en', help='Which job?')
o.add_option('-N', dest='nout', help='How many?')

opts,args = o.parse_args(sys.argv[1:])

if not opts.ofile == None:
    ofile=opts.ofile
else:
    ofile='empty.uvfits'
if not opts.en == None:
    en=int(opts.en)
else:
    en=1
if not opts.nout == None:
    nout=opts.nout
else:
    nout=1

#Import the sim and instrument settings from the file.
exec('from %s import sim_prms as _sim_params' % args[1].split('.')[0] )
exec('from %s import instr_prms as _instr_params' % args[1].split('.')[0] )
	
uvd = UVData()

default_attrs=[aa for aa in dir(uvd) if not aa.startswith('__')]

for key, val in _sim_params.iteritems():
    if key in default_attrs:
	uvp = getattr(uvd, key)
    else:
        uvp = UVProperty()
    uvp.value= val
    setattr(uvd,key,uvp)
for key, val in _instr_params.iteritems():
    if key in default_attrs:
	uvp = getattr(uvd, key)
    else:
        uvp = UVProperty()
    uvp.value= val
    setattr(uvd,key,uvp)


# Generate an antenna array from a cal file.
aa = a.cal.get_aa(args[0].split('.')[0], uvd.channel_width.value/1e9, uvd.sfreq.value/1e9, uvd.Nfreqs.value)  #Saved in Hz, AIPY expects GHz
Nants = len(aa)

uvd.Nants_telescope.value = Nants
uvd.Nants_data.value = Nants
uvd.latitude.value, uvd.longitude.value, uvd.altitude.value= aa.lat, aa.long, aa.elev    ### Altitude vs elevation?

uvd.channel_width.value = uvd.channel_width.value

polar = uvd.polarization_array.value

if polar == 'stokes':
    setattr(uvd.polarization_array, 'value', n.arange(1,5,1))
elif polar == 'circ':
    setattr(uvd.polarization_array, 'value', n.arange(-1,-5,-1))
elif polar == 'lin':
    setattr(uvd.polarization_array, 'value', n.arange(-5,-9,-1))


#Time -- The rest of this uses julian date. Take param from sim_prms (in uvd) and 
# convert to ephem object. Then use that to get julian date.

if uvd.time_format.value == 'julian':
    tzero = ephem.date(uvd.start_time.value - 2415020.)
    tfin = ephem.date(uvd.end_time.value - 2415020.)

tzero = ephem.julian_date(tzero)
tfin = ephem.julian_date(tfin)

uvd.dateobs.value = tzero
dt = (tfin - tzero)/ float(uvd.Ntimes.value * int(nout))

try:
   uvd.x_telescope.value, uvd.y_telescope.value, uvd.z_telescope.value = aa.get_xyz_telescope()
except AttributeError:
   pass

#Frequency array
uvd.freq_array.value = n.array([aa.get_afreqs()* 1e9])

#Baseline array:
ant_exclude=[]
if hasattr(uvd,'bad_antenna'):
	ant_exclude = uvd.bad_antenna          #List antennas to be excluded here, for whatever reason.
bls = n.array([uvd.antnums_to_baseline(j,i,attempt256=True) 
       for i in range(0,len(aa)) 
       for j in range(i,len(aa)) if not j in ant_exclude and not i in ant_exclude ])

uvd.baseline_array.value = n.tile(bls, uvd.Ntimes.value)

#number of baselines
nbl = len(bls)
uvd.Nbls.value = nbl

#Antennas
uvd.antenna_indices.value = n.arange(1,Nants+1,1)       #1 indexed, not 0
uvd.antenna_names.value = ["ANT"+str(i) for i in uvd.antenna_indices.value]
uvd.antenna_positions.value = n.array([ant.pos for ant in aa])
uvd.ant_1_array.value, uvd.ant_2_array.value = \
          uvd.baseline_to_antnums(uvd.baseline_array.value)

#Delays
uvd.delays = UVProperty()
uvd.delays.value=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
uvd.delays.form = (6,)

t0 = tzero+en*uvd.Ntimes.value*dt
#Data array
uvd.data_array.value = n.zeros((nbl * uvd.Ntimes.value, uvd.Nspws.value, uvd.Nfreqs.value,uvd.Npols.value), dtype=n.complex)

#Time array:
## Format -- Repeat times for each baseline number.

tims = n.arange(uvd.Ntimes.value, dtype=n.float) * dt + t0
uvd.time_array.value = n.sort(n.tile(tims, nbl))    #Should be of length Nblts, baseline fast time slow
uvd.Nblts.value = len(uvd.time_array.value)    #nbls * Ntimes

t = t0 + dt
#t0 = max(tims)
aa.set_jultime(t)
RA = str(aa.sidereal_time())
dec = str(aa.lat)
src = RA+'_'+dec

#Artificial point
#src="23:43:06.0_-25:57:51.84"
RA,dec = src.split('_')
print src
#sys.exit()


#uvd.phase_center_ra.value= n.rad2deg(float(repr(ephem.hours(RA))))
#uvd.phase_center_dec.value  = n.rad2deg(float(repr(ephem.degrees(dec))))
uvd.phase_center_ra.value= ephem.hours(RA)
uvd.phase_center_dec.value  = ephem.degrees(dec)
uvd.object_name.value= "zenith"
uvd.phase_center_epoch.value = ephem.J2000
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

uvd.uvw_array.value = n.array(uvw_array).T * a.const.len_ns / 100.


del uvw_array

#TODO --- Check line below. uvd needs a flag array. True = flagged
uvd.flag_array.value = n.zeros(shape=uvd.data_array.value.shape, dtype=n.bool)
uvd.nsample_array.value = n.ones(shape=uvd.data_array.value.shape, dtype=n.intc)

extra_attrs=[atr for atr in dir(uvd) if not atr.startswith('__') if not atr in default_attrs]
for attr in extra_attrs:
		delattr(uvd, attr)

#uvd.instrument.expected_type=str
#delattr(uvd, 'end_time')
uvd.set_lsts_from_time_array()
#uvd.lst_array.value = n.zeros(uvd.Ntimes.value*nbl)
if nout > 1:
     ofi = ofile.split(".")[0] +"_"  +  str(en) + ".uvfits"
print ofi
uvd.write_uvfits(ofi, spoof_nonessential=True)
