#!/bin/env python

from astropy.io import fits
from astropy.time import Time
import aipy as a
import sys, optparse, csv
import numpy as n
import ephem
from uvdata.uvbase import UVBase
import uvdata.parameter as uvp
import uvdata.utils as utils
from uvdata.uv import UVData

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
#default_attrs.append("delays")

for key, val in _sim_params.iteritems():
    if key in default_attrs:
	param = getattr(uvd, key)
    else:
        param = uvp.UVParameter(key, description=key)
    param = val
    setattr(uvd,key,param)

for key, val in _instr_params.iteritems():
    if key in default_attrs:
	param = getattr(uvd, key)
    else:
        param = uvp.UVParameter(key, description=key)
    param = val
    setattr(uvd,key,param)


# Generate an antenna array from a cal file.
aa = a.cal.get_aa(args[0].split('.')[0], uvd.channel_width/1e9, uvd.sfreq/1e9, uvd.Nfreqs)  #Saved in Hz, AIPY expects GHz
Nants = len(aa)

#Set telescope parameters
uvd.Nants_telescope = Nants
uvd.Nants_data = Nants
uvd.set_telescope_params()
uvd.antenna_numbers = n.arange(Nants)
#uvd.latitude, uvd.longitude, uvd.altitude = aa.lat, aa.long, aa.elev    ### Altitude vs elevation?

#uvd.channel_width = uvd.channel_width

uvd.is_phased = True

polar = uvd.polarization_array

if polar == 'stokes':
    uvd.polarization_array = n.arange(1,5,1)
elif polar == 'circ':
    uvd.polarization_array = n.arange(-1,-5,-1)
elif polar == 'lin':
    uvd.polarization_array =  n.arange(-5,-9,-1)


#Time -- The rest of this uses julian date. Take param from sim_prms (in uvd) and 
# convert to ephem object. Then use that to get julian date.

if uvd.time_format == 'julian':
    tzero = ephem.date(uvd.start_time - 2415020.)
    tfin = ephem.date(uvd.end_time - 2415020.)

tzero = ephem.julian_date(tzero)
tfin = ephem.julian_date(tfin)

uvd.dateobs = tzero
dt = (tfin - tzero)/ float(uvd.Ntimes * int(nout))

try:
   uvd.x_telescope, uvd.y_telescope, uvd.z_telescope = aa.get_xyz_telescope()
except AttributeError:
   pass

#Frequency array
uvd.freq_array = n.array([aa.get_afreqs()* 1e9])

#Baseline array:
ant_exclude=[]
if hasattr(uvd,'bad_antenna'):
	ant_exclude = uvd.bad_antenna          #List antennas to be excluded here, for whatever reason.
bls = n.array([uvd.antnums_to_baseline(j,i,attempt256=True) 
       for i in range(0,len(aa)) 
       for j in range(i,len(aa)) if not j in ant_exclude and not i in ant_exclude ])

uvd.baseline_array = n.tile(bls, uvd.Ntimes)

#number of baselines
nbl = len(bls)
uvd.Nbls = nbl

#Antennas
uvd.antenna_indices = n.arange(1,Nants+1,1)       #1 indexed, not 0
uvd.antenna_names = ["ANT"+str(i) for i in uvd.antenna_indices]
uvd.antenna_positions = n.array([ant.pos for ant in aa])
uvd.ant_1_array, uvd.ant_2_array = \
          uvd.baseline_to_antnums(uvd.baseline_array)

#Delays
if uvd.instrument == 'MWA':
	uvd.extra_keywords['delays'] = '0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0'

t0 = tzero+en*uvd.Ntimes*dt
#Data array
uvd.data_array = n.zeros((nbl * uvd.Ntimes, uvd.Nspws, uvd.Nfreqs,uvd.Npols), dtype=n.complex)

#Time array:
## Format -- Repeat times for each baseline number.

tims = n.arange(uvd.Ntimes, dtype=n.float) * dt + t0
uvd.time_array = n.sort(n.tile(tims, nbl))    #Should be of length Nblts, baseline fast time slow
uvd.Nblts = len(uvd.time_array)    #nbls * Ntimes

print uvd.Ntimes
print uvd.Nblts
sys.exit()

t = t0 # + dt
#t0 = max(tims)

###
#Temporary
t = tzero
###

aa.set_jultime(t)
RA = str(aa.sidereal_time())
dec = str(aa.lat)
src = RA+'_'+dec

#Artificial point
#src="23:43:06.0_-25:57:51.84"
#src="23:31:31.24_-30:43:17.5"
RA,dec = src.split('_')
print src

uvd.phase_center_ra= ephem.hours(RA)
uvd.phase_center_dec  = ephem.degrees(dec)
uvd.object_name= "zenith"
uvd.phase_center_epoch = 36525.0
uvd.history = ''

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

uvd.uvw_array = n.array(uvw_array).T * a.const.len_ns / 100.

del uvw_array

#TODO --- Check line below. uvd needs a flag array. True = flagged
uvd.flag_array = n.zeros(shape=uvd.data_array.shape, dtype=n.bool)
uvd.nsample_array = n.ones(shape=uvd.data_array.shape, dtype=n.intc)

extra_attrs=[atr for atr in dir(uvd) if not atr.startswith('__') if not atr in default_attrs]

for attr in extra_attrs:
		delattr(uvd, attr)

#uvd.instrument.expected_type=str
#delattr(uvd, 'end_time')
uvd.set_lsts_from_time_array()
if nout > 1:
     ofi = ofile.split(".")[0] +"_"  +  str(en) + ".uvfits"
print ofi
uvd.write_uvfits(ofi, spoof_nonessential=True)
