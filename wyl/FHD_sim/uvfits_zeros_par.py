#!/bin/env python

from astropy.io import fits
from astropy.time import Time
import aipy as a
import sys, optparse, csv
import numpy as n
import ephem
from uvdata.uv import UVData, UVParameter


d1 = '0,2,4,6,0,2,4,6,0,2,4,6,0,2,4,6'
d2 = '0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3'
d3 = '0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0'
d4 = '3,2,1,0,3,2,1,0,3,2,1,0,3,2,1,0'
d5 = '6,4,2,0,6,4,2,0,6,4,2,0,6,4,2,0'
d6 = '9,6,3,0,9,6,3,0,9,6,3,0,9,6,3,0'
delaylist = [d1]*16+[d2]*15+[d3]*15+[d4]*15+[d5]*15+[d6]*18
RAlist = [355.7747416798699,
 356.2764481557374,
 356.7781830979642,
 357.3133366733206,
 357.8150733065082,
 358.3167830879389,
 358.8185209542662,
 359.3536780792277,
 359.8554085924123,
 0.3571325760772602,
 0.8588570209470606,
 1.3940298708495,
 1.895754640387602,
 2.397479889200486,
 2.899204944525079,
 3.434379000398047,
 356.2622658470658,
 356.7639825860167,
 357.2657277679912,
 357.8008922468393,
 358.3026387849658,
 358.8043590054173,
 359.3061070086476,
 359.8412781183658,
 0.3430114008402328,
 0.8447447894349475,
 1.346478500841574,
 1.881661637666322,
 2.383395829809858,
 2.885130148169588,
 3.386864667899797,
 356.3038382582308,
 356.8055857896703,
 357.3073065338226,
 357.80905542086,
 358.3442243890356,
 358.8459748440076,
 359.34769868155,
 359.8494338190942,
 0.3846190644417725,
 0.8863557388661116,
 1.388092505092876,
 1.88982950720896,
 2.425015931293939,
 2.926753343844366,
 3.428491012285335,
 356.3120171014037,
 356.8471900834547,
 357.348929291076,
 357.8506697826101,
 358.3523832181399,
 358.8875873394217,
 359.3893022997333,
 359.8910411529602,
 0.3927744324354322,
 0.9279567460803007,
 1.429690487899068,
 1.931424621328283,
 2.433158800085693,
 2.968342116862803,
 3.470076475698106,
 356.2979731701596,
 356.7996748960408,
 357.3348659453688,
 357.8365692898507,
 358.3383008024638,
 358.8400328929332,
 359.3751998961305,
 359.8769181779305,
 0.3786421132368082,
 0.8803663258642238,
 1.415539214586308,
 1.917264077161869,
 2.418988933198115,
 2.92071413861356,
 3.45588779305723,
 356.1598432935251,
 356.66155591122,
 357.1632420147617,
 357.6984171544369,
 358.2001048831038,
 358.7018206908715,
 359.2035376031143,
 359.7386881100908,
 0.2403946022632056,
 0.7421034552278809,
 1.243812339405027,
 1.778969056294653,
 2.280678446367146,
 2.782388133879986,
 3.284098155471878,
 3.819255570503198,
 4.320965585445152,
 4.822675828634868]
DEClist = [-25.96440923165652,
 -25.96445713758903,
 -25.96450043872224,
 -25.96453832505396,
 -25.96456710382853,
 -25.96459150336854,
 -25.96460969720478,
 -25.96462225526444,
 -25.9646276175542,
 -25.96462677005856,
 -25.96461959363735,
 -25.96460533564646,
 -25.96458715420021,
 -25.96456093921331,
 -25.96453034571866,
 -25.96449075695374,
 -26.58037965217246,
 -26.58042234560487,
 -26.58046042697623,
 -26.58049273884796,
 -26.58051640402648,
 -26.58053544645961,
 -26.58054668921665,
 -26.58055376829496,
 -26.58055375844145,
 -26.58054765402211,
 -26.58053533621952,
 -26.58051534437507,
 -26.58049189625293,
 -26.58046041229938,
 -26.5804245478187,
 -26.78347202115501,
 -26.78351435268676,
 -26.78355059766967,
 -26.78358234490085,
 -26.78360765045253,
 -26.78362485500215,
 -26.78363755342135,
 -26.78364391874284,
 -26.78364239177672,
 -26.78363591308383,
 -26.78362333899589,
 -26.78360455119787,
 -26.78357936577551,
 -26.78354774106804,
 -26.78351161569386,
 -26.58038369893068,
 -26.58042860240856,
 -26.58046407755443,
 -26.58049493536289,
 -26.58051799556159,
 -26.58053585791104,
 -26.5805477900263,
 -26.58055350965854,
 -26.58055289673568,
 -26.58054563531007,
 -26.58053230088981,
 -26.58051263528075,
 -26.58048687806837,
 -26.58045438153886,
 -26.58041567633605,
 -25.96445774613692,
 -25.96449909364679,
 -25.96453843022946,
 -25.96456708359223,
 -25.96459112252485,
 -25.9646090742525,
 -25.964619786562,
 -25.96462478875087,
 -25.96462370095679,
 -25.9646181114369,
 -25.9646034792171,
 -25.9645833517346,
 -25.96455872538269,
 -25.96452777555406,
 -25.96448634585424,
 -24.91387505207776,
 -24.91391685215951,
 -24.91395416672596,
 -24.91398713871044,
 -24.91401164950478,
 -24.91402996251302,
 -24.91404195574739,
 -24.91404986125075,
 -24.91404915642058,
 -24.91404383778373,
 -24.91403072905496,
 -24.91401161393318,
 -24.91398739910824,
 -24.91395698673738,
 -24.9139202602099,
 -24.91387438688075,
 -24.91382669291075,
 -24.91377098766836]
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
exec('from %s import sim_prms as _sim_params' % args[0].split('.')[0] )
exec('from %s import instr_prms as _instr_params' % args[0].split('.')[0] )
	
uvd = UVData()

default_attrs=[aa for aa in dir(uvd) if not aa.startswith('__')]
#default_attrs.append("delays")

for key, val in _sim_params.iteritems():
    if key in default_attrs:
	uvp = getattr(uvd, key)
    else:
        uvp = UVParameter(key)
    uvp = val
    setattr(uvd,key,uvp)
for key, val in _instr_params.iteritems():
    if key in default_attrs:
	uvp = getattr(uvd, key)
    else:
        uvp = UVParameter(key)
    uvp = val
    setattr(uvd,key,uvp)

print 'check:', uvd.channel_width, uvd.sfreq
# Generate an antenna array from a cal file.
aa = a.cal.get_aa(args[1].split('.')[0], uvd.channel_width/1e9, uvd.sfreq/1e9, uvd.Nfreqs)  #Saved in Hz, AIPY expects GHz
Nants = len(aa)

uvd.Nants_telescope = Nants
uvd.Nants_data = Nants
uvd.latitude, uvd.longitude, uvd.altitude = aa.lat, aa.long, aa.elev    ### Altitude vs elevation?

uvd.channel_width = uvd.channel_width

polar = uvd.polarization_array

if polar == 'stokes':
    uvd.polarization_array=n.arange(1,5,1)
elif polar == 'circ':
    uvd.polarization_array=n.arange(-1,-5,-1)
elif polar == 'lin':
    uvd.polarization_array=n.arange(-5,-9,-1)


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
	uvd.extra_keywords['delays'] = d3

t0 = tzero+en*uvd.Ntimes*dt
#Data array
uvd.data_array = n.zeros((nbl * uvd.Ntimes, uvd.Nspws, uvd.Nfreqs, uvd.Npols), dtype=n.complex)

#Time array:
## Format -- Repeat times for each baseline number.

tims = n.arange(uvd.Ntimes, dtype=n.float) * dt + t0
uvd.time_array = n.sort(n.tile(tims, nbl))    #Should be of length Nblts, baseline fast time slow
uvd.Nblts = len(uvd.time_array)    #nbls * Ntimes

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
uvd.phase_center_ra= ephem.hours(RA)
uvd.phase_center_dec  = ephem.degrees(dec)
uvd.object_name= "zenith"
uvd.phase_center_epoch = ephem.J2000
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
#uvd.lst_array.value = n.zeros(uvd.Ntimes.value*nbl)
if nout > 1:
     ofi = ofile.split(".")[0] +"_"  +  str(en) + ".uvfits"
print ofi
uvd.write_uvfits(ofi, spoof_nonessential=True)

#hdu = fits.open(ofi,mode='update')
#hdu[0].header['DELAYS'] = d3 # delaylist[en-1]
#hdu[0].header['RA'] = RAlist[en-1]
#hdu[0].header['DEC'] = DEClist[en-1]
#hdu[0].header['RAPHASE'] = 0.0
#hdu[0].header['DECPHASE'] = -27.0
#hdu.flush()
#hdu.close()
