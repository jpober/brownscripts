import numpy as np
import pyuvdata.uvdata as uvd
import aipy as a, ephem
import sys, optparse

def decdeg2dms(dd):
    if dd < 0: dp = -dd
    else: dp = dd
    mnt,sec=divmod(dp*3600,60)
    deg,mnt=divmod(mnt,60)
    if dd<0: return str(-int(deg))+':'+str(int(mnt))+':%.2f'%(sec)
    else: return str(int(deg))+':'+str(int(mnt))+':%.2f'%(sec)


o = optparse.OptionParser()

o.set_usage('python uvfits_sim.py [options] obsid')
o.add_option('-o', dest='outdir', default='/users/wl42/data/wl42/perfectpos/',help='output directory')
o.add_option('-i', dest='indir', default='/users/wl42/data/wl42/NOISE/',help='directory to data uvfits')

opts,args = o.parse_args(sys.argv[1:])
obsid = args[0]
print 'loading...'
uv = uvd.UVData()
uv.read_uvfits(opts.indir+obsid+'.uvfits',run_check=False,run_check_acceptability=False)

exec('from PhaseII_cal import *')
aa = a.cal.get_aa('MWA_PhaseII_cal', uv.channel_width/1e9, uv.freq_array[0]/1e9, uv.Nfreqs)
RA = ephem.hours(uv.phase_center_ra)
DEC = decdeg2dms(uv.phase_center_dec_degrees)
src = str(RA)+'_'+DEC
RA,dec = src.split('_')
print src
srclist,cutoff,catalogs = a.scripting.parse_srcs(src, 'helm')
src = a.cal.get_catalog('MWA_PhaseII_cal', srclist, cutoff, catalogs).values()[0]
uvw_array = []
for ii in range(uv.Nblts):
    aa.set_jultime(uv.time_array[ii])
    src.compute(aa)
    i = uv.ant_1_array[ii]
    j = uv.ant_2_array[ii]
    uvw_array.append(aa.get_baseline(i,j,src='z'))
uv.uvw_array = -np.array(uvw_array)*a.const.len_ns/100
print 'writing...'
uv.write_uvfits(opts.outdir+obsid+'_sim.uvfits')
