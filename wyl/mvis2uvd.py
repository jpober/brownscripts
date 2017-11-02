import sys,optparse,aipy,glob
import numpy as np, mp2cal 
import pyuvdata.uvdata as uvd

o = optparse.OptionParser()
o.set_usage('mvis2uvd.py [options] obsid') #only takes 1 obsid
o.set_description(__doc__)
o.add_option('-d',dest='datpath',default='/users/wl42/data/wl42/FHD_out/fhd_MWA_PhaseII_EoR0/',type='string', help='Path to data. Include final / in path.')
o.add_option('-s',dest='solpath',default='/users/wl42/data/wl42/Nov2016EoR0/mdl_sol/',type='string', help='Path to omnical solutions. Include final / in path.')
o.add_option('-o',dest='outpath',default='/users/wl42/data/wl42/MDLVIS/',type='string', help='Path to save uvfits. Include final / in path.')
opts,args = o.parse_args(sys.argv[1:])
exec('from PhaseII_cal import *')
obsid = args[0]
uv = uvd.UVData()
print '     Loading data'
fhdlist = glob.glob(opts.datpath+'vis_data/'+obsid+'*') + glob.glob(opts.datpath+'metadata/'+obsid+'*')
uv.read_fhd(fhdlist,run_check=False,run_check_acceptability=False)
print '     Loading mdlvis'
npz_x = np.load(opts.solpath+obsid+'.xx.omni.npz')
npz_y = np.load(opts.solpath+obsid+'.yy.omni.npz')
ant = []
for k in npz_x.keys():
    if k[0].isdigit(): ant.append(int(k[0:-1]))
ant.sort()
mdvis = {'xx':{}, 'yy':{}}
info = mp2cal.wyl.pos_to_info(antpos,ants=ant)
a1, a2 = [], []
for ii in range(info.nAntenna):
    for jj in range(ii,info.nAntenna):
        a1.append(ii)
        a2.append(jj)
ant_dict = {}
for a in ant: ant_dict[info.ant_index(a)] = a
reds = info.get_reds()
reds_ind = {}
ubls = []
chix = npz_x['chisq2']
chiy = npz_y['chisq2']
maskx = chix > 1.2
masky = chiy > 1.2
flag = {}
flag['xx'] = np.logical_or(npz_x['flags'], maskx)
flag['yy'] = np.logical_or(npz_y['flags'], masky)
for key in npz_x.keys():
    if key.startswith('<'):
        bl,pol = key.split()
        bl = tuple(map(int,bl[1:-1].split(',')))
        mdvis[pol][bl] = npz_x[key]
        ubls.append(bl)
for key in npz_y.keys():
    if key.startswith('<'):
        bl,pol = key.split()
        bl = tuple(map(int,bl[1:-1].split(',')))
        mdvis[pol][bl] = npz_y[key]
for r in reds:
    ubl = None
    for bl in ubls:
        if bl in r: ubl = bl
    if ubl is None: continue
    for b in r: reds_ind[b] = ubl
Nbls0 = uv.Nbls
Nbls1 = info.nAntenna*(info.nAntenna+1)/2
b0 = 128*uv.ant_1_array[:Nbls0] + uv.ant_2_array[:Nbls0]
uv.Nbls = Nbls1
uv.Nants_data = info.nAntenna
uv.Nants_telescope = info.nAntenna
uv.Nblts = uv.Ntimes*uv.Nbls
times = np.resize(np.unique(uv.time_array),(uv.Nbls,uv.Ntimes)).T
uv.time_array = np.resize(times,(times.size))
lsts = np.resize(np.unique(uv.lst_array),(uv.Nbls,uv.Ntimes)).T
uv.lst_array = np.resize(lsts,(lsts.size))
uvw = np.zeros((uv.Nblts,3))
nsample = np.zeros((uv.Nblts,uv.Nspws,uv.Nfreqs,uv.Npols))
uv.ant_1_array = np.array(a1*uv.Ntimes)
uv.ant_2_array = np.array(a2*uv.Ntimes)
b1 = 128*uv.ant_1_array[:Nbls1] + uv.ant_2_array[:Nbls1]
for ii in range(uv.Nbls):
    i = b1[ii]/128
    j = b1[ii]%128
    ai = ant_dict[i]
    aj = ant_dict[j]
    try:
        ind = np.where(b0 == 128*ai + aj)[0][0]
        uvw[ii::Nbls1] = uv.uvw_array[ind::Nbls0]
    except:
        ind = np.where(b0 == 128*aj + ai)[0][0]
        uvw[ii::Nbls1] = -uv.uvw_array[ind::Nbls0]
    nsample[ii::Nbls1] = uv.nsample_array[ind::Nbls0]
uv.uvw_array = uvw
uv.nsample_array = nsample
uv.data_array = np.zeros((uv.Nblts,uv.Nspws,uv.Nfreqs,uv.Npols),dtype=np.complex64)
uv.flag_array = np.zeros((uv.Nblts,uv.Nspws,uv.Nfreqs,uv.Npols),dtype=bool)
uv.baseline_array = uv.antnums_to_baseline(uv.ant_1_array,uv.ant_2_array)
uv.antenna_positions = np.zeros((info.nAntenna,3))
uv.antenna_numbers = np.arange(info.nAntenna)
uv.antenna_names = []
for ii in range(info.nAntenna): uv.antenna_names.append(str(ii))

for pp in ['xx','yy']:
    pn = aipy.miriad.str2pol[pp]
    pid = np.where(uv.polarization_array==pn)[0][0]
    for ii in range(uv.Nbls):
        i = a1[ii]
        j = a2[ii]
        ai = ant_dict[i]
        aj = ant_dict[j]
        if (ai,aj) in reds_ind.keys(): vis = mdvis[pp][reds_ind[(ai,aj)]]
        elif (aj,ai) in reds_ind.keys(): vis = mdvis[pp][reds_ind[(aj,ai)]].conj()
        else: continue
        uv.data_array[:,0][:,:,pid][ii::uv.Nbls] = vis
        uv.flag_array[:,0][:,:,pid][ii::uv.Nbls] = flag[pp]
outuvfits = opts.outpath + obsid + '_mvis.uvfits'
print '     Writing ' + outuvfits
uv.write_uvfits(outuvfits,spoof_nonessential=True)


