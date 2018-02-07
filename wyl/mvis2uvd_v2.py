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
obsid = args[0]
uv = uvd.UVData()
print '     Loading data'
fhdlist = glob.glob(opts.datpath+'vis_data/'+obsid+'*') + glob.glob(opts.datpath+'metadata/'+obsid+'*')
uv.read_fhd(fhdlist,run_check=False,run_check_acceptability=False)
print '     Loading mdlvis'
npz_x = np.load(opts.solpath+obsid+'.xx.omni.npz')
npz_y = np.load(opts.solpath+obsid+'.yy.omni.npz')
ant = []
a1, a2 = [], []
mdvis = {'xx':{}, 'yy':{}}
ubls = []
chix = npz_x['chisq2']
chiy = npz_y['chisq2']
sample = np.ones(chix.shape)*16
for ii in range(384):
    if ii%16 == 8: sample[:,ii] = 8
maskx = chix > 1.2
masky = chiy > 1.2
flag = {}
flag['xx'] = np.logical_or(npz_x['flags'], maskx)
flag['yy'] = np.logical_or(npz_y['flags'], masky)
mask = {'xx': npz_x['flags'], 'yy': npz_y['flags']}
for key in npz_x.keys():
    if key.startswith('<'):
        bl,pol = key.split()
        bl = tuple(map(int,bl[1:-1].split(',')))
        mdvis[pol][bl] = npz_x[key]
        i,j = bl
        if not i in ant: ant.append(i)
        if not j in ant: ant.append(j)
        ubls.append(bl)
for key in npz_y.keys():
    if key.startswith('<'):
        bl,pol = key.split()
        bl = tuple(map(int,bl[1:-1].split(',')))
        mdvis[pol][bl] = npz_y[key]
ant.sort()
ubls.sort()
uv.antenna_numbers = np.arange(len(ant))
ant_dict = {}
for ii in range(len(ant)): ant_dict[ant[ii]] = ii
for bl in ubls:
    i,j = bl
    a1.append(ant_dict[i])
    a2.append(ant_dict[j])
Nbls0 = uv.Nbls
Nbls1 = len(ubls)
b0 = 128*uv.ant_1_array[:Nbls0] + uv.ant_2_array[:Nbls0]
uv.Nbls = Nbls1
uv.Nants_data = len(ant)
uv.Nants_telescope = len(ant)
uv.Nblts = uv.Ntimes*uv.Nbls
times = np.resize(np.unique(uv.time_array),(uv.Nbls,uv.Ntimes)).T
uv.time_array = np.resize(times,(times.size))
lsts = np.resize(np.unique(uv.lst_array),(uv.Nbls,uv.Ntimes)).T
uv.lst_array = np.resize(lsts,(lsts.size))
uvw = np.zeros((uv.Nblts,3))
uv.ant_1_array = np.array(a1*uv.Ntimes)
uv.ant_2_array = np.array(a2*uv.Ntimes)
b1 = 128*uv.ant_1_array[:Nbls1] + uv.ant_2_array[:Nbls1]
for ii in range(uv.Nbls):
    i = b1[ii]/128
    j = b1[ii]%128
    try:
        ind = np.where(b0 == 128*ant[i] + ant[j])[0][0]
        uvw[ii::Nbls1] = uv.uvw_array[ind::Nbls0]
    except:
        ind = np.where(b0 == 128*ant[j] + ant[i])[0][0]
        uvw[ii::Nbls1] = -uv.uvw_array[ind::Nbls0]
uv.uvw_array = uvw
uv.nsample_array = np.zeros((uv.Nblts,uv.Nspws,uv.Nfreqs,uv.Npols))
uv.data_array = np.zeros((uv.Nblts,uv.Nspws,uv.Nfreqs,uv.Npols),dtype=np.complex64)
uv.flag_array = np.ones((uv.Nblts,uv.Nspws,uv.Nfreqs,uv.Npols),dtype=bool)
uv.baseline_array = uv.antnums_to_baseline(uv.ant_1_array,uv.ant_2_array)
uv.antenna_positions = np.zeros((len(ant),3))
uv.antenna_names = []
for ii in ant: uv.antenna_names.append(str(ii))
for pp in ['xx','yy']:
    pn = aipy.miriad.str2pol[pp]
    pid = np.where(uv.polarization_array==pn)[0][0]
    for ii in range(uv.Nbls):
        i = a1[ii]
        j = a2[ii]
        uv.data_array[:,0][:,:,pid][ii::uv.Nbls] = mdvis[pp][(ant[i],ant[j])]
        uv.flag_array[:,0][:,:,pid][ii::uv.Nbls] = flag[pp]
        uv.nsample_array[:,0][:,:,pid][ii::uv.Nbls] = sample*np.logical_not(mask[pp])
outuvfits = opts.outpath + obsid + '_mvis_unique.uvfits'
print '     Writing ' + outuvfits
uv.write_uvfits(outuvfits,spoof_nonessential=True)


