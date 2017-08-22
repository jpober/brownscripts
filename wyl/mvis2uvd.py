import sys,optparse,aipy
import numpy as np, mp2cal 
import pyuvdata.uvdata as uvd

o = optparse.OptionParser()
o.set_usage('mvis2uvd.py [options] obsid') #only takes 1 obsid
o.set_description(__doc__)
o.add_option('-d',dest='datpath',default='/users/wl42/data/wl42/Nov2016EoR0/',type='string', help='Path to data. Include final / in path.')
o.add_option('-s',dest='solpath',default='/users/wl42/data/wl42/Nov2016EoR0/mdl_sol/',type='string', help='Path to omnical solutions. Include final / in path.')
o.add_option('-o',dest='outpath',default='/users/wl42/data/wl42/MDLVIS/',type='string', help='Path to save uvfits. Include final / in path.')
opts,args = o.parse_args(sys.argv[1:])
exec('from PhaseII_cal import *')
obsid = args[0]
uv = uvd.UVData()
print '     Loading data'
uv.read_uvfits(opts.datpath+obsid+'.uvfits',run_check=False,run_check_acceptability=False)
print '     Loading mdlvis'
npz_x = np.load(opts.solpath+obsid+'.xx.omni.npz')
npz_y = np.load(opts.solpath+obsid+'.yy.omni.npz')
mdvis = {'xx':{}, 'yy':{}}
reds = mp2cal.wyl.cal_reds_from_pos(antpos)
reds_ind = {}
ubls = []
chix = npz_x['chisq']
chiy = npz_y['chisq']
maskx = np.zeros(chix.shape)
masky = np.zeros(chiy.shape)
indx = np.where(chix < 1.2)
indy = np.where(chiy < 1.2)
maskx[indx] = 1
masky[indy] = 1
for key in npz_x.keys():
    if key.startswith('<'):
        bl,pol = key.split()
        bl = tuple(map(int,bl[1:-1].split(',')))
        mdvis[pol][bl] = npz_x[key]*maskx
        ubls.append(bl)
for key in npz_y.keys():
    if key.startswith('<'):
        bl,pol = key.split()
        bl = tuple(map(int,bl[1:-1].split(',')))
        mdvis[pol][bl] = npz_y[key]*masky
uv.data_array *= 0
for r in reds:
    ubl = None
    for bl in ubls:
        if bl in r: ubl = bl
    if ubl is None: continue
    for b in r: reds_ind[b] = ubl
a1 = uv.ant_1_array[:uv.Nbls]
a2 = uv.ant_2_array[:uv.Nbls]
for pp in ['xx','yy']:
    pn = aipy.miriad.str2pol[pp]
    pid = np.where(uv.polarization_array==pn)[0][0]
    for ii in range(uv.Nbls):
        i = a1[ii]
        j = a2[ii]
        if (i,j) in reds_ind.keys(): vis = mdvis[pp][reds_ind[(i,j)]]
        elif (j,i) in reds_ind.keys(): vis = mdvis[pp][reds_ind[(j,i)]].conj()
        else: continue
        uv.data_array[:,0][:,:,pid][ii::uv.Nbls] = vis
ind = np.where(uv.data_array == 0)
uv.flag_array[ind] = True
outuvfits = opts.outpath + obsid + '_mvis.uvfits'
print '     Writing ' + outuvfits
uv.write_uvfits(outuvfits,spoof_nonessential=True)


