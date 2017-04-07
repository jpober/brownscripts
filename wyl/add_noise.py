import numpy, pyuvdata.uvdata as uvd, sys, optparse, aipy
import glob, matplotlib.pyplot as plt
from scipy.io.idl import readsav
from IPython import embed

o = optparse.OptionParser()
o.set_usage('add_noise.py [options] obs')
o.add_option('-o', dest='outpath', help='Destination directory',default='/users/wl42/CALSIM/')
o.add_option('-f', dest='fhdpath', help='FHD directory', default='/users/wl42/FHD_out/fhd_PhaseII_TESTSET_EoR0/')
o.add_option('--gain', dest='gain', default=False, action='store_true', help='')
opts,args = o.parse_args(sys.argv[1:])
obs = args[0]
Trec = 30.
Tsky = 180.
fhdpath = opts.fhdpath
fn = glob.glob(fhdpath+'vis_data/'+obs+'*') + glob.glob(fhdpath+'metadata/'+obs+'*')
uv = uvd.UVData()
uv.read_fhd(fn,use_model=True)
dt = uv.integration_time
df = uv.channel_width
fqs = uv.freq_array[0]/1e6
Tsys = float(opts.Tref)*np.power(fqs/(180.),-2.6) + float(opts.Trec)*np.ones(fqs.shape)
Area = (198000.-215000.)/(200.*200.-150.*150.)*(fqs*fqs-150.*150.)+215000.
sigs =  k*Tsys/(Area*np.sqrt(df*dt))*1e23/np.sqrt(2)

print '   adding noise:'
for ff in range(uv.Nfreqs):
    noise = (np.random.normal(0,sigs[ff],(uv.Nblts,uv.Nspws,uv.Npols))+1j*np.random.normal(0,sigs[ff],(uv.Nblts,uv.Nspws,uv.Npols)))*np.logical_not(uv.flag_array[:,:,ff])
    uv.data_array[:,:,ff] += noise
if opts.gain:
    print '   apply gains:'
    cal = readsav(fhdpath+'calibration/'+obs+'_cal.sav',python_dict=True)
    a1 = uv.ant_1_array[:uv.Nbls]
    a2 = uv.ant_2_array[:uv.Nbls]
    g = {'x':[],'y':[]}
    for i in range(uv.Nbls):
        g['x'].append(cal['cal']['GAIN'][0][0][a1[i]])
        g['y'].append(cal['cal']['GAIN'][0][1][a2[i]])
    g['x'] = g['x']*uv.Ntimes
    g['y'] = g['y']*uv.Ntimes
    g['x'] = np.array(g['x'])
    g['y'] = np.array[g['y']]
    for i in uv.Npols:
        p1,p2 = aipy.miriad.pol2str[uv.polarization_array[i]]
        uv.data_array[:,0][:,:,i] *= (g[p1]*g[p2].conj())
print '   writing...'
uv.write_uvfits(outpath+obs+'.uvfits',spoof_nonessential=True)




