import numpy as np, sys
import pyuvdata.uvdata as uvd
obs = sys.argv[1]
uv = uvd.UVData()
print 'loading...'
uv.read_uvfits(obs+'.uvfits',run_check=False,run_check_acceptability=False)
print 'averaging..'
data = (uv.data_array*np.logical_not(uv.flag_array)).reshape(uv.Nblts,uv.Nspws,uv.Nfreqs/2,2,uv.Npols)
wgts = np.logical_not(uv.flag_array.reshape(uv.Nblts,uv.Nspws,uv.Nfreqs/2,2,uv.Npols))
data = np.sum(data,axis=3)
wgts = np.sum(wgts,axis=3)
uv.flag_array = (wgts<=0)
uv.nsample_array = np.float32(wgts*8)
ind = np.where(wgts==0)
wgts[ind] = -1
uv.data_array = data/wgts
fq = uv.freq_array.reshape(uv.Nspws,uv.Nfreqs/2,2)
uv.freq_array = np.mean(fq,axis=2)
uv.Nfreqs /= 2
uv.channel_width *= 2
del data, wgts
print 'writing...'
uv.write_uvfits('./uvfitsdata/+'obs+'.uvfits',spoof_nonessential=True)

