import numpy as n
import uvdata
import aipy as a
import pylab as pl
from glob import glob
import optparse,sys

#obs = glob('./even_LST/*.uvH')

o = optparse.OptionParser()
opts, args = o.parse_args(sys.argv[1:])
uv = a.miriad.UV(args[0])
aa = a.cal.get_aa('psa6240_FHD',uv['sdf'],uv['sfreq'],uv['nchan'])
del(uv)

uv = uvdata.miriad.Miriad()
uv.read_miriad(args[0])

bls = []
blnum = []
for i in range(0,64):
    for j in range(0,64):
        bls.append(str(j)+'_'+str(i))
        blnum.append(uv.antnums_to_baseline(i,j))
#print len(bls)
#print bls
#bls = n.array(bls)
#bls = ('41_49','3_49','47_49','25_49','1_49','48_49','4_49','18_49','29_49','24_49','17_49')
bl_len = []
for b in bls:
    i,j = b.split('_')
    ix = aa.get_arr_params()['antpos'][int(i)]['top_x']/100.0
    iy = aa.get_arr_params()['antpos'][int(i)]['top_y']/100.0
    jx = aa.get_arr_params()['antpos'][int(j)]['top_x']/100.0
    jy = aa.get_arr_params()['antpos'][int(j)]['top_y']/100.0
    bl_len.append(round(n.sqrt((ix - jx)**2 + (iy-jy)**2),1))
bl_len = n.array(bl_len)
#print bls[bl_len==40.0]
blnum = n.array(blnum)
uni_bl = n.unique(bl_len)
uni_bl = uni_bl[uni_bl > 25.0]
#uv = uvdata.miriad.Miriad()
#uv.read_miriad(args[0])

bsls = []
for i in uni_bl:
    print (i==bl_len)
    blnum_id = blnum[i == bl_len]
    avg_all = n.zeros((14,203))
    for j in blnum_id:
        bsls = uv.data_array[uv.baseline_array==j,0,:,0]
        if n.shape(bsls) == (0,203):
            continue
        else:
            avg_all = n.mean((bsls,avg_all),axis=0)
    for k in blnum_id:
        if avg_all.shape == (0,203) or uv.data_array[uv.baseline_array==k,0,:,0].shape == (0,203):
            continue
        else:
            print avg_all.shape
            print uv.data_array[uv.baseline_array==k,0,:,0].shape
            uv.data_array[uv.baseline_array==k,0,:,0] = avg_all

uv.write_miriad('AveragedOverSameSeps.uv')
