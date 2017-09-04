import numpy as np
import matplotlib.pyplot as plt
import heracal, mp2cal, sys, optparse
from astropy.io import fits
from scipy.io.idl import readsav

o = optparse.OptionParser()
o.set_usage('redplot.py [options] obsid')
o.add_option('-f', dest='f', help='frequency channel',default='356')
opts,args = o.parse_args(sys.argv[1:])
obsid = args[0]

fqs = int(opts.f)
fhd_path = '/users/wl42/data/wl42/FHD_out/fhd_PhaseII_Longrun_EoR0/calibration/'
omn_path = '/users/wl42/data/wl42/Nov2016EoR0/omni_test/'

def baseline_to_antnums(baseline): # from pyuvdata
    if np.min(baseline) > 2**16:
        ant2 = (baseline - 2**16) % 2048 - 1
        ant1 = (baseline - 2**16 - (ant2 + 1)) / 2048 - 1
    else:
        ant2 = (baseline) % 256 - 1
        ant1 = (baseline - (ant2 + 1)) / 256 - 1
    return np.int32(ant1), np.int32(ant2)

g_fhd = {'x':{}, 'y':{}}
g_omn = {'x':{}, 'y':{}}
fhd_cal = readsav(fhd_path+obsid+'_cal.sav',python_dict=True)
for a in range(fhd_cal['cal']['N_TILE'][0]):
    g_fhd['x'][a] = fhd_cal['cal']['GAIN'][0][0][a] + fhd_cal['cal']['GAIN_RESIDUAL'][0][0][a]
    g_fhd['y'][a] = fhd_cal['cal']['GAIN'][0][1][a] + fhd_cal['cal']['GAIN_RESIDUAL'][0][1][a]
ox_cal = np.load(omn_path+obsid+'.xx.omni.npz')
oy_cal = np.load(omn_path+obsid+'.yy.omni.npz')
for k in ox_cal.keys():
    if k[0].isdigit(): g_omn['x'][int(k[:-1])] = ox_cal[k]
for k in oy_cal.keys():
    if k[0].isdigit(): g_omn['y'][int(k[:-1])] = oy_cal[k]

F = fits.open(obsid+'.uvfits')
D = F[0]
try:
    a1 = np.int32(D.data.field('ANTENNA1')) - 1
    a2 = np.int32(D.data.field('ANTENNA2')) - 1
except(KeyError):
    bl_input_array = np.int64(D.data.field('BASELINE'))
    a1,a2 = baseline_to_antnums(bl_input_array)
baseline_ind = 256*a1 + a2
Nbls = len(np.unique(baseline_ind))
Ntimes = baseline_ind.size/Nbls
Nfreqs = D.header['NAXIS4']
Npols = D.header['NAXIS3']
print '   (Nbls,Ntimes,Nfreqs,Npols): ',(Nbls,Ntimes,Nfreqs,Npols)
if D.header['NAXIS'] == 7:
    data = (D.data.field('DATA')[:, 0, 0, :, :, :, 0] + 1j * D.data.field('DATA')[:, 0, 0, :, :, :, 1])[:,0].reshape(Ntimes,Nbls,Nfreqs,Npols)
    wegt = (D.data.field('DATA')[:, 0, 0, :, :, :, 2] > 0)[:,0].reshape(Ntimes,Nbls,Nfreqs,Npols)
else:
    data = (D.data.field('DATA')[:, 0, 0, :, :, 0] + 1j * D.data.field('DATA')[:, 0, 0, :, :, 1]).reshape(Ntimes,Nbls,Nfreqs,Npols)
    wegt = (D.data.field('DATA')[:, 0, 0, :, :, 2] > 0).reshape(Ntimes,Nbls,Nfreqs,Npols)
data = data*wegt
wegt = np.sum(wegt,axis=0)
data = np.sum(data,axis=0)/(wegt+1e-10)
datx = data[:,fqs][:,0]
daty = data[:,fqs][:,1]
a1 = a1[:Nbls]
a2 = a2[:Nbls]
exec("from PhaseII_cal import *")
reds = mp2cal.wyl.cal_reds_from_pos(antpos,ex_ubls=[(57,58)])
group = {}
print "Generating reds:"
for red in reds:
    if len(red) < 40: continue
    group[red[0]] = []
    for ii in range(Nbls):
        bl1 = (a1[ii],a2[ii])
        bl2 = (a2[ii],a1[ii])
        if bl1 in red:
            group[red[0]].append((ii,0))
        elif bl2 in red:
            group[red[0]].append((ii,1))
print len(group.keys())

colorlist = ['red','blue','black','silver','gold','green','violet','lime','brown','indigo','yellow','cyan','magenta','darksalmon','bisque','sandybrown','tan','skyblue','pink','purple']
markerlist = ['o','v','^','<','>','1','2','3','4','8','s','p','*','h','H','+','x','D','d','|','_']
l0 = len(colorlist)
l1 = len(markerlist)

print "Plotting:"
fig = plt.figure()

p0x = fig.add_subplot(2,3,1)
p0x.set_title('Raw Data (xx)',size = 10.0)
#p0x.set_xlabel('vis.real (Jy)')
p0x.set_ylabel('vis.imag (Jy)',size = 10.0)

p1x = fig.add_subplot(2,3,2)
p1x.set_title('After omnical (xx)',size = 10.0)
#p1x.set_xlabel('vis.real (Jy)')
#p1x.set_ylabel('vis.imag (Jy)')

p2x = fig.add_subplot(2,3,3)
p2x.set_title('After FHD (xx)',size = 10.0)
#p2x.set_xlabel('vis.real (Jy)')
#p2x.set_ylabel('vis.imag (Jy)')

p0y = fig.add_subplot(2,3,4)
p0y.set_title('Raw Data (yy)',size = 10.0)
p0y.set_xlabel('vis.real (Jy)',size = 10.0)
p0y.set_ylabel('vis.imag (Jy)',size = 10.0)

p1y = fig.add_subplot(2,3,5)
p1y.set_title('After omnical (yy)',size = 10.0)
p1y.set_xlabel('vis.real (Jy)',size = 10.0)
#p1y.set_ylabel('vis.imag (Jy)')

p2y = fig.add_subplot(2,3,6)
p2y.set_title('After FHD (yy)',size = 10.0)
p2y.set_xlabel('vis.real (Jy)',size = 10.0)
#p2y.set_ylabel('vis.imag (Jy)')

nn=0
maxcx, maxcy, mincx, mincy = 0,0,0,0
for r in group.keys():
    d0x,d0y,d1x,d1y,d2x,d2y = [],[],[],[],[],[]
    for ii in group[r]:
        i,j = (a1[ii[0]],a2[ii[0]])
        if ii[1]==0:
            d0x.append(datx[ii[0]])
            d0y.append(daty[ii[0]])
            d1x.append(datx[ii[0]]/(g_omn['x'][i][fqs]*g_omn['x'][j][fqs].conj()))
            d1y.append(daty[ii[0]]/(g_omn['y'][i][fqs]*g_omn['y'][j][fqs].conj()))
            d2x.append(datx[ii[0]]/(g_fhd['x'][i][fqs]*g_fhd['x'][j][fqs].conj()))
            d2y.append(daty[ii[0]]/(g_fhd['y'][i][fqs]*g_fhd['y'][j][fqs].conj()))
        else:
            d0x.append(datx[ii[0]].conj())
            d0y.append(daty[ii[0]].conj())
            d1x.append(datx[ii[0]].conj()/(g_omn['x'][i][fqs]*g_omn['x'][j][fqs].conj()))
            d1y.append(daty[ii[0]].conj()/(g_omn['y'][i][fqs]*g_omn['y'][j][fqs].conj()))
            d2x.append(datx[ii[0]].conj()/(g_fhd['x'][i][fqs]*g_fhd['x'][j][fqs].conj()))
            d2y.append(daty[ii[0]].conj()/(g_fhd['y'][i][fqs]*g_fhd['y'][j][fqs].conj()))
    d0x = np.array(d0x)
    d0y = np.array(d0y)
    d1x = np.array(d1x)
    d1y = np.array(d1y)
    d2x = np.array(d2x)
    d2y = np.array(d2y)
    tx = np.concatenate((d1x.real,d1x.imag,d2x.real,d2x.imag))
    ty = np.concatenate((d1y.real,d1y.imag,d2y.real,d2y.imag))
    if np.min(tx) < mincx: mincx = np.min(tx)
    if np.max(tx) > maxcx: maxcx = np.max(tx)
    if np.min(ty) < mincy: mincy = np.min(ty)
    if np.max(ty) > maxcy: maxcy = np.max(ty)
    p0x.scatter(d0x.real, d0x.imag, marker=markerlist[nn/l0], color=colorlist[nn%l0],s=0.5)
    p0y.scatter(d0y.real, d0y.imag, marker=markerlist[nn/l0], color=colorlist[nn%l0],s=0.5)
    p1x.scatter(d1x.real, d1x.imag, marker=markerlist[nn/l0], color=colorlist[nn%l0],s=0.5)
    p1y.scatter(d1y.real, d1y.imag, marker=markerlist[nn/l0], color=colorlist[nn%l0],s=0.5)
    p2x.scatter(d2x.real, d2x.imag, marker=markerlist[nn/l0], color=colorlist[nn%l0],s=0.5)
    p2y.scatter(d2y.real, d2y.imag, marker=markerlist[nn/l0], color=colorlist[nn%l0],s=0.5)
    nn += 1

limx = (mincx-0.1*(maxcx-mincx),maxcx+0.1*(maxcx-mincx))
limy = (mincy-0.1*(maxcy-mincy),maxcy+0.1*(maxcy-mincy))
p0x.axis('equal')
p0x.grid(True)
p0x.tick_params(labelsize=8)
p0y.axis('equal')
p0y.grid(True)
p0y.tick_params(labelsize=8)
p1x.set_xlim(limx)
p1x.set_ylim(limx)
#p1x.axis('equal')
p1x.grid(True)
p1x.tick_params(labelsize=8)
p1y.set_xlim(limy)
p1y.set_ylim(limy)
#p1y.axis('equal')
p1y.grid(True)
p1y.tick_params(labelsize=8)
p2x.set_xlim(limx)
p2x.set_ylim(limx)
#p2x.axis('equal')
p2x.grid(True)
p2x.tick_params(labelsize=8)
p2y.set_xlim(limy)
p2y.set_ylim(limy)
#p2y.axis('equal')
p2y.grid(True)
p2y.tick_params(labelsize=8)


plt.subplots_adjust(top=0.9,bottom=0.1,left=0.1,right=0.9)
plt.show()


    
