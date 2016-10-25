import numpy as np
import matplotlib.pyplot as plt
import sys, math, glob
import capo.omni as omni
import capo.wyl as wyl

type = sys.argv[2]
if type == 'uvfits':
    uv = wyl.data_uvfits()
    uv.read_data_only(sys.argv[1]+'.uvfits')
elif type == 'fhd':
    fn = glob.glob(sys.argv[1]+'*')
    uv = wyl.data_fhd()
    uv.read_data_only(fn)
a1 = uv.ant_1_array
a2 = uv.ant_2_array
dat = uv.data_array[:,0][:,0][:,0]
exec("from PhaseII_cal_pos import antpos as _antpos")
info = omni.pos_to_info(_antpos)
reds = info.get_reds()
group = {}
print "Generating reds:"
for red in reds:
    if len(red) < 3: continue
    group[red[0]] = []
    for ii in range(0,a1.shape[0]):
        bl1 = (a1[ii],a2[ii])
        bl2 = (a2[ii],a1[ii])
        if bl1 in red: group[red[0]].append((ii,0))
        elif bl2 in red: group[red[0]].append((ii,1))
print len(group.keys())

Nbl = a1.shape[0]
colorlist = ['red','blue','black','silver','gold','green','violet','lime','brown','indigo','yellow']
markerlist = ['o','v','^','<','>','1','2','3','4','8','s','p','*','h','H','+','x','D','d','|','_']
l0 = len(colorlist)
l1 = len(markerlist)
nn = 0
print "Plotting:"
for key in group.keys():
    d = []
    for ii in group[key]:
        if ii[1]==0: d.append(dat[ii[0]])
        else: d.append(dat[ii[0]].conj())
#        if np.abs(dat[ii[0]])>100000: print (a1[ii[0]],a2[ii[0]])
    d = np.array(d)
    if True:#key.len == bl_len[0]:
        ll = plt.scatter(d.real, d.imag, marker=markerlist[nn/l0], color=colorlist[nn%l0])
    nn += 1
#    if nn == l0: break
plt.xlabel('vis.real (Jy)')
plt.ylabel('vis.imag (Jy)')
plt.title('omni_calibrated')
plt.axes().set_aspect('equal', 'datalim')
plt.show()

    
    
