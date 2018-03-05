import numpy as np, glob, sys, copy
exec('from plot_vis import *')
obs = sys.argv[1]
pol = sys.argv[2]
fn = glob.glob('./sol_chi/'+obs+'*'+pol+'.omni.npz')
fn.sort()
n = 0
#m = np.zeros((56,384),dtype=bool)
#m[0]=True;m[-1]=True;m[-2]=True
#for ii in range(384):
#    if ii%16 in [0,15]: m[:,ii]=True
for f in fn:
    d = np.load(f)
    if n == 0:
        dd = d['chisq2']
        mm = d['flags']
    else:
        dd = np.concatenate((dd,d['chisq2']))
        mm = np.concatenate((mm,d['flags']))
    n += 1
plot_chisq(dd,mm)

