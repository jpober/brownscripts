import numpy as np, glob, sys, copy
import matplotlib.pyplot as plt
#exec('from plot_vis import *')
obs = sys.argv[1]
pol = sys.argv[2]
fn = glob.glob('./sol_chi/'+obs+'*'+pol+'.omni.npz')
fn.sort()
n = 0
for f in fn:
    d = np.load(f)
    if n == 0:
        dd = d['chisq2']
        mm = d['flags']
    else:
        dd = np.concatenate((dd,d['chisq2']))
        mm = np.concatenate((mm,d['flags']))
    n += 1
ind = np.where(mm<1)
x = dd[ind]
#ind = np.where(x>1.2)
#x[ind] = 1.2
n, bins, patches = plt.hist(x, 500, normed=1, facecolor='green', alpha=0.75)
#print x.shape, n , bins
plt.xlabel('chi-square/DoF')
plt.ylabel('normed-counts')
plt.show()
