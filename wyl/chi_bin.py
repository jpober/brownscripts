import numpy as np, glob, sys, copy, optparse
import matplotlib.pyplot as plt
#exec('from plot_vis import *')
#obs = sys.argv[1]
#pol = sys.argv[2]
#fn = glob.glob('./mdl_sol/'+obs+'*'+pol+'.omni.npz')
o = optparse.OptionParser()
o.set_usage('chi_bin.py [options] obsfile')
o.set_description(__doc__)
o.add_option('-o',dest='omnipath',type='string',default='/users/wl42/data/wl42/OBS0/sol_chi/',help='path to omnical solutions')
o.add_option('-p',dest='pol',type='string',default='xx',help='polarization')
opts,args = o.parse_args(sys.argv[1:])
exec('from plot_vis import *')
obs = args[0]
pol = opts.pol
fn = glob.glob(opts.omnipath+obs+'*'+pol+'.omni.npz')

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
ind = np.where(x>5)
x[ind] = 5.
n, bins, patches = plt.hist(x, 500, normed=1, facecolor='green', alpha=0.75)
#print x.shape, n , bins
plt.xlabel('chi-square/DoF')
plt.ylabel('normed-counts')
plt.show()
