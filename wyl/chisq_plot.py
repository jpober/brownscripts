import numpy as np, glob, sys, copy, optparse

o = optparse.OptionParser()
o.set_usage('chisq_plot.py [options] obsfile')
o.set_description(__doc__)
o.add_option('-o',dest='omnipath',type='string',default='/users/wl42/data/wl42/OBS0/sol_chi/',help='path to omnical solutions')
o.add_option('-p',dest='pol',type='string',default='xx',help='polarization')
opts,args = o.parse_args(sys.argv[1:])
exec('from plot_vis import *')
obs = args[0]
pol = opts.pol
fn = glob.glob(opts.omnipath+obs+'*'+pol+'.omni.npz')
#print opts.omnipath+obs+'*'+pol+'.omni.npz'
#print len(fn)
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

