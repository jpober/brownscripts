import sys, optparse
import numpy as np

o = optparse.OptionParser()
o.set_usage('eva_chi.py [options] obsfile')
o.set_description(__doc__)
o.add_option('-u',dest='up',default=1.15,type='float',help='upper threshold for chi-square/dof')
o.add_option('-l',dest='low',default=0.9,type='float',help='lower threshold for chi-square/dof')
o.add_option('-o',dest='omnipath',type='string',default='/users/wl42/data/wl42/OBS0/sol_chi/',help='path to omnical solutions')
o.add_option('-n',dest='nsamp',default=50,type='int',help='Number of bad samples threshold for flagging obs')
opts,args = o.parse_args(sys.argv[1:])

fn = open(args[0],'rb')
obs = []
for line in fn: obs.append(line.strip())
bad_obs = []
for o in obs:
    try: 
         npzx = np.load(opts.omnipath+o+'.xx.omni.npz')
         npzy = np.load(opts.omnipath+o+'.yy.omni.npz')
    except:
         print "obs "+o+" failed"
         bad_obs.append(o)
         continue
    indx = np.where(npzx['flags']==False)
    indy = np.where(npzy['flags']==False)
    chix = npzx['chisq2'][indx]
    chiy = npzy['chisq2'][indy]
    nx = np.where(chix>opts.up)[0].size
    ny = np.where(chiy>opts.up)[0].size
    if nx > opts.nsamp or ny > opts.nsamp: 
        bad_obs.append(o)
        print o, 'high',nx, ny
#    chix += npzx['flags']
#    chiy += npzy['flags']
    nx = np.where(chix<opts.low)[0].size
    ny = np.where(chiy<opts.low)[0].size
    if nx > opts.nsamp or ny > opts.nsamp:
        if not o in bad_obs: bad_obs.append(o)
        print o, 'low', nx, ny
print 'Number of flagged obs:', len(bad_obs)
for o in bad_obs: print o
#print bad_obs


