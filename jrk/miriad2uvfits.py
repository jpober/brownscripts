import uvdata
import numpy as n
import optparse,sys
import matplotlib.pyplot as mp
import os
import ephem

#parser = OptionParser()
#parser.add_option("-f", "--file", dest="filename",help="Observation to be converted.", metavar="FILE",action='store')
#(options,args) = parser.parse_args()
o = optparse.OptionParser()
opts,args = o.parse_args(sys.argv[1:])
cwd = os.getcwd()
for filename in args:
    print filename
    #import observation
    a = uvdata.miriad.Miriad()
    try:
        a.read_miriad(filename)
    except KeyError:
        print 'KeyError'

    #time=a.time_array[n.argsort(a.time_array)[len(a.time_array)/2]]
    #a.phase_to_time(a.time_array.mean())
    a.data_array = n.nan_to_num(a.data_array)
    a.flag_array[a.data_array==0j] = True
    ### flag 166:203 and maybe 0:15
    
    #a.flag_array[:,:,0:24,:] = True
    #a.flag_array[:,:,166:203,:] = True 
    
    # flag problem channels
    #for j in range(0,shp[3]):
    #    for i in range(0,shp[2]):
    #        if n.isnan(a.data_array[:,0,i,j].sum())==True:
    #            a.flag_array[:,0,i,j] = True
    
    #Test this method of flagging
    #for i in range(0,shp[3]):
    #    for j in range(0,shp[2]):
    #        for k in range(0,shp[0]):
    #            if n.isnan(a.data_array[k,0,j,i])==True:
    #                a.data_array[k,0,j,i] = n.nan_to_num(a.data_array[k,0,j,i])
    #                a.flag_array[k,0,j,i] = 1   
    # Make sure NaNs don't go unflagged
    #NaNindx = n.where(n.isnan(a.data_array)==True)
    #a.flag_array[NaNindx] = True
    #a.flag_array[a.data_array==0j] = True

    
    os.chdir("/users/jkerriga/data/jkerriga/AnalysisOutput")
    a.phase(dec=a.zenith_dec[0],ra=a.zenith_ra[0],epoch=2000.0)
    #save as uvfits
    print 'Saving...'+filename.split('/')[-1:][0]+'--->'+'P'+filename.split('/')[-1:][0]+'.uvfits'
    a.write_uvfits('P'+filename.split('/')[-1:][0]+'.uvfits',spoof_nonessential=True,force_phase=False)
    del(a)
    os.chdir(cwd)
