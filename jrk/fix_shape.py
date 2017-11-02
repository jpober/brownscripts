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

for filename in args:
    #import observation
    a = uvdata.miriad.Miriad()
    a.read_miriad(filename)
    
    print 'Saving...'+filename.split('/')[-1:][0]+'--->'+filename.split('/')[-1:][0]+'.uvfits'
    a.write_miriad(filename.split('/')[-1:][0])
