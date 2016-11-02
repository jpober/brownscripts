import numpy as n
import uvdata
from glob import glob
import sys,os,optparse

o = optparse.OptionParser()
opts,args = o.parse_args(sys.argv[1:])

#print args
#obsDirty = glob('Pzen*HPAB')
#obsResidual = glob('Pzen*SPAB')
obsSum = {}
res={}
dirty={}
chans = ['30_50','95_115','135_155']
for i in args:
    print i
    a = uvdata.miriad.Miriad()
    try:
        a.read_miriad(i)
    except KeyError:
        for j in chans:
            ch = j.split('_')
            if i.split('.uvcRREcACOTUc')[1].split('PA')[0]=='H':
                dirty[j] = a.data_array[:,0,int(ch[0]):int(ch[1]),0].sum()
            else:
                res[j] = a.data_array[:,0,int(ch[0]):int(ch[1]),0].sum()
        continue
res['time'] = 24.0*n.mean(a.lst_array)/(2*n.pi)
dirty['time'] = 24.0*n.mean(a.lst_array)/(2*n.pi)

obsSum['res'] = res
obsSum['dirty'] = dirty
n.save(i.split('Pzen.')[1].split('.uv')[0],obsSum)
