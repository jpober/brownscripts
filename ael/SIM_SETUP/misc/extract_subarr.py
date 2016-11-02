#!/usr/bin/python

import aipy as a
import numpy as n
import sys

aa = a.cal.get_aa(sys.argv[1], .1, .1, 1)

## Get rows of array
y = 0
rows = []
m = 0
for ant in aa:
    pt = ant.pos
#    print y - pt[1]
    if n.abs(y - pt[1]) > 100:
         rows.append([])
         rows[-1].append(m)
    else:
         rows[-1].append(m)
    m += 1
    y = pt[1]

#print rows
#print rows[-1][0:3]   #3
#print rows[-2][0:4]   #4
#print rows[-3][0:5]   #5
#print rows[-4][0:4]   #4
#print rows[-1][0:3]   #3

ar19 = n.array( rows[-1][0:3] + rows[-2][0:4] + rows[-3][0:5] + rows[-4][1:5] + rows[-5][2:5] ).flatten()
ar37 = n.array( rows[-1][0:4] + rows[-2][0:5] + rows[-3][0:6] + rows[-4][0:7] 
			+ rows[-5][1:7] + rows[-6][2:7] + rows[-7][3:7] ).flatten()

heracal = __import__(sys.argv[1])
prms = heracal.prms['antpos']

ar19_antpos = {i: prms[ant] for (i,ant) in enumerate(ar19)}
ar37_antpos = {i: prms[ant] for (i,ant) in enumerate(ar37)}
print len(ar19_antpos), len(ar37_antpos)


print ar19_antpos

print ' '

print ar37_antpos

