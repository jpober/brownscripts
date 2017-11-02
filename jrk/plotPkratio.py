import numpy as n
import pylab as pl

raw = n.load('/users/jkerriga/brownscripts/jrk/LSTRawNoCov-2DayTest/95_115/I/pspec_LSTRawNoCov-2DayTest_95_115_I.npz')
dirty=n.load('/users/jkerriga/brownscripts/jrk/LSTDirtyNoCov-2DayTest/95_115/I/pspec_LSTDirtyNoCov-2DayTest_95_115_I.npz')

rawk3pk=raw['k3pk']
dirtyk3pk=dirty['k3pk']

ratiok3pk = dirtyk3pk/rawk3pk
k = raw['k']

pl.plot(k,ratiok3pk)
pl.xlabel('k')
pl.ylabel('P(k)/k^3')
pl.savefig('ratio.png')
