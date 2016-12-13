#!/bin/env python

from scipy.io import readsav
import matplotlib.pyplot as plt
import numpy as np
import optparse, aipy as a
import sys

o = optparse.OptionParser()
o.set_usage('plot_sav.py [options] <obs_id>')
o.set_description(__doc__)
a.scripting.add_standard_options(o, src=True,ant=True, pol=True, chan=True, dec=True,
    cmap=True, max=True, drng=True, cal=True)
o.add_option('-m', '--mode', dest='mode', default='log',
    help='Plot mode can be log (logrithmic), lin (linear), phs (phase), real, or imag.')
o.add_option('-d', '--delay', dest='delay', action='store_true',
    help='Take FFT of frequency axis to go to delay (t) space.')


opts, args = o.parse_args(sys.argv[1:])

xx, yy = {}, {}
jdates = {}
pairs, nants = [], 0
for f in args:
	print 'Reading '+str(f)
	xx[f] = readsav(f+"_vis_XX.sav")
	yy[f] = readsav(f+"_vis_YY.sav")
	bl_info = xx[f]['obs']['BASELINE_INFO'][0]
	jdates[f]=bl_info['jdate']
	if len(pairs) == 0:
		ant_1_array  = bl_info['TILE_A'][0] - 1
		ant_2_array  = bl_info['TILE_B'][0] - 1
		nants = max(ant_1_array) + 1
		pairs = np.array(zip(ant_1_array, ant_2_array))
#a.scripting.parse_ants(opts.ant,nants=nants)
bls = []
searches={}
for bl in opts.ant.split(','):
	bstr = bl
	bl = map(int, bl.split('_'))
	assert(len(bl) in [1,2])
	if len(bl) == 1: continue
	else: bls.append(bstr)
	test = pairs==tuple(sorted(bl,reverse=True))
	searches[bstr] = np.array([all(en) for en in test])

nbls=len(bls)
#print searches

#Sort files into temporal order
j0s = np.array([[f,j[0][0]] for f,j in jdates.iteritems()])
j0s = j0s[j0s[:,1].argsort()]
args = j0s[:,0]
j0s = j0s[:,1]

#Stitch together the data from the various files.

wfalls = {'xx':{}, 'yy': {}}
for b in searches.keys():
   dx, dy = np.array([]), np.array([])
   for f in args:
	sel_xx = xx[f].vis_ptr[searches[b],:]
	sel_yy = yy[f].vis_ptr[searches[b],:]
	if opts.delay:
	   sel_xx=np.fft.ifft(sel_xx)
	   sel_xx=np.fft.fftshift(sel_xx,axes=1)
	   sel_yy=np.fft.ifft(sel_yy)
	   sel_yy=np.fft.fftshift(sel_yy,axes=1)
	if dx.size == 0: dx, dy = sel_xx, sel_yy
	else:
	   dx = np.concatenate((dx, sel_xx),axis=0)
	   dy = np.concatenate((dy, sel_yy),axis=0)
   wfalls['xx'][b] = dx
   wfalls['yy'][b] = dy

npols=2
pols = ['xx','yy']

fig, axes = plt.subplots(nrows=nbls, ncols=npols,squeeze=False)
for bli,bl in enumerate(bls):
   for pol in range(npols):
	ax = axes[bli,pol]
	p = pols[pol]
	if opts.mode == 'lin': ax.imshow(np.abs(wfalls[p][bl]),interpolation='nearest', aspect='auto')
	elif opts.mode == 'phs': ax.imshow(np.angle(wfalls[p][bl]),interpolation='nearest', aspect='auto')
	else: ax.imshow(np.log10(np.abs(wfalls[p][bl])),interpolation='nearest', aspect='auto')
	ax.set_title(bl+p)
	if bli < nbls-1: ax.get_xaxis().set_visible(False)
	if pol > 0: ax.get_yaxis().set_visible(False)
	#plt.colorbar(shrink=0.5)

plt.show()#block=True)

#test = np.array([all(en) for en in test])
#
#
#xx = xx.vis_ptr[test,:]
#yy = yy.vis_ptr[test,:]
#
##nblts=8256
#
##sel_xx=xx[100::nblts,:]
##sel_yy=yy[100::nblts,:]
#
#fig = plt.figure()
#
#ax1 = fig.add_subplot(121)
#ax2 = fig.add_subplot(122)
#
###if opts.delay:
#xx=np.fft.ifft(xx)
#xx=np.fft.fftshift(xx,axes=1)
#yy=np.fft.ifft(yy)
#yy=np.fft.fftshift(yy,axes=1)
#
#ax1.imshow(np.log10(np.abs(xx)),interpolation='nearest',aspect='auto')
#ax2.imshow(np.log10(np.abs(yy)),interpolation='nearest',aspect='auto')
#
#plt.show()
