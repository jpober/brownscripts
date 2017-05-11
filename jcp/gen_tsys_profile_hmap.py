#! /usr/bin/env python

import numpy as np
import aipy as a
import pylab as plt
from mpl_toolkits.basemap import Basemap

beamfile = '/Users/jpober/science/hera/beams/GX4Y2H_4900_150.hmap'
#skyfile = '/Users/jpober/science/hera/global/dOC_gsm_150MHz.fits'
skyfile = '/Users/jpober/Dropbox/manuscripts/delayspec/plots/haslam408.fits'

beam = a.map.Map(fromfits=beamfile)
sky = a.map.Map(fromfits=skyfile)

#sky map mapping
im = a.img.Img(size=100, res=.5)
tx,ty,tz = im.get_top()
invalid = tx.mask.copy()

#beam map mapping
coords = Basemap(projection='ortho',lat_0=90,lon_0=180,rsphere=1.)
lons,lats,x,y = coords.makegrid(200,200, returnxy=True)
lons = 360 - lons
lats *= a.img.deg2rad; lons *= a.img.deg2rad
y,x,z = a.coord.radec2eq(np.array([lons.flatten(), lats.flatten()]))
beammap = np.real(beam[x,y,z]*np.conj(beam[x,y,z]))
#beammap = np.real((beam[x,y,z]/beam[0,0,1])*np.conj(beam[x,y,z]/beam[0,0,1]))
beammap.shape = (200,200)
beammap = np.where(a.img.recenter(invalid,(100,100)), 0, beammap)
plt.imshow(np.log10(beammap));plt.colorbar();plt.show()
#plt.imshow(beammap);plt.colorbar();plt.show()
#beammap = np.ones_like(beammap)

ha_range = np.arange(0,2*np.pi, .01)
tsky_profile = []
for ha in ha_range:
    ex, ey, ez = im.get_eq(ra=ha, dec=-0.53581608)
    ex = ex.filled(1).flatten() 
    ey = ey.filled(0).flatten()
    ez = ez.filled(0).flatten()
    tsky = sky[ex,ey,ez]*(408/150)**2.55 
    tsky.shape = invalid.shape
    tsky = np.where(invalid, 0, tsky)
    tsky = a.img.recenter(tsky,(100,100))
    #plt.imshow(tsky*beammap);plt.colorbar();plt.show()

    tsky_profile.append(np.sum(tsky*beammap)/np.sum(beammap))

plt.plot(ha_range,tsky_profile)
plt.show()
