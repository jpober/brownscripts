import pyuvdata, pylab as plt, numpy as np, sys, aipy as a,
from mpl_toolkits.basemap import Basemap
from scipy import interpolate

args = sys.argv[1:]
x1= 64
x2= 64
#for filename in sys.argv:[]
#print filename

#def makeplot(x1, x2, heradat):
data = numpy.array([])
lst = numpy.array([])

for filename in args:
    uv = pyuvdata.UVData() #uvdata is the package name and UV creates the image
    #uv.read_uvfits(filename)
    uv.read_miriad(filename, run_check=False, run_check_acceptability=False)
    waterfall = (uv.data_array[np.where(uv.baseline_array == uv.antnums_to_baseline(x1, x2))].squeeze().real)
    wslice = waterfall[:, 500]
    data = np.append(data, wslice)
    lstslice = (uv.lst_array[np.where(uv.baseline_array == uv.antnums_to_baseline(x1, x2))].squeeze().real)
    lst= np.append(lst, lstslice)

print data.shape


#download basemap

#download file and replace path
#beamfile = '/Users/jpober/science/hera/beams/GX4Y2H_4900_150.hmap'
beamfile = '/Users/katievasquez/Desktop/pober/GX4Y2H_4900_150.hmap'
#skyfile = '/Users/jpober/science/hera/global/dOC_gsm_150MHz.fits'
#skyfile = '/Users/jpober/Dropbox/manuscripts/delayspec/plots/haslam408.fits'
skyfile = '/Users/katievasquez/Desktop/pober/haslam408.fits'

#pyGSM global sky model, file drop and replace for given sky model.
# Also need to get beam file form github repository
#GX vs GY for different dipoles
#start at 150m Hz -- hard coded

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
#beammap = np.abs(beam[x,y,z]/beam[0,0,1])
beammap = np.real((beam[x,y,z]/beam[0,0,1])*np.conj(beam[x,y,z]/beam[0,0,1]))
beammap.shape = (200,200)
beammap = np.where(a.img.recenter(invalid,(100,100)), 0, beammap)
#plt.imshow(np.log10(beammap));plt.colorbar();plt.show()
#beammap = np.ones_like(beammap)
#plot data curve lst_array -time for x axis- corresponding down select -
# data array and lst array to match then plot two curves on same plot- on my code

ha_range = np.arange(0,2*np.pi, .01)
tsky_profile = []
for ha in ha_range:
    ex, ey, ez = im.get_eq(ra=ha, dec=-0.53581608)
    ex = ex.filled(1).flatten()
    ey = ey.filled(0).flatten()
    ez = ez.filled(0).flatten()
    tsky = sky[ex,ey,ez]*(408/150)**2.55 #scaling factor for mHz
    tsky.shape = invalid.shape
    tsky = np.where(invalid, 0, tsky)
    tsky = a.img.recenter(tsky,(100,100))
    #plt.imshow(tsky*beammap);plt.colorbar();plt.show()

    tsky_profile.append(np.sum(tsky*beammap)/np.sum(beammap))


#create interpolation function here
tsky_profile = np.array(tsky_profile) * 2
interpfit = interpolate.interp1d(ha_range, tsky_profile)
ratio = ((interpfit(lst) + 100)/data)
ratio = np.mean(ratio) #turns array into scalar
print ratio

plt.plot(ha_range,tsky_profile)
plt.plot(lst, ratio*data) #data times number I got from interpolation then ratio then average
plt.show()
