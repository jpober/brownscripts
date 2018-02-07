import numpy as np
import matplotlib.gridspec as gridspec
import time

from matplotlib.pyplot import *
from mpl_toolkits.axes_grid1 import make_axes_locatable

M = int(np.sqrt(20000))
N = int(np.sqrt(20000))

xs = np.linspace(0,10,M)
ys = np.linspace(0,10,N)
extent = [xs.min(), xs.max(), ys.min(), ys.max()]

ks_x = np.fft.fftshift(np.fft.fftfreq(len(xs), d=np.mean(np.diff(xs))))
ks_y = np.fft.fftshift(np.fft.fftfreq(len(ys), d=np.mean(np.diff(ys))))
extent_ks = [ks_x.min(), ks_x.max(), ks_y.min(), ks_y.max()]

X,Y = np.meshgrid(xs, ys)
Z = np.exp(-((X-5)**2+(Y-5)**2))

Xd,KX = np.meshgrid(xs, ks_x)
Yd,KY = np.meshgrid(ys, ks_y)

DFT_x = np.exp(-1j*2*np.pi*Xd*KX)
DFT_y = np.exp(-1j*2*np.pi*Yd*KY)

start = time.time()
ftZ = np.dot(np.dot(DFT_x, Z), DFT_y.T)
stop = time.time()
mytime = stop - start

start = time.time()
npftZ = np.fft.fftshift(np.fft.fft2(Z))
stop = time.time()
nptime = stop - start

print 'My FFT: %.1e' %mytime
print 'np FFT: %.1e' %nptime

fig = figure(figsize=(16,4))
gs = gridspec.GridSpec(1,4)

fmt = 'o-'
imgs = []

ax = subplot(gs[0])
ogim = ax.imshow(Z,
                            interpolation='nearest',
                            origin='center',
                            extent=extent)
ax.set_title('original')
ax.set_xlabel('x', size=16)
ax.set_ylabel('y', size=16)

myax = subplot(gs[1])
myim = myax.imshow(np.abs(ftZ),
                                 interpolation='nearest',
                                 origin='center',
                                 extent=extent_ks)
myax.set_title('my fft')

npax = subplot(gs[2])
npim = npax.imshow(np.abs(npftZ),
                                interpolation='nearest',
                                origin='center',
                                extent=extent_ks)
npax.set_title('numpy fft')

dax = subplot(gs[3])
dim = dax.imshow(np.abs(ftZ)-np.abs(npftZ),
                            interpolation='nearest',
                            origin='center',
                            extent=extent_ks)
dax.set_title('my fft - numpy fft')

for a in fig.axes[1:]:
    a.set_xlabel('k_x', size=16)
    a.set_ylabel('k_y', size=16)

imgs = [ogim, myim, npim, dim]

for i,axe in enumerate(fig.axes):
    divider = make_axes_locatable(axe)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    colorbar(imgs[i], cax=cax)


gs.tight_layout(fig)

show()
