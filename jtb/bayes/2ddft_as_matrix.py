import numpy as np
import matplotlib.gridspec as gridspec
import time

from matplotlib.pyplot import *
from mpl_toolkits.axes_grid1 import make_axes_locatable

M = 100
N = 100

# Construct sky
ls = np.linspace(0,10,M)
ms = np.linspace(0,10,N)
N_im = ls.size*ms.size
extent_lm = [ls.min(), ls.max(), ms.min(), ms.max()]

L, M = np.meshgrid(ls, ms)
I = np.exp(-(L-5)**2-(M-5)**2)

# Construct position vectors for N_im pixels
I_vec = np.reshape(I, (N_im,), order='F')
ls_vec, ms_vec = np.zeros(0), np.zeros(0)
for l in ls:
    for m in ms:
        ls_vec = np.append(ls_vec, l)
        ms_vec = np.append(ms_vec, m)

# Construct uv points
us = np.fft.fftshift(np.fft.fftfreq(len(ls), d=np.mean(np.diff(ls))))
vs = np.fft.fftshift(np.fft.fftfreq(len(ms), d=np.mean(np.diff(ms))))
N_vis = us.size*vs.size
extent_uv = [us.min(), us.max(), vs.min(), vs.max()]

# Construct position vectors for N_vis pixels
us_vec, vs_vec = np.zeros(0), np.zeros(0)
for u in us:
    for v in vs:
        us_vec = np.append(us_vec, u)
        vs_vec = np.append(vs_vec, v)

# Construct DFT matrix
Ls_m, Us_m = np.meshgrid(ls_vec, us_vec)
Ms_m, Vs_m = np.meshgrid(ms_vec, vs_vec)
DFT = np.exp(-1j*2*np.pi*(Us_m*Ls_m+Vs_m*Ms_m))

start = time.time()
Vs_vec = np.dot(DFT, I_vec)
Vs = np.reshape(Vs_vec, (us.size, vs.size), order='F')
stop = time.time()
mytime = stop - start

print I.shape

start = time.time()
npVs = np.fft.fftshift(np.fft.fft2(I))
stop = time.time()
nptime = stop - start

print npVs.shape

# print 'DFT: ',
# print DFT.shape
#
# print 'Vs: ',
# print Vs.shape
#
# print 'npVs: ',
# print npVs.shape

# print 'My FFT: %.1e' %mytime
# print 'np FFT: %.1e' %nptime

fig = figure(figsize=(16,4))
gs = gridspec.GridSpec(1,4)

fmt = 'o-'
imgs = []
aspect = 'equal'

ax = subplot(gs[0])
ogim = ax.imshow(I,
                            interpolation='nearest',
                            origin='center',
                            extent=extent_lm,
                            aspect=aspect)
ax.set_title('original')
ax.set_xlabel('l', size=16)
ax.set_ylabel('m', size=16, rotation=0)

myax = subplot(gs[1])
myim = myax.imshow(np.abs(Vs),
                                 interpolation='nearest',
                                 origin='center',
                                 extent=extent_uv,
                                 aspect=aspect)
myax.set_title('my fft')

npax = subplot(gs[2])
npim = npax.imshow(np.abs(npVs),
                                interpolation='nearest',
                                origin='center',
                                extent=extent_uv,
                                aspect=aspect)
npax.set_title('numpy fft')

dax = subplot(gs[3])
dim = dax.imshow(np.abs(Vs)-np.abs(npVs),
                            interpolation='nearest',
                            origin='center',
                            extent=extent_uv,
                            aspect=aspect)
dax.set_title('my fft - numpy fft')

for a in fig.axes[1:]:
    a.set_xlabel('u', size=16)
    a.set_ylabel('v', size=16, rotation=0)

imgs = [ogim, myim, npim, dim]

for i,axe in enumerate(fig.axes):
    divider = make_axes_locatable(axe)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    colorbar(imgs[i], cax=cax)


gs.tight_layout(fig)

show()
