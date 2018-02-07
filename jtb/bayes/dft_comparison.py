import numpy as np
import optparse, sys

from matplotlib.pyplot import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec

# Option parser
o = optparse.OptionParser()
o.add_option('--d',
    type=int,
    help='Dimension of DFT.')
opts,args = o.parse_args(sys.argv[1:])

# x = range(10)
# y = range(10)
# X,Y = np.meshgrid(x,y)
# Z = np.zeros_like(X,Y)
# Z[4:7] = 1.0
# ftZ = np.zeros_like(Z)
#
# nx = len(x)
# ny = len(y)

fig = figure(figsize=(16,4))
gs = gridspec.GridSpec(1,4)

# Test 1d
if opts.d == 1:
    xs = np.linspace(-4,4,100)
    ys = np.exp(-xs**2)

    # Perform DFT
    DFT = lambda k, x: np.exp(-1j*2*np.pi*k*x)
    ks = np.sort(np.fft.fftfreq(len(xs), d=np.mean(np.diff(xs))))
    DFT_mat = np.zeros([len(xs),len(xs)])
    print xs.shape, ys.shape
    print DFT_mat.shape
    for i in range(ks.shape[0]):
        for j in range(ys.shape[0]):
            DFT_mat[i,j] = np.exp(-1j*2*np.pi*ks[i]*xs[j])
    fty = np.dot(DFT_mat, ys)

    print ks.shape, fty.shape

    npfty = np.fft.fft(ys)

    fmt = 'o-'

    ax = subplot(gs[0])
    ax.plot(xs, ys, fmt, ms=4)
    ax.set_title('original')
    ax.set_xlabel('x', size=16)
    ax.set_ylabel('y(x)', size=16)

    myax = subplot(gs[1])
    myax.plot(ks, fty, fmt, ms=4)
    myax.set_title('my fft')
    myax.set_ylim([-3,25])
    # myax.imshow(np.abs(fty), interpolation='nearest')
    # ax.colorbar()

    npax = subplot(gs[2])
    npax.plot(ks, np.abs(np.fft.fftshift(npfty)), fmt, ms=4)
    npax.set_title('numpy fft')
    npax.set_ylim([-3,25])
    # npax.imshow(np.abs(np.fft.fft(ys)), interploation='nearest')
    # myax.colorbar()

    dax = subplot(gs[3])
    dax.plot(ks, np.abs(fty) - np.abs(np.fft.fftshift(npfty)), fmt, ms=4)
    dax.set_title('my fft - numpy fft')

    for a in fig.axes[1:]:
        a.set_xlabel('k', size=16)
        a.set_ylabel('y_dft(k)', size=16)

# Test 2d
if opts.d == 2:
    xs = np.linspace(-4,4,25)
    ys = np.linspace(-4,4,25)
    extent = [xs.min(), xs.max(), ys.min(), ys.max()]
    X,Y = np.meshgrid(xs,ys)
    Z = np.exp(-X**2-Y**2)
    ftZ = np.zeros_like(Z)

    # Perform DFT
    DFT = lambda k, x: np.exp(-1j*2*np.pi*k*x)
    ks_x = np.sort(np.fft.fftfreq(len(xs), d=np.mean(np.diff(xs))))
    ks_y = np.sort(np.fft.fftfreq(len(ys), d=np.mean(np.diff(ys))))
    extent_ks = [ks_x.min(), ks_x.max(), ks_y.min(), ks_y.max()]

    for i,kx in enumerate(ks_x):
        print i
        for j,ky in enumerate(ks_y):
            # print str(kx) + ', ' + str(ky) + ':'
            for l,x in enumerate(xs):
                # ftZ[i,j] += Z[i,j]*np.sum(np.exp(-1j*2*np.pi*(kx*x + ky*ys)))
                for m,y in enumerate(ys):
                    ftZ[i,j] += Z[l,m]*np.exp(-1j*2*np.pi*(kx*x + ky*y))
            # ftZ[i,j] = Z[i,j]*np.sum(np.exp(-1j*2*np.pi*(kx*xs+ky*ys)))

    npftZ = np.fft.fft2(Z)

    fmt = 'o-'
    imgs = []

    ax = subplot(gs[0])
    ogim = ax.imshow(Z, interpolation='nearest', origin='center', extent=extent)
    ax.set_title('original')
    ax.set_xlabel('x', size=16)
    ax.set_ylabel('y', size=16)

    myax = subplot(gs[1])
    myim = myax.imshow(np.abs(ftZ), interpolation='nearest', origin='center', extent=extent_ks)
    myax.set_title('my fft')

    npax = subplot(gs[2])
    npim = npax.imshow(np.abs(np.fft.fftshift(npftZ)), interpolation='nearest', origin='center', extent=extent_ks)
    npax.set_title('numpy fft')

    dax = subplot(gs[3])
    dim = dax.imshow(np.abs(ftZ)-np.abs(np.fft.fftshift(npftZ)), interpolation='nearest', origin='center', extent=extent_ks)
    dax.set_title('my fft - numpy fft')

    for a in fig.axes[1:]:
        a.set_xlabel('k_x', size=16)
        a.set_ylabel('k_y', size=16)

    imgs = [ogim, myim, npim, dim]

    for i,axe in enumerate(fig.axes):
        divider = make_axes_locatable(axe)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        colorbar(imgs[i], cax=cax)

    # fig.colorbar(ogim, ax=ax)
    # fig.colorbar(npim, ax=npax)
    # fig.colorbar(myim, ax=myax)


gs.tight_layout(fig)

show()
