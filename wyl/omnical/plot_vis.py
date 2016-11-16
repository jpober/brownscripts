import matplotlib.pyplot as plt
import numpy as np

def plot_vis(dd,flag):
    for ii in range(0,dd.shape[0]):
        for jj in range(0,dd.shape[1]):
            if flag[ii][jj]: dd[ii][jj] = np.nan
    fig=plt.figure()
    p1=fig.add_subplot(1,2,1)
    p1.set_title("real")
    i1=p1.imshow(dd.real, interpolation='nearest',aspect='auto')
    fig.colorbar(i1)
    p2=fig.add_subplot(1,2,2)
    p2.set_title("imag")
    i2=p2.imshow(dd.imag, interpolation='nearest',aspect='auto')
    fig.colorbar(i2)
    plt.show()

def plot_vis2(dd,flag):
    for ii in range(0,dd.shape[0]):
        for jj in range(0,dd.shape[1]):
            if flag[ii][jj]: dd[ii][jj] = np.nan
    fig=plt.figure()
    p1=fig.add_subplot(1,2,1)
    p1.set_title("amplitude")
    i1=p1.imshow(np.abs(dd), interpolation='nearest',aspect='auto')
    fig.colorbar(i1)
    p2=fig.add_subplot(1,2,2)
    p2.set_title("phase")
    i2=p2.imshow(np.angle(dd), interpolation='nearest',aspect='auto')
    fig.colorbar(i2)
    plt.show()

def plot_tave_amp(dd,flag):
    for ii in range(0,dd.shape[0]):
        for jj in range(0,dd.shape[1]):
            if flag[ii][jj]: dd[ii][jj] = np.nan
    plt.scatter(np.abs(np.nanmean(dd,axis=0)),marker='+')
    plt.show()

def plot_tave_phs(dd,flag):
    for ii in range(0,dd.shape[0]):
        for jj in range(0,dd.shape[1]):
            if flag[ii][jj]: dd[ii][jj] = np.nan
    plt.scatter(np.angle(np.nanmean(dd,axis=0)),marker='+')
    plt.show()

