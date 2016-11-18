import matplotlib.pyplot as plt
import numpy as np, copy

def plot_vis(data,flag):
    dd=copy.copy(data)
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

def plot_vis2(data,flag):
    dd=copy.copy(data)
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

def plot_tave_amp(data,flag):
    dd=copy.copy(data)
    for ii in range(0,dd.shape[0]):
        for jj in range(0,dd.shape[1]):
            if flag[ii][jj]: dd[ii][jj] = np.nan
    plt.scatter(np.abs(np.nanmean(dd,axis=0)),marker='+')
    plt.show()

def plot_tave_phs(data,flag):
    dd=copy.copy(data)
    for ii in range(0,dd.shape[0]):
        for jj in range(0,dd.shape[1]):
            if flag[ii][jj]: dd[ii][jj] = np.nan
    plt.scatter(np.angle(np.nanmean(dd,axis=0)),marker='+')
    plt.show()

def plot_vis3(data,flag):
    dd=copy.copy(data)
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
    i2=p2.imshow(np.cos(np.angle(dd)), interpolation='nearest',aspect='auto')
    fig.colorbar(i2)
    plt.show()

def low_pass_filter(arr,tau=1e-7,nres=1024,fqstart=1.6715501e+08,fqend=1.9779501e+08):
    fq=np.linspace(fqstart,fqend,nres)
    tfq=np.fft.fftfreq(nres,fq[1]-fq[0])
    farr=np.fft.fft(arr,axis=1,n=nres)
    for ii in range(0,nres):
        if abs(tfq[ii])>tau: farr[:,ii]=0
    return np.fft.ifft(farr,axis=1,n=arr.shape[1])

def avoid_flag_filter(arr,mask,tau=1e-7,nres=1024,fqstart=1.6715501e+08,fqend=1.9779501e+08):
    ar2=copy.copy(arr)
    m=np.sum(1-mask,axis=0)
    m=(1-m.astype(bool)).astype(bool)
    for ii in range(0,m.size):
        if m[ii]:
            nn=0
            while(m[(ii+nn)%m.size] and m[ii-nn]):
                nn+=1
                if not m[ii-nn]: ar2[:,ii]=ar2[:,ii-nn]
                elif not m[(ii+nn)%m.size]: ar2[:,ii]=ar2[:,(ii+nn)%m.size]
    return low_pass_filter(ar2,tau=tau,nres=nres,fqstart=fqstart,fqend=fqend)
