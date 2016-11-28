import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np, copy

def plot_vis(data,flag,freq=None,title=''):
    dd=copy.copy(data)
    for ii in range(0,dd.shape[0]):
        for jj in range(0,dd.shape[1]):
            if flag[ii][jj]: dd[ii][jj] = np.nan
    if bool(type(freq)):
        x=[0,freq.size/3,freq.size*2/3,freq.size]
        lx=['%.2f'%(freq[0]/1e6),'%.2f'%(freq[freq.size/3]/1e6),'%.2f'%(freq[freq.size*2/3]/1e6),'%.2f'%(freq[-1]/1e6)]
    fig=plt.figure()
    p1=fig.add_subplot(1,2,1)
    p1.set_title(title+"real")
    i1=p1.imshow(dd.real, interpolation='nearest',aspect='auto')
    fig.colorbar(i1)
    if bool(type(freq)):
        p1.set_xticks(x)
        p1.set_xticklabels(lx)
        p1.set_xlabel('frequency(MHz)')
    p2=fig.add_subplot(1,2,2)
    p2.set_title(title+"imag")
    i2=p2.imshow(dd.imag, interpolation='nearest',aspect='auto')
    fig.colorbar(i2)
    if bool(type(freq)):
        p2.set_xticks(x)
        p2.set_xticklabels(lx)
        p2.set_xlabel('frequency(MHz)')
    plt.show()

def plot_vis2(data,flag,freq=None,title='',log_plot=False):
    dd=copy.copy(data)
    for ii in range(0,dd.shape[0]):
        for jj in range(0,dd.shape[1]):
            if flag[ii][jj]: dd[ii][jj] = np.nan
    if bool(type(freq)):
        x=[0,freq.size/3,freq.size*2/3,freq.size]
        lx=['%.2f'%(freq[0]/1e6),'%.2f'%(freq[freq.size/3]/1e6),'%.2f'%(freq[freq.size*2/3]/1e6),'%.2f'%(freq[-1]/1e6)]
    fig=plt.figure()
    p1=fig.add_subplot(1,2,1)
    p1.set_title(title+"amplitude")
    if log_plot:
        i1=p1.imshow(np.abs(dd), interpolation='nearest',aspect='auto',norm=colors.LogNorm())
    else:
        i1=p1.imshow(np.abs(dd), interpolation='nearest',aspect='auto')
    fig.colorbar(i1)
    if bool(type(freq)):
        p1.set_xticks(x)
        p1.set_xticklabels(lx)
        p1.set_xlabel('frequency(MHz)')
    p2=fig.add_subplot(1,2,2)
    p2.set_title(title+"phase")
    i2=p2.imshow(np.angle(dd), interpolation='nearest',aspect='auto')
    fig.colorbar(i2)
    if bool(type(freq)):
        p2.set_xticks(x)
        p2.set_xticklabels(lx)
        p2.set_xlabel('frequency(MHz)')
    plt.show()

def plot_vis_xy(data,flag,freq=None,title='',log_plot=False):
    dd=copy.copy(data)
    for ii in range(0,dd.shape[0]):
        for jj in range(0,dd.shape[1]):
            if flag[ii][jj]: dd[ii][jj] = np.nan
    if not freq is None:
        x=[0,freq.size/3,freq.size*2/3,freq.size]
        lx=['%.2f'%(freq[0]/1e6),'%.2f'%(freq[freq.size/3]/1e6),'%.2f'%(freq[freq.size*2/3]/1e6),'%.2f'%(freq[-1]/1e6)]
    fig=plt.figure()
    p1=fig.add_subplot(2,2,1)
    p1.set_title(title+"_xx_amplitude")
    if log_plot:
        i1=p1.imshow(np.abs(dd[:,:,0]), interpolation='nearest',aspect='auto',norm=colors.LogNorm())
    else:
        i1=p1.imshow(np.abs(dd[:,:,0]), interpolation='nearest',aspect='auto')
    fig.colorbar(i1)
    if not freq is None:
        p1.set_xticks(x)
        p1.set_xticklabels(lx)
    #    p1.set_xlabel('frequency(MHz)')
    p2=fig.add_subplot(2,2,2)
    p2.set_title(title+"_xx_phase")
    i2=p2.imshow(np.angle(dd[:,:,0]), interpolation='nearest',aspect='auto')
    fig.colorbar(i2)
    if not freq is None:
        p2.set_xticks(x)
        p2.set_xticklabels(lx)
    #    p2.set_xlabel('frequency(MHz)')
    p3=fig.add_subplot(2,2,3)
    p3.set_title(title+"_yy_amplitude")
    if log_plot:
        i3=p3.imshow(np.abs(dd[:,:,1]), interpolation='nearest',aspect='auto',norm=colors.LogNorm())   
    else:
        i3=p3.imshow(np.abs(dd[:,:,1]), interpolation='nearest',aspect='auto')
    fig.colorbar(i3)
    if not freq is None:
        p3.set_xticks(x)
        p3.set_xticklabels(lx)
        p3.set_xlabel('frequency(MHz)')
    p4=fig.add_subplot(2,2,4)
    p4.set_title(title+"_yy_phase")
    i4=p4.imshow(np.angle(dd[:,:,1]), interpolation='nearest',aspect='auto')
    fig.colorbar(i4)
    if not freq is None:
        p4.set_xticks(x)
        p4.set_xticklabels(lx)
        p4.set_xlabel('frequency(MHz)')
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
    if len(mask.shape)==2:
        m=np.sum(1-mask,axis=0)
        m=(1-m.astype(bool)).astype(bool)
    elif len(mask.shape)==1:
        m=mask
    for ii in range(0,m.size):
        if m[ii]:
            nn=0
            while(m[(ii+nn)%m.size] and m[ii-nn]):
                nn+=1
                if not m[ii-nn]: ar2[:,ii]=ar2[:,ii-nn]
                elif not m[(ii+nn)%m.size]: ar2[:,ii]=ar2[:,(ii+nn)%m.size]
    return low_pass_filter(ar2,tau=tau,nres=nres,fqstart=fqstart,fqend=fqend)
