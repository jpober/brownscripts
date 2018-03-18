#!/bin/env python

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
from matplotlib.colors import SymLogNorm
import numpy as np, copy

def plot_vis(data,flag=None,freq=None,title=''):
    if flag is None: flag=np.zeros(data.shape,dtype=bool)
    dd=copy.copy(data)
    for ii in range(0,dd.shape[0]):
        for jj in range(0,dd.shape[1]):
            if flag[ii][jj]: dd[ii][jj] = np.nan
    if not freq is None:
        x=[0,freq.size/3,freq.size*2/3,freq.size]
        lx=['%.2f'%(freq[0]/1e6),'%.2f'%(freq[freq.size/3]/1e6),'%.2f'%(freq[freq.size*2/3]/1e6),'%.2f'%(freq[-1]/1e6)]
    fig=plt.figure()
    p1=fig.add_subplot(1,2,1)
    p1.set_title(title+"real")
    i1=p1.imshow(dd.real, interpolation='nearest',aspect='auto')
    fig.colorbar(i1)
    if not freq is None:
        p1.set_xticks(x)
        p1.set_xticklabels(lx)
        p1.set_xlabel('frequency(MHz)')
    p2=fig.add_subplot(1,2,2)
    p2.set_title(title+"imag")
    i2=p2.imshow(dd.imag, interpolation='nearest',aspect='auto')
    fig.colorbar(i2)
    if not freq is None:
        p2.set_xticks(x)
        p2.set_xticklabels(lx)
        p2.set_xlabel('frequency(MHz)')
    plt.show()

def plot_vis2(data,flag=None,freq=None,title='',log_plot=False):
    if flag is None: flag=np.zeros(data.shape,dtype=bool)
    dd=copy.copy(data)
    for ii in range(0,dd.shape[0]):
        for jj in range(0,dd.shape[1]):
            if flag[ii][jj]: dd[ii][jj] = np.nan
    if not freq is None:
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
    if not freq is None:
        p1.set_xticks(x)
        p1.set_xticklabels(lx)
        p1.set_xlabel('frequency(MHz)')
    p2=fig.add_subplot(1,2,2)
    p2.set_title(title+"phase")
    i2=p2.imshow(np.angle(dd), interpolation='nearest',aspect='auto')
    fig.colorbar(i2)
    if not freq is None:
        p2.set_xticks(x)
        p2.set_xticklabels(lx)
        p2.set_xlabel('frequency(MHz)')
    plt.show()

def plot_chisq(data,flag=None,freq=np.linspace(167.075,197.715,384),title='',cut=1.2,log_plot=False):
    if flag is None: flag=np.zeros(data.shape,dtype=bool)
    dd=copy.copy(data)
    ind = np.where(dd>cut)
    dd[ind]=cut
    for ii in range(0,dd.shape[0]):
        for jj in range(0,dd.shape[1]):
            if flag[ii][jj]: dd[ii][jj] = np.nan
    if not freq is None:
        x=[0,freq.size/3,freq.size*2/3,freq.size-1]
        lx=['%.2f'%(freq[0]),'%.2f'%(freq[freq.size/3]),'%.2f'%(freq[freq.size*2/3]),'%.2f'%(freq[-1])]
    fig=plt.figure()
    p1=fig.add_subplot(1,1,1)
    p1.set_title("chi-square/DOF")
    if not freq is None:
        p1.set_xticks(x)
        p1.set_xticklabels(lx)
        p1.set_xlabel('Frequency(MHz)')
    p1.set_ylabel('Time steps')
    if log_plot:
        i1=p1.imshow(np.abs(dd), interpolation='nearest',aspect='auto',norm=colors.LogNorm())
    else:
        i1=p1.imshow(np.abs(dd), interpolation='nearest',aspect='auto')
    fig.colorbar(i1)
    plt.show()

def plot_phase(data,flag=None,freq=np.linspace(167.075,197.715,384),title='',log_plot=False):
    if flag is None: flag=np.zeros(data.shape,dtype=bool)
    dd=copy.copy(data)
#    ind = np.where(dd>cut)
#    dd[ind]=cut
    for ii in range(0,dd.shape[0]):
        for jj in range(0,dd.shape[1]):
            if flag[ii][jj]: dd[ii][jj] = np.nan
    if not freq is None:
        x=[0,freq.size/3,freq.size*2/3,freq.size-1]
        lx=['%.2f'%(freq[0]),'%.2f'%(freq[freq.size/3]),'%.2f'%(freq[freq.size*2/3]),'%.2f'%(freq[-1])]
    fig=plt.figure()
    p1=fig.add_subplot(1,1,1)
    p1.set_title(title)
    if not freq is None:
        p1.set_xticks(x)
        p1.set_xticklabels(lx)
        p1.set_xlabel('Frequency(MHz)')
    p1.set_ylabel('Time steps')
    if log_plot:
        i1=p1.imshow(np.angle(dd), interpolation='nearest',aspect='auto',norm=colors.LogNorm())
    else:
        i1=p1.imshow(np.angle(dd), interpolation='nearest',aspect='auto')
    fig.colorbar(i1)
    plt.show()

def plot_vis_xy(data,flag=None,freq=None,title='',log_plot=False):
    if flag is None: flag=np.zeros(data[:,:,0].shape,dtype=bool)
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


def plot_tave_amp(data,flag=None):
    if flag is None: flag=np.zeros(data.shape,dtype=bool)
    dd=copy.copy(data)
    for ii in range(0,dd.shape[0]):
        for jj in range(0,dd.shape[1]):
            if flag[ii][jj]: dd[ii][jj] = np.nan
    plt.scatter(np.abs(np.nanmean(dd,axis=0)),marker='+')
    plt.show()

def plot_tave_phs(data,flag=None):
    if flag is None: flag=np.zeros(data.shape,dtype=bool)
    dd=copy.copy(data)
    for ii in range(0,dd.shape[0]):
        for jj in range(0,dd.shape[1]):
            if flag[ii][jj]: dd[ii][jj] = np.nan
    plt.scatter(np.angle(np.nanmean(dd,axis=0)),marker='+')
    plt.show()

def plot_vis3(data,flag=None):
    if flag is None: flag=np.zeros(data.shape,dtype=bool)
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

def plot_vis_3(d1,d2,d3):
    fig = plt.figure()
    p1=fig.add_subplot(3,2,1)
    p1.set_title("data_amp")
    i1=p1.imshow(np.abs(d1), interpolation='nearest',aspect='auto')
    fig.colorbar(i1)
    p1.xaxis.set_ticklabels([])
    p2=fig.add_subplot(3,2,2)
    p2.set_title("data_phs")
    i2=p2.imshow(np.angle(d1),interpolation='nearest',aspect='auto')
    fig.colorbar(i2)
    p2.xaxis.set_ticklabels([])
    p3=fig.add_subplot(3,2,3)
    p3.set_title("model_amp")
    i3=p3.imshow(np.abs(d2),interpolation='nearest',aspect='auto')
    fig.colorbar(i3)
    p3.xaxis.set_ticklabels([])
    p4=fig.add_subplot(3,2,4)
    p4.set_title("model_phs")
    i4=p4.imshow(np.angle(d2),interpolation='nearest',aspect='auto')
    fig.colorbar(i4)
    p4.xaxis.set_ticklabels([])
    p5=fig.add_subplot(3,2,5)
    p5.set_title("sol*model_amp")
    i5=p5.imshow(np.abs(d3),interpolation='nearest',aspect='auto')
    fig.colorbar(i5)
    p6=fig.add_subplot(3,2,6)
    p6.set_title("sol*model_phs")
    i6=p6.imshow(np.angle(d3),interpolation='nearest',aspect='auto')
    fig.colorbar(i6)
    plt.show()

def plot_vis_ri(d1,d2,d3):
    rmax = np.nanmax(d2.real)
    rmin = np.nanmin(d2.real)
    imax = np.nanmax(d2.imag)
    imin = np.nanmin(d2.imag)
    fig = plt.figure()
    p1=fig.add_subplot(3,2,1)
    p1.set_title("d1 real")
    i1=p1.imshow(np.clip(d1.real,rmin,rmax), interpolation='nearest',aspect='auto')
    c1=fig.colorbar(i1)
    p1.xaxis.set_ticklabels([])
    p2=fig.add_subplot(3,2,2)
    p2.set_title("d1 imag")
    i2=p2.imshow(np.clip(d1.imag,imin,imax),interpolation='nearest',aspect='auto')
    fig.colorbar(i2)
    p2.xaxis.set_ticklabels([])
    p3=fig.add_subplot(3,2,3)
    p3.set_title("d2 real")
    i3=p3.imshow(d2.real,interpolation='nearest',aspect='auto')
    fig.colorbar(i3)
    p3.xaxis.set_ticklabels([])
    p4=fig.add_subplot(3,2,4)
    p4.set_title("d2 imag")
    i4=p4.imshow(d2.imag,interpolation='nearest',aspect='auto')
    fig.colorbar(i4)
    p4.xaxis.set_ticklabels([])
    p5=fig.add_subplot(3,2,5)
    p5.set_title("diff real")
    i5=p5.imshow(d3.real,interpolation='nearest',aspect='auto')
    fig.colorbar(i5)
    p6=fig.add_subplot(3,2,6)
    p6.set_title("diff imag")
    i6=p6.imshow(d3.imag,interpolation='nearest',aspect='auto')
    fig.colorbar(i6)
    plt.show()

def delayps(data,bl,tau,log_plot=False):
    if log_plot: plt.imshow(data, aspect="auto", extent=(bl[0],bl[-1], tau[1], tau[-1]), interpolation='nearest', norm=LogNorm(),cmap='jet')
    else: plt.imshow(data, aspect="auto", extent=(bl[0],bl[-1], tau[1], tau[-1]), interpolation='nearest',cmap='jet')
    plt.xscale('log')
    plt.yscale('log')
    plt.colorbar()
    plt.show()

def delaydiff(diff,bl,tau,linthresh=1.):
    plt.imshow(diff,aspect="auto", extent=(bl[0],bl[-1], tau[1], tau[-1]), interpolation='nearest', norm=colors.SymLogNorm(linthresh=linthresh),cmap='jet') 
    plt.xscale('log')
    plt.yscale('log')
    plt.colorbar()
    plt.show()   

def organize_power(sav_dict,power_key):
    kx = sav_dict['kx_mpc']
    ky = sav_dict['ky_mpc']
    kz = sav_dict['kz_mpc']
    power = sav_dict[power_key]
    nx,ny,nz = kx.size,ky.size,kz.size
    dat = np.zeros((nz,nx,nx),dtype=power.dtype)
    for zi in range(nz):
        arr1 = power[zi]
        arr2 = np.rot90(arr1,2).conj()
        for ii in range(nx/2):
            dat[zi][ii] = arr2[ii]
            dat[zi][ii+nx/2+1] = arr1[ii+1]
        for jj in range(nx/2):
            dat[zi][nx/2][jj] = arr2[nx/2][jj]
            dat[zi][nx/2][jj+nx/2+1] = arr1[0][jj+nx/2+1]
    return dat

def plot_power_zslice(dat,z,kx,clim=None):
    plt.imshow(dat[z],aspect='equal',interpolation='nearest',extent=(kx[0],kx[-1],kx[0],kx[-1]),cmap='jet',norm=LogNorm(),clim=clim)
    plt.colorbar()
    plt.show()

