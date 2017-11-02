#! /usr/bin/env python
"""
Filter visibilites per-baseline using a delay transform.
"""

import aipy as a, numpy as n, os, sys, optparse 
import capo as C
import pyuvdata

o = optparse.OptionParser()
o.set_usage('pspec_prep.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, cal=True, ant=True, pol=True)
o.add_option('--nophs', dest='nophs', action='store_true',
    help='Do not phase to zenith bin.')
o.add_option('--nogain', dest='nogain', action='store_true',
    help='Do not normalize gain.')
o.add_option('--window', dest='window', default='none',
    help='DSP window to use.  Default: none')
o.add_option('--horizon', dest='horizon', type=float, default=0.,
    help='An additional additive term (in ns) applied to the baseline length to determine horizon cutoff.  Default is 0.')
o.add_option('--clean', dest='clean', type='float', default=1e-5,
    help='Deconvolve delay-domain data by the response that results from flagged data.  Specify a tolerance for termination.  Default 1e-5')
o.add_option('--model', dest='model', action='store_true',
    help='Return the foreground model summed with the residuals (in Fourier space).')
opts, args = o.parse_args(sys.argv[1:])

def findEdges(data):
    ct1 = 0
    ct2 = len(data)-1
    while n.isnan(data[ct1]) == True or data[ct1] == 0:
        ct1+=1
        if ct1==len(data)-1:
            el = 0
            break
    el = ct1-1
    while n.isnan(data[ct2]) == True or data[ct2] == 0:
        ct2-=1
        if ct2==0:
            er = 202
            break
    er = ct2+1
    return el,er

#uv = a.miriad.UV(args[0])
#aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
#filters = gen_skypass_delay(aa, uv['sdf'], uv['nchan'], max_bl_add=opts.horizon)

def WBCLEAN(d,f,clean,ant1,ant2,WIN,rl,ll):
    clean = float(clean)
    window = a.dsp.gen_window(rl-ll,window=WIN)
    if n.average(f) < .5: return n.ones_like(f)
    d_res = n.zeros_like(d)
    dw = d
    fw = f
    dw[:,ll:rl] = d[:,ll:rl]*window
    fw[:,ll:rl] = f[:,ll:rl]*window
    uthresh,lthresh = filters[(ant1,ant2)]
    print lthresh,uthresh
    area = n.ones(dw[0,:].size, dtype=n.int)
    area[uthresh:lthresh] = 0
    for i in range(dw.shape[0]):
        _d = n.fft.ifft(dw[i,:])
        _w = n.fft.ifft(fw[i,:])
        _d_cl, info = a.deconv.clean(_d, _w, tol=clean, area=area, stop_if_div=False, maxiter=100)
        d_mdl = n.fft.fft(_d_cl)
        d_res[i,:] = d[i,:] - d_mdl*f[i,:]
    d_res[:,ll:rl] = d_res[:,ll:rl]/window
    return d_res

def get_filters(obs,hor_ext,chans):
    uva = a.miriad.UV(obs)
    aa = a.cal.get_aa('psa6240_FHD', uva['sdf'], uva['sfreq'], chans)
    filters = C.dspec.wedge_width_by_bl(aa, uva['sdf'], chans, offset=hor_ext)
    return filters

for uvfile in args:
    print uvfile
    
    uvd = pyuvdata.miriad.Miriad()
    uvd.read_miriad(uvfile)
    filters = get_filters(uvfile,15,203)
#    for b in n.unique(uvd.baseline_array):
    for b in filters.keys():
        #ant1,ant2 = uvd.baseline_to_antnums(b)
        ant1 = b[0]
        ant2 = b[1]
        idx = uvd.baseline_array == uvd.antnums_to_baseline(ant1,ant2)
        data = uvd.data_array[idx,0,:,0]
        leftLim,rightLim = findEdges(n.sum(data,0))
        if rightLim-leftLim <0:
            #print 'Continued...',ant1,ant2
            continue
        try:
            clnData = WBCLEAN(data,n.logical_not(uvd.flag_array[idx,0,:,0]).astype(float),opts.clean,ant1,ant2,opts.window,rightLim,leftLim)
            uvd.data_array[idx,0,:,0] = clnData
        except KeyError:
            print 'Key not found.',ant1,ant2
        except ValueError:
            print 'Wrong shape.'

    uvd.write_miriad(uvfile+'B')
    del(uvd)
