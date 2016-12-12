#!/gpfs/runtime/opt/python/2.7.3/bin/python

import aipy as a, numpy as n, pylab as p, sys, optparse

# Loop through UV files collecting relevant data
plot_x = {}
plot_t = {'jd':[], 'lst':[], 'cnt':[]}
times = []


for uvfile in args:
    print 'Reading', uvfile
    uv = a.miriad.UV(uvfile)
    else: aa = None
    # Only select data that is needed to plot
    # Read data from a single UV file
    for (uvw,t,(i,j)),d in uv.all():
        bl = '%d,%d,%d' % (i,j,uv['pol'])
        if len(times) == 0 or times[-1] != t:
            times.append(t)
            # Implement time selection
            use_this_time = time_sel(t, (len(times)-1) / opts.decimate)
            if use_this_time:
                if aa == None: lst = uv['lst']
                else:
                    aa.set_jultime(t)
                    lst = aa.sidereal_time()
                plot_t['lst'].append(lst)
                plot_t['jd'].append(t)
                plot_t['cnt'].append(len(times)-1)
        if not use_this_time: continue
        d = d.take(chans)
        #apply cal phases
        if not opts.cal is None:
            aa.set_jultime(t)
            if not opts.src is None:
                src.compute(aa)
                d = aa.phs2src(d, src, i, j)
            #else: took out this mode because it's not used, and prefer not to phase.
            #    d *= n.exp(-1j*n.pi*aa.get_phs_offset(i,j))
        # Do delay transform if required
        if opts.delay:
            w = a.dsp.gen_window(d.shape[-1], window=opts.window)
            if opts.unmask:
                flags = n.ones(d.shape, dtype=n.float)
                d = d.data
            else:
                flags = n.logical_not(d.mask).astype(n.float)
                d = d.filled(0)
            d = n.fft.ifft(d*w)
            ker = n.fft.ifft(flags*w)
            gain = a.img.beam_gain(ker)
            if not opts.clean is None and not n.all(d == 0):
                d, info = a.deconv.clean(d, ker, tol=opts.clean)
                d += info['res'] / gain
            d = n.ma.array(d)
            d = n.fft.fftshift(d, axes=0)
        elif opts.unmask: d = d.data
        d.shape = (1,) + d.shape
        if not plot_x.has_key(bl): plot_x[bl] = []
        plot_x[bl].append(d)
    del(uv)

