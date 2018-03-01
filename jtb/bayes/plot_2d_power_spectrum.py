import numpy as np
import matplotlib.gridspec as gridspec
import optparse, sys, os

from matplotlib.pyplot import *
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Option parser
o = optparse.OptionParser()
o.add_option('--uvdata',
    type = str,
    help = 'Filename for input uvw binary .npy file.')
o.add_option('--l_offset',
    type = float,
    default = 0.0,
    help = 'Moves source in l-direction by l_offset*ls.max().  Must be between 0 and 1.')
o.add_option('--m_offset',
    type = float,
    default = 0.0,
    help = 'Moves source in m-direction by m_offset*ms.max().  Must be between 0 and 1.')
opts,args = o.parse_args(sys.argv[1:])

pixel_nums = [25,50,100,150]
msizes = [13,11,9,7]
j = 0

for pixel_num in pixel_nums:
    ## ----------------- Construct sky ----------------- ##
    ls = np.cos(np.linspace(-np.pi, np.pi, pixel_num))
    ms = np.sin(np.linspace(-np.pi, np.pi, pixel_num))
    N_im = ls.size*ms.size
    extent_lm = [ls.min(), ls.max(), ms.min(), ms.max()]

    L, M = np.meshgrid(ls, ms)
    I = 0*L
    if len(ls) % 2 == 0:
        mid_l = len(ls)/2
    else:
        mid_l = int(len(ls)/2. - 0.5)
    if len(ms) % 2 == 0:
        mid_m = len(ms)/2
    else:
        mid_m = int(len(ms)/2. - 0.5)

    l_off = int(mid_l*opts.l_offset)
    m_off = int(mid_m*opts.m_offset)

    I[mid_m + m_off, mid_l + l_off] = 1.0
    N_freq = 1
    pixel_area = (ls[1]-ls[0])*(ms[1]-ms[0])

    # Construct position vectors for N_im pixels
    I_vec = np.reshape(I, (N_im*N_freq,), order='F')
    ls_vec, ms_vec = np.zeros(0), np.zeros(0)
    for l in ls:
        for m in ms:
            ls_vec = np.append(ls_vec, l)
            ms_vec = np.append(ms_vec, m)



    ## ----------------- Construct uv ----------------- ##

    # Read in u,v sampling points
    uvws = np.load(opts.uvdata)
    us = np.concatenate((-uvws[0,0,:], uvws[0,0,:]))
    vs = np.concatenate((-uvws[0,1,:], uvws[0,1,:]))
    N_vis = us.size*vs.size
    extent_uv = [us.min(), us.max(), vs.min(), vs.max()]


    ## ----------------- Construct nudft ----------------- ##

    # Construct DFT matrix using outer product
    DFT = np.exp(-1j*2*np.pi*(np.outer(us,ls_vec) + np.outer(vs, ms_vec)))

    # Compute visibilities from DFT
    Vs = np.dot(DFT, I_vec) #*pixel_area

    # Average over redundant baselines
    Vs_unique = np.zeros(0)
    uvs_unique, loc_inds, inverse_inds = np.unique(np.stack((us, vs), axis=1),
                                                                          axis = 0,
                                                                          return_index = True,
                                                                          return_inverse = True)
    for u,v in uvs_unique:
        avg_inds = np.where(np.logical_and(us == u, vs == v))[0]
        Vs_unique = np.append(Vs_unique, np.mean(Vs[avg_inds]))


    ## ----------------- Compute Power Spectrum ----------------- ##
    # P(k_u,k_v) = sum_r(Visibilities(u,v)**2) s.t. sqrt(u**2+v**2) = r

    rs = np.sqrt(uvs_unique[:,0]**2 + uvs_unique[:,1]**2)
    rs_unique = np.unique(rs)
    power_spec = np.zeros_like(rs_unique)*1.0j

    for i,r in enumerate(rs_unique):
        inds = np.where(rs == r)[0]
        power_spec[i] = np.mean(np.abs(Vs_unique[inds])**2)

    # plot power spectrum
    plot(rs_unique, power_spec, 'o', markersize = msizes[j], label = 'Numerical, %s'%(str(pixel_num)))
    j += 1


## ----------------- Analytic solution comparison ----------------- ##

# Point source, Flat beam
Vs_func = lambda u,v: np.exp(-2*np.pi*1j*(u*ls[mid_l + l_off] + v*ms[mid_m + m_off]))
Vs_analytic = Vs_func(uvs_unique[:,0], uvs_unique[:,1]) #*pixel_area
power_spec_analytic = np.zeros_like(rs_unique)*1.0j
for i,r in enumerate(rs_unique):
    inds = np.where(rs == r)[0]
    power_spec_analytic[i] = np.mean(np.abs(Vs_analytic[inds])**2)

plot(rs_unique, power_spec_analytic, 'o', label = 'Analytic', markersize = 5)


## ----------------- Plotting ----------------- ##

xlabel('\"k\"', size = 16)
ylabel('Power', size = 16)
legend(loc='upper left', fontsize = 16, ncol=2, bbox_to_anchor=(1.,1.))

tight_layout()

show()
