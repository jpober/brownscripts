#! /usr/bin/env python
'''After omnical calibration antennas 50 and 19 were found to be bad.
   Omnical flagged the data instead removing the antennas.
   This script removes those two antennas from the data set.'''
import aipy as a, numpy as n, os, sys, glob

def grid_jd(jds, temps, binsize=120):
    jdbin = binsize * a.ephem.second
    nbins = int((jds[-1] - jds[0]) / jdbin)
    wgts,bins = n.histogram(jds, bins=nbins)
    dats,bins = n.histogram(jds, weights=temps, bins=nbins)
    return dats / wgts, bins

def conv_tempfilename(f):
    f = os.path.basename(f)
    ftime = '%s/%s/%s %s:%s' % (f[:4],f[4:6],f[6:8],f[9:11],f[11:13])
    return a.phs.ephem2juldate(a.ephem.date(ftime))

def filter_files(filenames, start_jd, end_jd):
    jds = [conv_tempfilename(f) for f in filenames]
    return [f for f,jd in zip(filenames,jds) if jd >= start_jd and jd <= end_jd]

def get_temps(tempdir, file_jd, uv_filetime=600, temp_filetime=3600,
        t_offset=-2*a.ephem.hour):
    print '    Folding in temp data from', tempdir
    RECVR,ANT58,GOM_E = 0,1,2
    files = glob.glob(tempdir + '/201*.txt'); files.sort()
    start_jd = file_jd - a.ephem.second * (uv_filetime + temp_filetime) - t_offset
    end_jd = file_jd + a.ephem.second * (uv_filetime + temp_filetime) - t_offset
    files = filter_files(files, start_jd, end_jd)
    lines = [L.split() for f in files for L in open(f).readlines()]
    jds = n.array([a.phs.ephem2juldate(a.ephem.date(' '.join(L[:2]))) for L in lines])
    jds += t_offset # Adjust for local time being 2 hrs ahead of UTC
    dat = n.array([map(float,L[2:]) for L in lines])
    T_r, bins = grid_jd(jds, dat[:,RECVR])
    T_c, bins = grid_jd(jds, dat[:,ANT58])
    T_b, bins = grid_jd(jds, dat[:,ANT58])
    T_l, bins = grid_jd(jds, dat[:,GOM_E])
    return bins, (T_r,T_c,T_b,T_l)
    

import optparse
o = optparse.OptionParser()
#o.add_option('-t', '--tempdir', default='/data3/paper/psa/psalive',
o.add_option('-t', '--tempdir', 
    help='Directory containing temperature data from the (labjack) gainometer.')
opts,args = o.parse_args(sys.argv[1:])

aa = a.phs.ArrayLocation(('-30:43:17.5', '21:25:41.9')) # Karoo, ZAR, GPS. #elevation=1085m

#rewire = { 0:49,1:10,2:9,3:22,4:29,5:24,6:17,7:5,
#           8:47,9:25,10:1,11:35,12:34,13:27,14:56,15:30,
#           16:41,17:3,18:58,19:61,20:28,21:55,22:13,23:32,
#           24:19,25:48,26:4,27:18,28:51,29:57,30:59,31:23}
polswapant = []
polswap = {}
rmant = [19,50]

for filename in sys.argv[1:]:
    print filename, '->', filename+'c'
    if os.path.exists(filename+'c'):
        print '    File exists... skipping.'
        continue
    file_jd = float('.'.join(filename.split('.')[1:-1]))
    if opts.tempdir != None:
        bins, (T_r,T_c,T_b,T_l) = get_temps(opts.tempdir, file_jd)
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(filename+'c', status='new')
    curtime = 0
    def _mfunc(uv, p, d, f):
        global curtime
        crd,t,(i,j) = p
        p1,p2 = a.miriad.pol2str[uv['pol']]
        if i in rmant or j in rmant: return p, None, None
        if i in polswapant: p1 = polswap[p1]
        if j in polswapant: p2 = polswap[p2]
        # prevent multiple entries arising from xy and yx on autocorrelations
        if i == j and (p1,p2) == ('y','x'): return p, None, None

        #i,j = pol2to1[i][p1], pol2to1[j][p2]
        #i,j = rewire[i],rewire[j]
        if i > j: i,j,d = j,i,n.conjugate(d)
        #if i == 8: d = -d  # I think the dipole for this antenna is rotated 180 deg
        #if j == 8: d = -d
        uvo['pol'] = a.miriad.str2pol[p1+p2]

        if t != curtime:
            aa.set_jultime(t)
            uvo['lst'] = uvo['ra'] = uvo['obsra'] = aa.sidereal_time()
            curtime = t
        
        return (crd,t,(i,j)),d,f

    override = {
        'lst': aa.sidereal_time(),
        'ra': aa.sidereal_time(),
        'obsra': aa.sidereal_time(),
        #'sdf': sdf,
        'latitud': aa.lat,
        'dec': aa.lat,
        'obsdec': aa.lat,
        'longitu': aa.long,
        #'sfreq': sfreq,
        #'freq': sfreq,
        #'inttime': 5.37,
        #'nchan': 1024,
        #'nants': 64,
        #'ngains': 128,
        #'nspect0': 64,
        'pol':a.miriad.str2pol['xx'],
        #'telescop':'PAPER',
        #'bandpass': n.ones((33,1024), dtype=n.complex64).flatten(),
        'antpos': n.transpose(n.array([
            [0., 0., 0.], #0
            [0., 0., 0.], #1 
            [0., 0., 0.], #2
            [0., 0., 0.], #3
            [0., 0., 0.], #4
            [0., 0., 0.], #5
            [0., 0., 0.], #6
            [0., 0., 0.], #7
            [0., 0., 0.], #8
            [0., 0., 0.], #9
            [0., 0., 0.], #10
            [0., 0., 0.], #11
            [0., 0., 0.], #12
            [0., 0., 0.], #13
            [0., 0., 0.], #14
            [0., 0., 0.], #15
            [0., 0., 0.], #16
            [0., 0., 0.], #17
            [0., 0., 0.], #18
            [0., 0., 0.], #19
            [0., 0., 0.], #20
            [0., 0., 0.], #21
            [0., 0., 0.], #22
            [0., 0., 0.], #23
            [0., 0., 0.], #24
            [0., 0., 0.], #25
            [0., 0., 0.], #26
            [0., 0., 0.], #27
            [0., 0., 0.], #28
            [0., 0., 0.], #29
            [0., 0., 0.], #30
            [0., 0., 0.], #31
            [0., 0., 0.], #32
            [0., 0., 0.], #33
            [0., 0., 0.], #34
            [0., 0., 0.], #35
            [0., 0., 0.], #36
            [0., 0., 0.], #37
            [0., 0., 0.], #38
            [0., 0., 0.], #39
            [0., 0., 0.], #40
            [0., 0., 0.], #41
            [0., 0., 0.], #42
            [0., 0., 0.], #43
            [0., 0., 0.], #44
            [0., 0., 0.], #45
            [0., 0., 0.], #46
            [0., 0., 0.], #47
            [0., 0., 0.], #48
            [0., 0., 0.], #49
            [0., 0., 0.], #50
            [0., 0., 0.], #51
            [0., 0., 0.], #52
            [0., 0., 0.], #53
            [0., 0., 0.], #54
            [0., 0., 0.], #55
            [0., 0., 0.], #56
            [0., 0., 0.], #57
            [0., 0., 0.], #58
            [0., 0., 0.], #59
            [0., 0., 0.], #60
            [0., 0., 0.], #61
            [0., 0., 0.], #62
            [0., 0., 0.], #63
        ])).flatten(),
    }

    uvo.init_from_uv(uvi, override=override)

    # Add GoM information
    uvo.add_var('t_recvr', 'r')
    uvo.add_var('t_cable', 'r')
    uvo.add_var('t_balun', 'r')
    uvo.add_var('t_load', 'r')
    if opts.tempdir != None:
        curtime = None
        def mfunc(uv, p, d, f):
            global curtime, uvo
            (crd, t, (i,j)) = p
            if curtime != t:
                curtime = t
                b = int(n.floor((t - bins[0]) / (bins[1] - bins[0])))
                if not n.isnan(T_r[b]): uvo['t_recvr'] = T_r[b]
                if not n.isnan(T_c[b]): uvo['t_cable'] = T_c[b]
                if not n.isnan(T_b[b]): uvo['t_balun'] = T_b[b]
                if not n.isnan(T_l[b]): uvo['t_load']  = T_l[b]
            return _mfunc(uv, p, d, f)
    else: mfunc = _mfunc

    uvo.pipe(uvi, mfunc=mfunc, raw=True,
        append2hist='\nCORRECT: '+' '.join(sys.argv)+'\n')
    del(uvo)
