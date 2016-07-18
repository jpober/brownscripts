"""
Calfile for data taken from 2456240. - 2456378.
"""
import aipy as a, numpy as n


class AntennaArray(a.pol.AntennaArray):
    def __init__(self, *args, **kwargs):
        a.pol.AntennaArray.__init__(self, *args, **kwargs)
        self.ant_layout = kwargs.pop('ant_layout')
    def update(self):
        a.pol.AntennaArray.update(self)
    def get_params(self, ant_prms={'*':'*'}):
        try: prms = a.pol.AntennaArray.get_params(self, ant_prms)
        except(IndexError): return {}
        for k in ant_prms:
            if k == 'aa':
                if not prms.has_key('aa'): prms['aa'] = {}
                for val in ant_prms[k]:
                    if   val == 'tau_ns': prms['aa']['tau_ns'] = self.tau_ns
                    elif val == 'tau_ew': prms['aa']['tau_ew'] = self.tau_ew
                    elif val == 'gain': prms['aa']['gain'] = self.gain
            else:
                try: top_pos = n.dot(self._eq2zen, self[int(k)].pos)
                # XXX should multiply this by len_ns to match set_params.
                except(ValueError): continue
                if ant_prms[k] == '*':
                    prms[k].update({'top_x':top_pos[0], 'top_y':top_pos[1], 'top_z':top_pos[2]})
                else:
                    for val in ant_prms[k]:
                        if   val == 'top_x': prms[k]['top_x'] = top_pos[0]
                        elif val == 'top_y': prms[k]['top_y'] = top_pos[1]
                        elif val == 'top_z': prms[k]['top_z'] = top_pos[2]
        return prms
    def set_params(self, prms):
        changed = a.pol.AntennaArray.set_params(self, prms)
        for i, ant in enumerate(self):
            ant_changed = False
            top_pos = n.dot(self._eq2zen, ant.pos)
            try:
                top_pos[0] = prms[str(i)]['top_x']
                ant_changed = True
            except(KeyError): pass
            try:
                top_pos[1] = prms[str(i)]['top_y']
                ant_changed = True
            except(KeyError): pass
            try:
                top_pos[2] = prms[str(i)]['top_z']
                ant_changed = True
            except(KeyError): pass
            if ant_changed: ant.pos = n.dot(n.linalg.inv(self._eq2zen), top_pos) / a.const.len_ns * 100.
            changed |= ant_changed
        if changed: self.update()
        return changed



prms = {

### Careful -- roughly calculated from coordinates. Low-precision -- 116.67085 (long), -26.70331 (lat)
    'loc': ('-26:42:11.916', '116:40:15.06'), # MWA,  (GPS)
    'antpos': {
        
        #EastHex
        0:  {'top_x':  86.352,  'top_y':  287.221067,  'top_z':    375.032},
        1:  {'top_x': 100.352,  'top_y':  287.221067,  'top_z':    375.032},
        2:  {'top_x': 114.352,  'top_y':  287.221067,  'top_z':    375.032},
        3:  {'top_x': 128.352,  'top_y':  287.221067,  'top_z':    375.032},
        
        4:  {'top_x':  79.352,  'top_y':  275.0967113,  'top_z':    375.032},
        5:  {'top_x':  93.352,  'top_y':  275.0967113,  'top_z':    375.032},
        6:  {'top_x': 107.352,  'top_y':  275.0967113,  'top_z':    375.032},
        7:  {'top_x': 121.352,  'top_y':  275.0967113,  'top_z':    375.032},
        8:  {'top_x': 135.352,  'top_y':  275.0967113,  'top_z':    375.032},
        
        9:  {'top_x':  72.352,  'top_y':  262.9723557,  'top_z':    375.032},
        10: {'top_x':  86.352,  'top_y':  262.9723557,  'top_z':    375.032},
        11: {'top_x': 100.352,  'top_y':  262.9723557,  'top_z':    375.032},
        12: {'top_x': 114.352,  'top_y':  262.9723557,  'top_z':    375.032},
        13: {'top_x': 128.352,  'top_y':  262.9723557,  'top_z':    375.032},
        14: {'top_x': 142.352,  'top_y':  262.9723557,  'top_z':    375.032},
        
        15: {'top_x':  65.352,  'top_y': 250.848,  'top_z':    375.032},
        16: {'top_x':  79.352,  'top_y': 250.848,  'top_z':    375.032},
        17: {'top_x':  93.352,  'top_y': 250.848,  'top_z':    375.032},
        18: {'top_x': 121.352,  'top_y': 250.848,  'top_z':    375.032},
        19: {'top_x': 135.352,  'top_y': 250.848,  'top_z':    375.032},
        20: {'top_x': 149.352,  'top_y': 250.848,  'top_z':    375.032},
        
        21: {'top_x':  72.352,  'top_y': 238.7236444,  'top_z':    375.032},
        22: {'top_x':  86.352,  'top_y': 238.7236444,  'top_z':    375.032},
        23: {'top_x': 100.352,  'top_y': 238.7236444,  'top_z':    375.032},
        24: {'top_x': 114.352,  'top_y': 238.7236444,  'top_z':    375.032},
        25: {'top_x': 128.352,  'top_y': 238.7236444,  'top_z':    375.032},
        26: {'top_x': 142.352,  'top_y': 238.7236444,  'top_z':    375.032},
        
        27: {'top_x':  79.352,  'top_y': 226.5992887,  'top_z':    375.032},
        28: {'top_x':  93.352,  'top_y': 226.5992887,  'top_z':    375.032},
        29: {'top_x': 107.352,  'top_y': 226.5992887,  'top_z':    375.032},
        30: {'top_x': 121.352,  'top_y': 226.5992887,  'top_z':    375.032},
        31: {'top_x': 135.352,  'top_y': 226.5992887,  'top_z':    375.032},
        
        32: {'top_x':  86.352,  'top_y': 214.474933,  'top_z':    375.032},
        33: {'top_x': 100.352,  'top_y': 214.474933,  'top_z':    375.032},
        34: {'top_x': 114.352,  'top_y': 214.474933,  'top_z':    375.032},
        35: {'top_x': 128.352,  'top_y': 214.474933,  'top_z':    375.032},
        
        #south hex
        
        36: {'top_x':   -3.49,  'top_y': 156.590067,  'top_z':    376.351},
        37: {'top_x':   10.51,  'top_y': 156.590067,  'top_z':    376.351},
        38: {'top_x':   24.51,  'top_y': 156.590067,  'top_z':    376.351},
        39: {'top_x':   38.51,  'top_y': 156.590067,  'top_z':    376.351},

        40: {'top_x':  -10.49,  'top_y': 144.4657113,  'top_z':    376.351},
        41: {'top_x':    3.51,  'top_y': 144.4657113,  'top_z':    376.351},
        42: {'top_x':   17.51,  'top_y': 144.4657113,  'top_z':    376.351},
        43: {'top_x':   31.51,  'top_y': 144.4657113,  'top_z':    376.351},
        44: {'top_x':   45.51,  'top_y': 144.4657113,  'top_z':    376.351},
            
        45: {'top_x':  -17.49,  'top_y': 132.3413557,  'top_z':    376.351},
        46: {'top_x':   -3.49,  'top_y': 132.3413557,  'top_z':    376.351},
        47: {'top_x':   10.51,  'top_y': 132.3413557,  'top_z':    376.351},
        48: {'top_x':   24.51,  'top_y': 132.3413557,  'top_z':    376.351},
        49: {'top_x':   38.51,  'top_y': 132.3413557,  'top_z':    376.351},
        50: {'top_x':   52.51,  'top_y': 132.3413557,  'top_z':    376.351},
            
        51: {'top_x':  -24.49,  'top_y': 120.217,  'top_z':    376.351},
        52: {'top_x':  -10.49,  'top_y': 120.217,  'top_z':    376.351},
        53: {'top_x':    3.51,  'top_y': 120.217,  'top_z':    376.351},
        54: {'top_x':   31.51,  'top_y': 120.217,  'top_z':    376.351},
        55: {'top_x':   45.51,  'top_y': 120.217,  'top_z':    376.351},
        56: {'top_x':   59.51,  'top_y': 120.217,  'top_z':    376.351},
            
        57: {'top_x':  -17.49,  'top_y': 108.0926444,  'top_z':    376.351},
        58: {'top_x':   -3.49,  'top_y': 108.0926444,  'top_z':    376.351},
        59: {'top_x':   10.51,  'top_y': 108.0926444,  'top_z':    376.351},
        60: {'top_x':   24.51,  'top_y': 108.0926444,  'top_z':    376.351},
        61: {'top_x':   38.51,  'top_y': 108.0926444,  'top_z':    376.351},
        62: {'top_x':   52.51,  'top_y': 108.0926444,  'top_z':    376.351},
            
        63: {'top_x':  -10.49,  'top_y': 95.96828869,  'top_z':    376.351},
        64: {'top_x':    3.51,  'top_y': 95.96828869,  'top_z':    376.351},
        65: {'top_x':   17.51,  'top_y': 95.96828869,  'top_z':    376.351},
        66: {'top_x':   31.51,  'top_y': 95.96828869,  'top_z':    376.351},
        67: {'top_x':   45.51,  'top_y': 95.96828869,  'top_z':    376.351},
            
        68: {'top_x':   -3.49,  'top_y': 83.84393304,  'top_z':    376.351},
        69: {'top_x':   10.51,  'top_y': 83.84393304,  'top_z':    376.351},
        70: {'top_x':   24.51,  'top_y': 83.84393304,  'top_z':    376.351},
        71: {'top_x':   38.51,  'top_y': 83.84393304,  'top_z':    376.351},
    
    
        72: {'top_x':   -149.785,  'top_y': 265.814,  'top_z':    377.011},
        73: {'top_x':   -95.365,   'top_y': 270.188,  'top_z':    377.285},
        74: {'top_x':   -88.512,   'top_y': 266.021,  'top_z':    377.309},
        75: {'top_x':   -78.698,   'top_y': 258.431,  'top_z':    377.291},
        76: {'top_x':   -86.940,   'top_y': 248.452,  'top_z':    377.339},
        77: {'top_x':   -102.325,  'top_y': 244.088,  'top_z':    377.258},
        78: {'top_x':   -93.346,   'top_y': 242.279,  'top_z':    377.270},
        79: {'top_x':   -148.33,   'top_y': 220.271,  'top_z':    377.214},
        80: {'top_x':   -132.939,  'top_y': 282.611,  'top_z':    377.002},
        81: {'top_x':   -139.867,  'top_y': 275.498,  'top_z':    377.018},
        82: {'top_x':   -123.667,  'top_y': 284.854,  'top_z':    377.050},
        83: {'top_x':   -92.889,   'top_y': 281.358,  'top_z':    377.254},
        84: {'top_x':   -77.835,   'top_y': 282.839,  'top_z':    377.253},
        85: {'top_x':   -83.476,   'top_y': 277.563,  'top_z':    377.300},
        86: {'top_x':   -76.988,   'top_y': 272.637,  'top_z':    377.298},
        87: {'top_x':   -176.006,  'top_y': 289.903,  'top_z':    376.623},
        88: {'top_x':   -23.620,   'top_y': 201.546,  'top_z':    376.805},
        89: {'top_x':   -23.802,   'top_y': 222.353,  'top_z':    376.872},
        90: {'top_x':   -40.489,   'top_y': 235.170,  'top_z':    376.961},
        91: {'top_x':   -54.527,   'top_y': 237.518,  'top_z':    377.123},
        92: {'top_x':   -64.851,   'top_y': 250.813,  'top_z':    377.239},
        93: {'top_x':   -70.075,   'top_y': 259.402,  'top_z':    377.266},
        94: {'top_x':   -60.057,   'top_y': 263.676,  'top_z':    377.180},
        95: {'top_x':   -50.676,   'top_y': 264.546,  'top_z':    377.100},
        96: {'top_x':   -67.990,   'top_y': 271.103,  'top_z':    377.236},
        97: {'top_x':   -60.058,   'top_y': 271.905,  'top_z':    377.204},
        98: {'top_x':   -70.885,   'top_y': 278.818,  'top_z':    377.252},
        99: {'top_x':   -45.540,   'top_y': 273.245,  'top_z':    377.066},
       100: {'top_x':   -53.156,   'top_y': 276.268,  'top_z':    377.090},
       101: {'top_x':   -22.659,   'top_y': 263.694,  'top_z':    376.948},
       102: {'top_x':   -55.097,   'top_y': 284.285,  'top_z':    377.181},
       103: {'top_x':   -14.900,   'top_y': 273.870,  'top_z':    376.910},
       104: {'top_x':   -105.281,  'top_y': 217.180,  'top_z':    377.255},
       105: {'top_x':   -98.53,    'top_y': 230.162,  'top_z':    377.260},
       106: {'top_x':   -81.845,   'top_y': 229.337,  'top_z':    377.270},
       107: {'top_x':   -79.541,   'top_y': 238.128,  'top_z':    377.270},
       108: {'top_x':   -75.120,   'top_y': 247.002,  'top_z':    377.255},
       109: {'top_x':   -71.017,   'top_y': 235.929,  'top_z':    377.239},
       110: {'top_x':   -62.567,   'top_y': 228.987,  'top_z':    377.210},
       111: {'top_x':   -50.691,   'top_y': 221.885,  'top_z':    377.078},
       112: {'top_x':   -160.465,  'top_y': 573.207,  'top_z':    374.773},
       113: {'top_x':   -128.046,  'top_y': 350.494,  'top_z':    376.570},
       114: {'top_x':   -69.771,   'top_y': 294.031,  'top_z':    377.192},
       115: {'top_x':   -78.731,   'top_y': 297.259,  'top_z':    377.196},
       116: {'top_x':   -103.941,  'top_y': 300.914,  'top_z':    377.050},
       117: {'top_x':   -100.172,  'top_y': 288.896,  'top_z':    377.157},
       118: {'top_x':   -395.452,  'top_y': 271.152,  'top_z':    374.515},
       119: {'top_x':   -263.547,  'top_y': 389.274,  'top_z':    375.322},
       120: {'top_x':    326.141,  'top_y': 203.608,  'top_z':    373.945},
       121: {'top_x':    173.425,  'top_y': 193.595,  'top_z':    375.032},
       122: {'top_x':   -2.515,    'top_y': 292.841,  'top_z':    376.892},
       123: {'top_x':   -46.554,   'top_y': 287.469,  'top_z':    377.115},
       124: {'top_x':   -56.141,   'top_y': 295.922,  'top_z':    377.156},
       125: {'top_x':   -33.705,   'top_y': 351.048,  'top_z':    376.976},
       126: {'top_x':    50.127,   'top_y': 536.159,  'top_z':    376.181},
       127: {'top_x':    84.880,   'top_y': 519.141,  'top_z':    376.314},
    #to be excluded
    #72: {'top_x':  0,  'top_y': 0,  'top_z':   0},
    },
    'ant_layout': n.array(
        [[  0,  1,  2,  3, -1, -1, -1, -1, -1, -1, -1, -1, -1, 36, 37, 38, 39, -1, -1, -1],
         [  4,  5,  6,  7,  8, -1, -1, -1, -1, -1, -1, -1, -1, 40, 41, 42, 43, 44, -1, -1],
         [  9, 10, 11, 12, 13, 14, -1, -1, -1, -1, -1, -1, -1, 45, 46, 47, 48, 49, 50, -1],
         [ 15, 16, 17, -1, 18, 19, 20, -1, -1, -1, -1, -1, -1, 51, 52, 53, -1, 54, 55, 56],
         [ -1, 21, 22, 23, 24, 25, 26, -1, -1, -1, -1, -1, -1, -1, 57, 58, 59, 60, 61, 62],
         [ -1, -1, 27, 28, 29, 30, 31, -1, -1, -1, -1, -1, -1, -1, -1, 63, 64, 65, 66, 67],
         [ -1, -1, -1, 32, 33, 34, 35, -1, -1, -1, -1, -1, -1, -1, -1, -1, 68, 69, 70, 71]]
         )
}
#

def get_aa(freqs):
    '''Return the AntennaArray to be used for simulation.'''
    location = prms['loc']
    antennas = []
    phsoff = {'x':[0.,0.], 'y':[0.,0.]}
    nants = len(prms['antpos'])
    for i in range(nants):
        bp_r = {'x':n.array([0.0]*9), 'y':n.array([0.0]*9)}
        bp_i = {'x':n.array([0.0]*3), 'y':n.array([0.0]*3)}
        antennas.append(a.pol.Antenna(0,0,0, a.fit.Beam(freqs),
             amp={'x':1.,'y':1.},bp_r=bp_r, bp_i=bp_i, phsoff=phsoff))
    aa = AntennaArray(prms['loc'], antennas, ant_layout=prms['ant_layout'])
#    for i in range(nants):
#        pos = prms['antpos-top'][i]
#        i = str(i)
#        aa.set_params({i:{'top_x':pos[0], 'top_y':pos[1], 'top_z':pos[2]}})
    for i in range(nants):
        pos = prms['antpos'][i]
        i = str(i)
        aa.set_params({i:pos})
    return aa

