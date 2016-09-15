import aipy as a, numpy as n, os

class AntennaArray(a.pol.AntennaArray):
    def __init__(self, *args, **kwargs):
        a.pol.AntennaArray.__init__(self, *args, **kwargs)
        self.array_params = {}
    def get_ant_params(self, ant_prms={'*':'*'}):
        prms = a.fit.AntennaArray.get_params(self, ant_prms)
        for k in ant_prms:
            top_pos = n.dot(self._eq2zen, self[int(k)].pos)
            if ant_prms[k] == '*':
                prms[k].update({'top_x':top_pos[0], 'top_y':top_pos[1], 'top_z':top_pos[2]})
            else:
                for val in ant_prms[k]:
                    if   val == 'top_x': prms[k]['top_x'] = top_pos[0]
                    elif val == 'top_y': prms[k]['top_y'] = top_pos[1]
                    elif val == 'top_z': prms[k]['top_z'] = top_pos[2]
        return prms
    def set_ant_params(self, prms):
        changed = a.fit.AntennaArray.set_params(self, prms)
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
            if ant_changed: ant.pos = n.dot(n.linalg.inv(self._eq2zen), top_pos)
            changed |= ant_changed
        return changed 
    def get_arr_params(self):
        return self.array_params
    def set_arr_params(self, prms):
        for param in prms:
            self.array_params[param] = prms[param]
            if param == 'dish_size_in_lambda':
                FWHM = 2.35*(0.45/prms[param]) #radians
                self.array_params['obs_duration'] = 60.*FWHM / (15.*a.const.deg)# minutes it takes the sky to drift through beam FWHM
            if param == 'antpos':
                bl_lens = n.sum(n.array(prms[param])**2,axis=1)**.5
        return self.array_params

#===========================ARRAY SPECIFIC PARAMETERS==========================

#Set antenna positions here; for regular arrays like Hera we can use an algorithm; otherwise antpos should just be a list of [x,y,z] coords in light-nanoseconds

xmean,ymean,zmean = 467217.65468,7046613.77961,375.42296875
antpos_dict = {
0: [86.352, 287.221067, 375.032],
1: [100.352, 287.221067, 375.032],
2: [114.352, 287.221067, 375.032],
3: [128.352, 287.221067, 375.032],
4: [79.352, 275.0967113, 375.032],
5: [93.352, 275.0967113, 375.032],
6: [107.352, 275.0967113, 375.032],
7: [121.352, 275.0967113, 375.032],
8: [135.352, 275.0967113, 375.032],
9: [72.352, 262.9723557, 375.032],
10: [86.352, 262.9723557, 375.032],
11: [100.352, 262.9723557, 375.032],
12: [114.352, 262.9723557, 375.032],
13: [128.352, 262.9723557, 375.032],
14: [142.352, 262.9723557, 375.032],
15: [65.352, 250.848, 375.032],
16: [79.352, 250.848, 375.032],
17: [93.352, 250.848, 375.032],
18: [121.352, 250.848, 375.032],
19: [135.352, 250.848, 375.032],
20: [149.352, 250.848, 375.032],
21: [72.352, 238.7236444, 375.032],
22: [86.352, 238.7236444, 375.032],
23: [100.352, 238.7236444, 375.032],
24: [114.352, 238.7236444, 375.032],
25: [128.352, 238.7236444, 375.032],
26: [142.352, 238.7236444, 375.032],
27: [79.352, 226.5992887, 375.032],
28: [93.352, 226.5992887, 375.032],
29: [107.352, 226.5992887, 375.032],
30: [121.352, 226.5992887, 375.032],
31: [135.352, 226.5992887, 375.032],
32: [86.352, 214.474933, 375.032],
33: [100.352, 214.474933, 375.032],
34: [114.352, 214.474933, 375.032],
35: [128.352, 214.474933, 375.032],
36: [-3.49, 156.590067, 376.351],
37: [10.51, 156.590067, 376.351],
38: [24.51, 156.590067, 376.351],
39: [38.51, 156.590067, 376.351],
40: [-10.49, 144.4657113, 376.351],
41: [3.51, 144.4657113, 376.351],
42: [17.51, 144.4657113, 376.351],
43: [31.51, 144.4657113, 376.351],
44: [45.51, 144.4657113, 376.351],
45: [-17.49, 132.3413557, 376.351],
46: [-3.49, 132.3413557, 376.351],
47: [10.51, 132.3413557, 376.351],
48: [24.51, 132.3413557, 376.351],
49: [38.51, 132.3413557, 376.351],
50: [52.51, 132.3413557, 376.351],
51: [-24.49, 120.217, 376.351],
52: [-10.49, 120.217, 376.351],
53: [3.51, 120.217, 376.351],
54: [31.51, 120.217, 376.351],
55: [45.51, 120.217, 376.351],
56: [59.51, 120.217, 376.351],
57: [-17.49, 108.0926444, 376.351],
58: [-3.49, 108.0926444, 376.351],
59: [10.51, 108.0926444, 376.351],
60: [24.51, 108.0926444, 376.351],
61: [38.51, 108.0926444, 376.351],
62: [52.51, 108.0926444, 376.351],
63: [-10.49, 95.96828869, 376.351],
64: [3.51, 95.96828869, 376.351],
65: [17.51, 95.96828869, 376.351],
66: [31.51, 95.96828869, 376.351],
67: [45.51, 95.96828869, 376.351],
68: [-3.49, 83.84393304, 376.351],
69: [10.51, 83.84393304, 376.351],
70: [24.51, 83.84393304, 376.351],
71: [38.51, 83.84393304, 376.351],
72: [-149.785, 265.814, 377.011],
73: [-95.365, 270.188, 377.285],
74: [-88.512, 266.021, 377.309],
75: [-78.698, 258.431, 377.291],
76: [-86.94, 248.452, 377.339],
77: [-102.325, 244.088, 377.258],
78: [-93.346, 242.279, 377.27],
79: [-148.33, 220.271, 377.214],
80: [-132.939, 282.611, 377.002],
81: [-139.867, 275.498, 377.018],
82: [-123.667, 284.854, 377.05],
83: [-92.889, 281.358, 377.254],
84: [-77.835, 282.839, 377.253],
85: [-83.476, 277.563, 377.3],
86: [-76.988, 272.637, 377.298],
87: [-176.006, 289.903, 376.623],
88: [-23.62, 201.546, 376.805],
89: [-23.802, 222.353, 376.872],
90: [-40.489, 235.17, 376.961],
91: [-54.527, 237.518, 377.123],
92: [-64.851, 250.813, 377.239],
93: [-70.075, 259.402, 377.266],
94: [-60.057, 263.676, 377.18],
95: [-50.676, 264.546, 377.1],
96: [-67.99, 271.103, 377.236],
97: [-60.058, 271.905, 377.204],
98: [-70.885, 278.818, 377.252],
99: [-45.54, 273.245, 377.066],
100: [-53.156, 276.268, 377.09],
101: [-22.659, 263.694, 376.948],
102: [-55.097, 284.285, 377.181],
103: [-14.9, 273.87, 376.91],
104: [-105.281, 217.18, 377.255],
105: [-98.53, 230.162, 377.26],
106: [-81.845, 229.337, 377.27],
107: [-79.541, 238.128, 377.27],
108: [-75.12, 247.002, 377.255],
109: [-71.017, 235.929, 377.239],
110: [-62.567, 228.987, 377.21],
111: [-50.691, 221.885, 377.078],
112: [-160.465, 573.207, 374.773],
113: [-128.046, 350.494, 376.57],
114: [-69.771, 294.031, 377.192],
115: [-78.731, 297.259, 377.196],
116: [-103.941, 300.914, 377.05],
117: [-100.172, 288.896, 377.157],
118: [-395.452, 271.152, 374.515],
119: [-263.547, 389.274, 375.322],
120: [326.141, 203.608, 373.945],
121: [173.425, 193.595, 375.032],
122: [-2.515, 292.841, 376.892],
123: [-46.554, 287.469, 377.115],
124: [-56.141, 295.922, 377.156],
125: [-33.705, 351.048, 376.976],
126: [50.127, 536.159, 376.181],
127: [84.88, 519.141, 376.314],
}

antpos = []
for pos in antpos_dict.values():
    antpos.append(n.array(pos)*100./a.const.len_ns)

#Set other array parameters here
prms = {
    'name': os.path.basename(__file__)[:-3], #remove .py from filename
#    'loc': ('38:25:59.24',  '-79:51:02.1'), # Green Bank, WV
    'loc': ('-26:42:11.916', '116:40:15.06'), # MWA,  (GPS)
    'antpos': antpos,
    'beam': a.fit.Beam2DGaussian,
    'dish_size_in_lambda': 2.32, #in units of wavelengths at 150 MHz = 2 meters; this will also define the observation duration
    'Trx': 1e5 #receiver temp in mK
}

#=======================END ARRAY SPECIFIC PARAMETERS==========================

def get_aa(freqs):
    '''Return the AntennaArray to be used for simulation.'''
    location = prms['loc']
    antennas = []
    nants = len(prms['antpos'])
    for i in range(nants):
        beam = prms['beam'](freqs, xwidth=(0.45/prms['dish_size_in_lambda']), ywidth=(0.45/prms['dish_size_in_lambda'])) #as it stands, the size of the beam as defined here is not actually used anywhere in this package, but is a necessary parameter for the aipy Beam2DGaussian object
        antennas.append(a.fit.Antenna(0, 0, 0, beam))
    aa = AntennaArray(prms['loc'], antennas)
    p = {}
    for i in range(nants):
        top_pos = prms['antpos'][i]
        p[str(i)] = {'top_x':top_pos[0], 'top_y':top_pos[1], 'top_z':top_pos[2]}
    aa.set_ant_params(p)
    aa.set_arr_params(prms) 
    return aa

def get_catalog(*args, **kwargs): return a.src.get_catalog(*args, **kwargs)
