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
    11:	[467101.60-xmean,	7046644.69-ymean,	377.02-zmean],
    12:	[467155.98-xmean,	7046649.21-ymean,	377.29-zmean],
    13:	[467162.84-xmean,	7046645.06-ymean,	377.32-zmean],
    14:	[467172.67-xmean,	7046637.50-ymean,	377.30-zmean],
    15:	[467164.45-xmean,	7046627.50-ymean,	377.34-zmean],
    16:	[467149.09-xmean,	7046623.10-ymean,	377.26-zmean],
    17:	[467158.07-xmean,	7046621.32-ymean,	377.28-zmean],
    18:	[467103.17-xmean,	7046599.18-ymean,	377.22-zmean],
    21:	[467118.39-xmean,	7046661.53-ymean,	377.01-zmean],
    22:	[467111.48-xmean,	7046654.40-ymean,	377.03-zmean],
    23:	[467127.65-xmean,	7046663.79-ymean,	377.06-zmean],
    24:	[467158.42-xmean,	7046660.38-ymean,	377.26-zmean],
    25:	[467173.47-xmean,	7046661.90-ymean,	377.26-zmean],
    26:	[467167.84-xmean,	7046656.61-ymean,	377.31-zmean],
    27:	[467174.34-xmean,	7046651.70-ymean,	377.30-zmean],
    28:	[467075.33-xmean,	7046668.70-ymean,	376.63-zmean],
    31:	[467227.86-xmean,	7046580.78-ymean,	376.81-zmean],
    32:	[467227.63-xmean,	7046601.58-ymean,	376.88-zmean],
    33:	[467210.92-xmean,	7046614.35-ymean,	376.97-zmean],
    34:	[467196.88-xmean,	7046616.66-ymean,	377.13-zmean],
    35:	[467186.53-xmean,	7046629.92-ymean,	377.24-zmean],
    36:	[467181.28-xmean,	7046638.49-ymean,	377.27-zmean],
    37:	[467191.28-xmean,	7046642.79-ymean,	377.19-zmean],
    38:	[467200.66-xmean,	7046643.68-ymean,	377.11-zmean],
    41:	[467183.34-xmean,	7046650.19-ymean,	377.24-zmean],
    42:	[467191.26-xmean,	7046651.01-ymean,	377.21-zmean],
    43:	[467180.42-xmean,	7046657.89-ymean,	377.26-zmean],
    44:	[467205.77-xmean,	7046652.39-ymean,	377.07-zmean],
    45:	[467198.15-xmean,	7046655.39-ymean,	377.10-zmean],
    46:	[467228.66-xmean,	7046642.90-ymean,	376.95-zmean],
    47:	[467196.19-xmean,	7046663.40-ymean,	377.19-zmean],
    48:	[467236.39-xmean,	7046653.09-ymean,	376.92-zmean],
    51:	[467196.14-xmean,	7046504.00-ymean,	376.80-zmean],
    52:	[467429.48-xmean,	7046336.50-ymean,	375.01-zmean],
    53:	[467230.09-xmean,	7046456.31-ymean,	376.35-zmean],
    54:	[467343.96-xmean,	7046299.11-ymean,	375.76-zmean],
    55:	[467406.28-xmean,	7046120.00-ymean,	376.48-zmean],
    56:	[467189.51-xmean,	7046361.89-ymean,	376.18-zmean],
    57:	[467134.12-xmean,	7046339.24-ymean,	376.06-zmean],
    58:	[467164.20-xmean,	7046504.31-ymean,	377.03-zmean],
    61:	[467146.20-xmean,	7046596.20-ymean,	377.26-zmean],
    62:	[467152.92-xmean,	7046609.19-ymean,	377.27-zmean],
    63:	[467169.60-xmean,	7046608.41-ymean,	377.28-zmean],
    64:	[467171.88-xmean,	7046617.20-ymean,	377.28-zmean],
    65:	[467176.27-xmean,	7046626.08-ymean,	377.26-zmean],
    66:	[467180.40-xmean,	7046615.03-ymean,	377.24-zmean],
    67:	[467188.87-xmean,	7046608.11-ymean,	377.21-zmean],
    68:	[467200.75-xmean,	7046601.04-ymean,	377.08-zmean],
    71:	[467055.62-xmean,	7046568.59-ymean,	377.17-zmean],
    72:	[467053.65-xmean,	7046473.99-ymean,	377.03-zmean],
    73:	[466968.19-xmean,	7046182.59-ymean,	374.91-zmean],
    74:	[466864.15-xmean,	7046334.40-ymean,	375.38-zmean],
    75:	[466799.58-xmean,	7046362.31-ymean,	375.13-zmean],
    76:	[466782.93-xmean,	7046563.28-ymean,	374.77-zmean],
    77:	[466887.18-xmean,	7046561.47-ymean,	375.75-zmean],
    78:	[466881.88-xmean,	7046644.80-ymean,	375.26-zmean],
    81:	[467090.13-xmean,	7046951.90-ymean,	374.80-zmean],
    82:	[467123.11-xmean,	7046729.39-ymean,	376.58-zmean],
    83:	[467181.50-xmean,	7046673.10-ymean,	377.20-zmean],
    84:	[467172.53-xmean,	7046676.31-ymean,	377.20-zmean],
    85:	[467147.33-xmean,	7046679.89-ymean,	377.06-zmean],
    86:	[467151.12-xmean,	7046667.89-ymean,	377.16-zmean],
    87:	[466855.78-xmean,	7046749.34-ymean,	374.54-zmean],
    88:	[466987.57-xmean,	7046767.80-ymean,	375.34-zmean],
    91:	[467577.44-xmean,	7046583.75-ymean,	373.96-zmean],
    92:	[467424.83-xmean,	7046573.35-ymean,	375.04-zmean],
    93:	[467248.72-xmean,	7046672.09-ymean,	376.90-zmean],
    94:	[467204.72-xmean,	7046666.60-ymean,	377.12-zmean],
    95:	[467195.12-xmean,	7046675.03-ymean,	377.16-zmean],
    96:	[467217.40-xmean,	7046730.18-ymean,	376.99-zmean],
    97:	[467300.71-xmean,	7046915.42-ymean,	376.20-zmean],
    98:	[467335.49-xmean,	7046898.50-ymean,	376.34-zmean],
   101:	[466732.58-xmean,	7047010.25-ymean,	372.66-zmean],
   102:	[466760.74-xmean,	7046982.40-ymean,	372.96-zmean],
   103:	[466675.66-xmean,	7046792.71-ymean,	373.41-zmean],
   104:	[466666.88-xmean,	7046276.41-ymean,	375.24-zmean],
   105:	[466598.49-xmean,	7046519.90-ymean,	374.50-zmean],
   106:	[466566.76-xmean,	7046589.69-ymean,	374.28-zmean],
   107:	[466172.72-xmean,	7046463.59-ymean,	377.10-zmean],
   108:	[466224.96-xmean,	7046950.99-ymean,	374.15-zmean],
   111:	[466914.12-xmean,	7048038.43-ymean,	371.29-zmean],
   112:	[467022.64-xmean,	7047276.49-ymean,	373.76-zmean],
   113:	[467076.35-xmean,	7047160.19-ymean,	374.27-zmean],
   114:	[466992.87-xmean,	7047270.29-ymean,	373.56-zmean],
   115:	[467055.10-xmean,	7047141.40-ymean,	374.14-zmean],
   116:	[467000.58-xmean,	7047110.09-ymean,	373.84-zmean],
   117:	[466949.63-xmean,	7047082.59-ymean,	373.56-zmean],
   118:	[466360.16-xmean,	7047464.69-ymean,	370.54-zmean],
   121:	[467964.39-xmean,	7047785.79-ymean,	374.80-zmean],
   122:	[467538.32-xmean,	7046998.10-ymean,	375.17-zmean],
   123:	[467488.20-xmean,	7047148.04-ymean,	375.86-zmean],
   124:	[467482.65-xmean,	7047199.30-ymean,	375.81-zmean],
   125:	[467311.30-xmean,	7047194.40-ymean,	376.05-zmean],
   126:	[467191.06-xmean,	7047189.37-ymean,	375.22-zmean],
   127:	[467284.75-xmean,	7047329.25-ymean,	376.35-zmean],
   128:	[467249.36-xmean,	7047791.05-ymean,	373.84-zmean],
   131:	[468088.95-xmean,	7047371.29-ymean,	372.83-zmean],
   132:	[468426.47-xmean,	7047240.63-ymean,	372.29-zmean],
   133:	[468563.51-xmean,	7046935.88-ymean,	370.51-zmean],
   134:	[467807.20-xmean,	7046625.20-ymean,	372.31-zmean],
   135:	[467822.32-xmean,	7046739.90-ymean,	372.55-zmean],
   136:	[467739.50-xmean,	7046805.43-ymean,	373.00-zmean],
   137:	[467778.93-xmean,	7046919.69-ymean,	372.83-zmean],
   138:	[467631.18-xmean,	7047032.92-ymean,	374.30-zmean],
   141:	[468351.99-xmean,	7046533.93-ymean,	369.16-zmean],
   142:	[468453.66-xmean,	7045939.48-ymean,	368.46-zmean],
   143:	[467751.05-xmean,	7046116.23-ymean,	373.25-zmean],
   144:	[467553.16-xmean,	7046331.25-ymean,	374.63-zmean],
   145:	[467622.55-xmean,	7046341.96-ymean,	374.35-zmean],
   146:	[467796.27-xmean,	7046289.70-ymean,	372.84-zmean],
   147:	[467758.54-xmean,	7046377.49-ymean,	372.82-zmean],
   148:	[467775.26-xmean,	7046456.21-ymean,	372.64-zmean],
   151:	[467908.11-xmean,	7045343.87-ymean,	371.19-zmean],
   152:	[467696.91-xmean,	7045545.66-ymean,	372.32-zmean],
   153:	[467194.52-xmean,	7045874.17-ymean,	372.93-zmean],
   154:	[467290.08-xmean,	7045974.65-ymean,	374.05-zmean],
   155:	[467485.64-xmean,	7046022.29-ymean,	375.28-zmean],
   156:	[467510.79-xmean,	7046013.80-ymean,	375.24-zmean],
   157:	[467578.34-xmean,	7046061.90-ymean,	375.05-zmean],
   158:	[467578.72-xmean,	7046074.09-ymean,	375.14-zmean],
   161:	[466840.80-xmean,	7045969.48-ymean,	374.00-zmean],
   162:	[467111.25-xmean,	7045958.99-ymean,	373.52-zmean],
   163:	[466774.07-xmean,	7045561.25-ymean,	375.79-zmean],
   164:	[466678.98-xmean,	7045514.57-ymean,	375.61-zmean],
   165:	[466579.96-xmean,	7045739.12-ymean,	375.41-zmean],
   166:	[466621.54-xmean,	7046132.79-ymean,	375.58-zmean],
   167:	[466619.51-xmean,	7046208.49-ymean,	375.65-zmean],
   168:	[466684.41-xmean,	7046120.40-ymean,	375.28-zmean],
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
