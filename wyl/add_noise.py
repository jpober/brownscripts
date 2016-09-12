import numpy, ephem, uvdata.uvdata as uvd, sys, math
import glob, matplotlib.pyplot as plt
from IPython import embed

antpos = {
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
}

for ii in range(0,72):
    v = numpy.array([antpos[ii]['top_x']+24.49,antpos[ii]['top_y'],antpos[ii]['top_z']])
    antpos[ii] = v

class bltype:
    def __int__(self):
        self.x = 0
        self.y = 0
        self.z = 0
        self.type = ((0,0,0))
    def read(self,vec):
        self.x = vec[0]
        self.y = vec[1]
        self.z = vec[2]
        self.type = ((round(vec[0]),round(vec[1]),round(vec[2])))
    def __eq__(self,other):
        if isinstance(other, self.__class__):
            thred = 0.5
            diff = (self.x-other.x)*(self.x-other.x)+(self.y-other.y)*(self.y-other.y)+(self.z-other.z)*(self.z-other.z)
            #diff2 = (self.x+other.x)*(self.x+other.x)+(self.y+other.y)*(self.y+other.y)+(self.z+other.z)*(self.z+other.z)
            if diff < thred: return True
            else: return False
        else:
            print('Classes do not match')
            return False
    def __hash__(self):
        xx = round(self.x)
        yy = round(self.y)
        zz = round(self.z)
        return hash((round(xx),round(yy),round(zz)))
        
def magnitude(vec):
    return vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]

#Area = 20.48 # m2
#sigma = kb*Tsys/(Area*math.sqrt(df*dt))*1e26/1.4142135623730951
def sigma(frequency):
    kb = 1.38064852e-23 # m2kg/s2/K
    Tsys = 210.  # K
    df = 80000.  # Hz
    dt = 2.  # s
    Area = 21.5+(19.8-21.5)*(frequency-1.5e8)/(2.0e8-1.5e8) # m2
    return kb*Tsys/(Area*math.sqrt(df*dt))*1e26/1.4142135623730951 # Jy

obs = sys.argv[1] + '*'
flist = glob.glob(obs)
print '   reading:'
print flist
uv = uvd.UVData()
uv.read_fhd(flist)
print '   adding noise:'
for kk in range(0,uv.Npols):
    data = uv.data_array[:,0][:,:,kk]
    for ii in range(0,uv.Nfreqs):
        sig = sigma(uv.freq_array[0][ii])
        print sig
        sre = numpy.random.normal(0, sig, uv.Nblts)
        sim = numpy.random.normal(0, sig, uv.Nblts)
        data[:,ii] += (sre + sim*1j)
ofn = sys.argv[1] + '.uvfits'
print '   writing...'
uv.write_uvfits(ofn,spoof_nonessential=True)




