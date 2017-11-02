import uvdata
import numpy as n
from numpy.random import randn

a = uvdata.miriad.Miriad()
try:
    a.read_miriad('/users/jkerriga/Pzen.2456242.30605.uvcRREcACOTUcHPA')
except:
    pass

shp = a.data_array.shape
noise = (10**7)*randn(shp[0],shp[1],shp[2],shp[3])+1j*(10**7)*randn(shp[0],shp[1],shp[2],shp[3])
NoiseSig = a.data_array+noise
a.epoch=2000.0
a.vis_units='Jy'
a.phase_center_epoch=2000.0
a.write_miriad('ZeroDayNoiseInjection.uv')
a.data_array = noise
a.write_miriad('JustNoise.uv')
