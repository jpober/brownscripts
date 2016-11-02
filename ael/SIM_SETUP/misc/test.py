#!/bin/env python

from uvdata.uv import UVData
import numpy as np
import pylab as plt

uvd0 = UVData()
uvd0.read_uvfits('ewbase_mwa_0.uvfits')

uvd1 = UVData()
uvd1.read_uvfits('ewbase_mwa_1.uvfits')

uvd2 = UVData()
uvd2.read_uvfits('ewbase_mwa_2.uvfits')


uvd0.set_lsts_from_time_array()
print uvd0.phase_center_ra, uvd0.lst_array[0]
#print uvd1.phase_center_ra, uvd1.phase_center_dec
#print uvd2.phase_center_ra, uvd2.phase_center_dec


if False:
	uvd1.unphase_to_drift()
	uvd0.unphase_to_drift()
	uvd2.unphase_to_drift()

#other
#plt.scatter(uvd0.uvw_array[1][1::3], uvd0.uvw_array[2][1::3], color='red')
#plt.scatter(uvd1.uvw_array[1][1::3], uvd1.uvw_array[2][1::3], color='blue')
#plt.scatter(uvd2.uvw_array[1][1::3], uvd2.uvw_array[2][1::3], color='green')


#plt.scatter(uvd0.uvw_array[1], uvd0.uvw_array[2], color='red')
#plt.scatter(uvd1.uvw_array[1], uvd1.uvw_array[2], color='blue')
#plt.scatter(uvd2.uvw_array[1], uvd2.uvw_array[2], color='green')


plt.figure()

plt.subplot(221)
plt.scatter(uvd0.time_array[1::3], uvd0.uvw_array[0][1::3], color='red')
plt.scatter(uvd1.time_array[1::3], uvd1.uvw_array[0][1::3], color='blue')
plt.scatter(uvd2.time_array[1::3], uvd2.uvw_array[0][1::3], color='green')
plt.subplot(222)

plt.scatter(uvd0.time_array[1::3], uvd0.uvw_array[1][1::3], color='red')
plt.scatter(uvd1.time_array[1::3], uvd1.uvw_array[1][1::3], color='blue')
plt.scatter(uvd2.time_array[1::3], uvd2.uvw_array[1][1::3], color='green')

plt.subplot(223)

plt.scatter(uvd0.time_array[1::3], uvd0.uvw_array[2][1::3], color='red')
plt.scatter(uvd1.time_array[1::3], uvd1.uvw_array[2][1::3], color='blue')
plt.scatter(uvd2.time_array[1::3], uvd2.uvw_array[2][1::3], color='green')




#plt.subplot(221)
#plt.scatter(uvd0.time_array, uvd0.uvw_array[0], color='red')
#plt.scatter(uvd1.time_array, uvd1.uvw_array[0], color='blue')
#plt.scatter(uvd2.time_array, uvd2.uvw_array[0], color='green')
#plt.subplot(222)
#
#plt.scatter(uvd0.time_array, uvd0.uvw_array[1], color='red')
#plt.scatter(uvd1.time_array, uvd1.uvw_array[1], color='blue')
#plt.scatter(uvd2.time_array, uvd2.uvw_array[1], color='green')
#
#plt.subplot(223)
#
#plt.scatter(uvd0.time_array, uvd0.uvw_array[2], color='red')
#plt.scatter(uvd1.time_array, uvd1.uvw_array[2], color='blue')
#plt.scatter(uvd2.time_array, uvd2.uvw_array[2], color='green')


plt.show()
