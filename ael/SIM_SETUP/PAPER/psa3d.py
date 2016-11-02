import numpy as np
from psa6622_v003_paper128 import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


exes, wise, zees = [], [], []
for i in prms['antpos']:
    dic = prms['antpos'][i]
    exes.append(dic['top_x'])
    wise.append(dic['top_y'])
    zees.append(dic['top_z'])

fig = plt.figure()
ax0 = fig.add_subplot(111, projection='3d')
#ax1 = fig.add_subplot(221)
#ax2 = fig.add_subplot(222)
#ax3 = fig.add_subplot(223)

ax0.scatter(exes,wise,zees)

#ax1.scatter(exes, wise)
#ax2.scatter(exes, zees)
#ax3.scatter(wise, zees)

plt.show(block=True)
