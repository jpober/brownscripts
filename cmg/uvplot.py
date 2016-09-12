from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

from matplotlib.collections import PolyCollection
from matplotlib.colors import colorConverter

from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np

import argparse
import uvdata


parser = argparse.ArgumentParser(description='Describe this thing here...')
parser.add_argument("data", nargs="+", help="path(s) to data")

parser.add_argument("-m", "--mode", default="amp", help="Plot modes: amp, phs, real, imag. Default is amp")
parser.add_argument("-a", "--antennas", nargs="*", help="Specify which antennas to include. Default is all")
parser.add_argument("-t", "--times", nargs="*", help="Specify which times to include. Default is all")
parser.add_argument("-f", "--freqs", nargs=2, help="Specify the range of frequencies to include. Default is all")

args = parser.parse_args()

print args



def create_antenna_list(antennas):
	ant_list = []
	for ant in antennas:
		s = ant.split("_")
		if len(s) == 2:
			ant_list.append((s[0], s[1]))
		elif len(s) == 1:
			ant_list.append(s[0])
		else:
			raise ValueError("Bad antenna input: {}".format(ant))
	return ant_list

def create_times_list(times):
	if not times:
		return []
	return [np.float64(t) for t in times]

def create_freqs_list(freqs):
	if not freqs:
		return []
	return [np.float64(f) for f in freqs]

def data_mode(data, mode):
	#plt.title(mode)
	if mode == "amp":
		return np.abs(data)
	if mode == "phs":
		return np.angle(data)
	if mode == "real":
		return np.real(data)
	if mode == "imag":
		return np.imag(data)

def find_min(data):
	mini = min(data[0])
	for i in range(1, len(data)):
		p = min(data[i])
		if p < mini:
			mini = p
	return mini

def find_max(data):
	maxi = max(data[0])
	for i in range(1, len(data)):
		p = max(data[i])
		if p > maxi:
			maxi = p
	return maxi



antenna_list = create_antenna_list(args.antennas)
times_list = create_times_list(args.times)
freqs_list = create_freqs_list(args.freqs)

obj = uvdata.UVData()
obj.read_miriad(args.data[0])
data = obj.select(antenna_list, times_list, freqs_list)
# data is in the form (dictionary, list of times in dictionary, list of freqs in dictionary)
#print data

m2 = int(np.sqrt(len(data[0].keys())))
m1 = int(np.ceil(float(len(data[0].keys())) / m2))

for i in range(len(data[0].keys())):
	b = data[0].keys()[i]
	Z = [[data[0][b][t][f] for f in range(len(data[2]))] for t in data[1]]
	Z = data_mode(Z, args.mode)
	print "{}\tMIN: {}".format(b, find_min(Z))
	print "{}\tMAX: {}".format(b, find_max(Z))
	plt.subplot(m2, m1, i+1)
	im = plt.imshow(Z, extent=[data[2][0], data[2][-1], data[1][0], data[1][-1]],
		            origin='upper', aspect='auto', interpolation='nearest', 
				    vmin=find_min(Z), vmax=find_max(Z), cmap=cm.jet)
	plt.colorbar(shrink=0.5)
	plt.title(obj.baseline_to_antnums(b))
	plt.xlabel("freq")
	plt.ylabel("time")

plt.show()

#fig = plt.figure()
#ax = fig.gca(projection='3d')
#X = data[1] # times
#Y = data[2] # freqs
#X, Y = np.meshgrid(X, Y)

#baseline = data[0].keys()[0]
#Z = []
#Z = [[np.abs(freq) for freq in data[0][baseline][time]] for time in data[0][baseline].keys()]
#for freq in range(len(Y)):
#	Z.append([np.abs(data[0][baseline][time][freq]) for time in data[0][baseline].keys()])

#print "X: {} by {}".format(len(X), len(X[0]))
#print "Y: {} by {}".format(len(Y), len(Y[0]))
#print "Z: {} by {}".format(len(Z), len(Z[0]))

#surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet,
#                       linewidth=0, antialiased=True)

#ax.set_zlim(find_min(Z) - 0.01, find_max(Z) + 0.01)

#ax.set_xlabel("time")
#ax.set_ylabel("freq")
#ax.set_zlabel("amp")

#ax.zaxis.set_major_locator(LinearLocator(10))
#ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

#fig.colorbar(surf, shrink=0.5, aspect=5)
