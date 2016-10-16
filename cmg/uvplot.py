#! /usr/bin/env python

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

parser.add_argument("-m", "--mode", default="log", help="Plot modes: log, lin, phs, real, imag. Default is log")
parser.add_argument("-a", "--antennas", nargs="*", help="Specify which antennas to include. Default is all")
parser.add_argument("-t", "--times", nargs="*", help="Specify which times to include. Default is all")
parser.add_argument("-f", "--freqs", nargs=2, help="Specify the range of frequencies to include. Default is all")
parser.add_argument("-fl", "--flags", help="Denote bad data in plot", action="store_true")
parser.add_argument("-e", "--extents", help="Define extents of plot from physical data", action="store_true")

args = parser.parse_args()

#print args

class UVData(uvdata.UVData):
    def select(self, antnums=[], times=[], freqs=[]):
        # For now, assumes the object has been read into already
        # freqs must be either empty, or have 2 values [lower bound, upper bound]
        baselines = self.antnum_list_to_baselines(antnums)
        if not baselines:
            baselines = set(self.baseline_array)
        
        baseline_map = {k:[] for k in baselines}
        for i in range(len(self.baseline_array)):
            if self.baseline_array[i] in baselines:
                baseline_map[self.baseline_array[i]].append(i)
    
        if not times:
            times = list(set(self.time_array))
        
        freq_ind = []
        if not freqs:
            freq_ind = range(len(self.freq_array[0]))
        else:
            for i in range(len(self.freq_array[0])):
                if freqs[0] <= self.freq_array[0][i] <= freqs[1]:
                    freq_ind.append(i)
        
        d = {bl:{self.time_array[k]:np.array([self.data_array[k][0][j][0]
                                             if not self.flag_array[k][0][j][0]
                                             else 0
                                             for j in freq_ind], dtype=np.complex64)
                for k in baseline_map[bl]
                if self.time_array[k] in times}
            for bl in baseline_map.keys()}

        #for bl in baseline_map.keys():
        #    for k in baseline_map[bl]:
        #        if self.time_array[k] in times:
        #            for j in freq_ind:
        #                self.data_array[k][0][j][0]

        return (d, sorted(times), sorted([self.freq_array[0][i] for i in freq_ind]))
    def antnum_list_to_baselines(self, antnums=[]):
        '''
        antnums will be a list of either tuples of strings, or strings
        this implementation allows the user to input both 0_1 and 1_0
        and it will return the expected baseline (0_1) in both cases
        '''
        antnums_in_data = set(self.ant_1_array)
        baselines = set()
        
        for i in antnums:
            if isinstance(i, tuple):
                ant1, ant2 = np.int64(i)
                
                if ant1 not in antnums_in_data:
                    raise ValueError('No antenna {} found in data.'.format(ant1))
                if ant2 not in antnums_in_data:
                    raise ValueError('No antenna {} found in data.'.format(ant2))
                
                baselines.add(self.antnums_to_baseline(min(ant1, ant2), max(ant1, ant2)))
            
            else:
                ant = np.int64(i)

                if ant not in antnums_in_data:
                    raise ValueError('No antenna {} found in data.'.format(ant))
                
                for j in antnums_in_data:
                    baselines.add(self.antnums_to_baseline(min(ant, j), max(ant, j)))

        return baselines


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
	if mode == "lin":
		return np.ma.absolute(data)
	if mode == "log":
		data = np.ma.masked_less_equal(np.ma.absolute(data), 0)
        return np.ma.log10(data)
	if mode == "phs":
		return np.angle(data)
	if mode == "real":
		return np.real(data)
	if mode == "imag":
		return np.imag(data)

def find_min(data):
	minimum = min(data[0])
	for i in range(1, len(data)):
		p = min(data[i])
		if p < minimum:
			minimum = p
	return minimum

def find_max(data):
	maximum = max(data[0])
	for i in range(1, len(data)):
		p = max(data[i])
		if p > maximum:
			maximum = p
	return maximum



antenna_list = create_antenna_list(args.antennas)
times_list = create_times_list(args.times)
freqs_list = create_freqs_list(args.freqs)

obj = UVData()
obj.read_miriad(args.data[0])
if args.flags:
	print "using flag array!"
data = obj.select(antenna_list, times_list, freqs_list, args.flags)
# data is in the form (dictionary, list of times in dictionary, list of freqs in dictionary)
#print data

m2 = int(np.sqrt(len(data[0].keys())))
m1 = int(np.ceil(float(len(data[0].keys())) / m2))

for i in range(len(data[0].keys())):
	b = data[0].keys()[i]
	Z = [[data[0][b][t][f] for f in range(len(data[2]))] for t in data[1]]

	Z = np.ma.masked_invalid(Z)
	Z = data_mode(Z, args.mode)

	# masked_invalid creates a masked array, masks where data is NaN or inf
	# in select function, data flagged by flag_array is replaced with NaN, if
	# useflags is set to True
	print "{}\tMIN: {}".format(b, find_min(Z))
	print "{}\tMAX: {}".format(b, find_max(Z))
	plt.subplot(m2, m1, i+1)
	cmap = plt.get_cmap('jet')
	#if args.flags:
	#	cmap.set_bad(color='white', alpha=1.0)
	
	if args.extents:
		im = plt.imshow(Z,
					extent=[data[2][0], data[2][-1], data[1][0], data[1][-1]],
		            origin='upper', aspect='auto', interpolation='nearest', 
				    vmin=find_min(Z), vmax=find_max(Z), cmap=cmap)
	else:
		im = plt.imshow(Z,
		            origin='upper', aspect='auto', interpolation='nearest', 
				    vmin=find_min(Z), vmax=find_max(Z), cmap=cmap)
	#if args.flags:
		#bad_data = np.ma.masked_where(~Z.mask, Z.mask)
		#print "BAD DATA MASK"
		#print bad_data.mask
		#plt.imshow(Z.mask, 
		#            origin='upper', aspect='auto', interpolation='none', 
		#		    vmin=find_min(Z), vmax=find_max(Z), cmap=cm.gray)
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

