import pyuvdata, numpy
from pylab import *

args = sys.argv[1:]
print(args)

lst = np.array([])

for lstfile in args:	
	uv = pyuvdata.UVData()
	uv.read_miriad(lstfile, run_check=False, run_check_acceptability=False)
	waterfall =(uv.nsample_array[np.where(uv.baseline_array == 264257)].squeeze())
	try: 
		data=np.append(data, waterfall, axis=0)
	except NameError:
		data = waterfall
	lstslice = (uv.lst_array[np.where(uv.baseline_array == 264257)].squeeze())
	lst= np.append(lst, lstslice)

indices = numpy.argsort(lst)
lst = lst[indices]
imshow(numpy.abs(data[indices]), interpolation='nearest', aspect='auto', extent=[uv.freq_array[0,0], uv.freq_array[0,-1], lst[-1], lst[0]])

xlabel('Frequency (Hz)')
ylabel('Local Sidereal Time (Radians)')
cb = colorbar()
cb.set_label('Number of Samples')
show()
