parser = argparse.ArgumentParser()
parser.add_argument('--uvf', '-u', 'action=store_true')
parser.add_argument('n11', type=int)
parser.add_argument('n12', type=int)
parser.add_argument('n21', type=int)
parser.add_argument('n22', type=int)
parser.add_argument('heradat', nargs=argparse.REMAINDER)

uv = uvdata.UVData()
a1 = numpy.empty([0,1024])
a2 = numpy.empty([0,1024])
for f in heradat:
	if uvf:
		uv.read_uvfits(f)
	else:
		uv.read_miriad(f)
	a1 = numpy.concatenate((a1, uv.data_array[numpy.where(uv.baseline_array == uv.antnums_to_baseline(n11, n12))].squeeze()))
	a2 = numpy.concatenate((a2, uv.data_array[numpy.where(uv.baseline_array == uv.antnums_to_baseline(n21, n22))].squeeze()))
a1 = numpy.absolute(a1)
a2 = numpy.absolute(a2)

pylab.subplot(131)
pylab.imshow(a1, interpolation='nearest', aspect='auto')
pylab.subplot(132)
pylab.imshow(a2, interpolation='nearest', aspect='auto')
pylab.imshow(a1/a2, interpolation='nearest', aspect='auto', vmin=-1, vmax=1)
pylab.show()
