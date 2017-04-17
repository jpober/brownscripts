parser = argparse.ArgumentParser()
parser.add_argument('--uvf', '-u', 'action=store_true')
parser.add_argument('n1', type=int)
parser.add_argument('n2', type=int)
parser.add_argument('heradat', nargs=argparse.REMAINDER)

uv = uvdata.UVData()
arr = numpy.empty([0,1024])
for f in heradat:
	if uvf:
		uv.read_uvfits(f)
	else:
		uv.read_miriad(f)
	arr = numpy.concatenate((arr, uv.data_array[numpy.where(uv.baseline_array == uv.antnums_to_baseline(n1, n2))].squeeze()))
pylab.imshow(numpy.log(numpy.absolute(arr)), interpolation='nearest', aspect='auto', vmin=-4, vmax=0)
pylab.colorbar()
pylab.show()
