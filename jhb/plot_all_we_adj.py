import uvdata, argparse, numpy, pylab

parser = argparse.ArgumentParser()
parser.add_argument('--uvf', '-u', action='store_true')
parser.add_argument('heradat', nargs=argparse.REMAINDER)
args = parser.parse_args()

uv = uvdata.UVData()

reda1 = [72, 112, 105, 22, 81, 88, 9, 20, 89, 64, 53, 31, 80, 104]
reda2 = [112, 97, 22, 81, 10, 9, 20, 89, 43, 53, 31, 65, 104, 96]
weredids = []
for i in range(len(reda1)):
	if reda1[i] > reda2[i]:
		weredids.append(uv.antnums_to_baseline(reda2[i], reda1[i]))
	else:
		weredids.append(uv.antnums_to_baseline(reda1[i], reda2[i]))

plots = []
for i in range(len(weredids)):
	plots.append(numpy.empty([0,1024]))

for f in args.heradat:
	if args.uvf:
		uv.read_uvfits(f)
	else:
		uv.read_miriad(f)
	for i in range(len(weredids)):
		plots[i] = numpy.concatenate((plots[i], uv.data_array[numpy.where(uv.baseline_array == weredids[i])].squeeze()))

for i in range(len(weredids)):
	pylab.subplot(3, 5, i)
	pylab.imshow(numpy.log(numpy.absolute(plots[i])), interpolation='nearest', aspect='auto', vmin=-4, vmax=0)
	pylab.colorbar()
pylab.show()