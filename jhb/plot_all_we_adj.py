parser = argparse.ArgumentParser()
parser.add_argument('--uvf', '-u', 'action=store_true')
parser.add_argument('heradat', nargs=argparse.REMAINDER)

uv = uvdata.UVData()

reda1 = [72, 112, 105, 22, 81, 88, 9, 20, 89, 64, 53, 31, 80, 104]
reda2 = [112, 97, 22, 81, 10, 9, 20, 89, 43, 53, 31, 65, 104, 96]
weredids = uv.antnums_to_baseline(reda1 reda2)
for i in range(len(reda1)):
	weredids[i] = uv.antnums_to_baseline(reda2[i] reda1[i]) if reda1[i] > reda2[i]

if uvf:
	uv.read_uvfits(heradat)
else:
	uv.read_miriad(heradat)
for i in range(len(weredids)):
	pylab.subplot(3 5 i)
	pylab.imshow(numpy.log(numpy.absolute(uv.data_array[numpy.where(uv.baseline_array == weredids[i])].squeeze())), interpolation='nearest', aspect='auto', vmin=-3, vmax=0)
	pylab.colorbar()
pylab.show()