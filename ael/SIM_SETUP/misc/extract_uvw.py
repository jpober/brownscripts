import numpy as np
import aipy as a


def extract_uvw(fn):
	uv = a.miriad.UV(fn)
	uvws = []
	for p, d in uv.all():
	    uvw, t, (i,j) = p
	    uvws.append(uvw)
	uvws = np.array(uvws)
	return uvws
