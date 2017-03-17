#!/bin/env python

from scipy.io import readsav
import numpy as np, sys

print sys.argv[1]

sav = readsav(sys.argv[1])
obs = sav['obs']

bi = obs['BEAM_INTEGRAL'][0]
freqs = np.linspace(1e8,2e8, 203)

fit = np.polyfit(freqs/1e9, bi[0], 8)

print ",".join(map(str,fit.tolist()))
