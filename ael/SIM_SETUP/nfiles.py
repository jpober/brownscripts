#!/bin/env python

import sys
from numpy import ceil
args = sys.argv[1:]

## Determine the number of files into which the simulation will be split.

exec('from %s import sim_prms as _sim_params' % args[0].split('.')[0] )
exec('from %s import instr_prms as _instr_params' % args[0].split('.')[0] )

inttime=_instr_params['integration_time']
end=_sim_params['end_time']
start=_sim_params['start_time']
Ntimes=_sim_params['Ntimes']
try:
	file_gap = _sim_params['file_gap']
except KeyError:
	file_gap=0.0
dayspersecond=1./(24.*60.**2)

print int(ceil((end - start)/((inttime*Ntimes+file_gap)*dayspersecond))) - 1
