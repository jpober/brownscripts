#!/bin/env python

import sys
from numpy import ceil
args = sys.argv[1:]

exec('from %s import sim_prms as _sim_params' % args[0].split('.')[0] )
exec('from %s import instr_prms as _instr_params' % args[0].split('.')[0] )

inttime=_instr_params['integration_time']
end=_sim_params['end_time']
start=_sim_params['start_time']
Ntimes=_sim_params['Ntimes']

dayspersecond=1./(24.*60.**2)

print int(ceil((end - start)/(inttime*Ntimes*dayspersecond))) - 1
