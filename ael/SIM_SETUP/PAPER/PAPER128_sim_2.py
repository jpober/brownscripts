
# Holds relevant information for uvfits_zeros_gen to setup a uvfits for a given simulation.
import numpy as n

sim_prms= {
    'start_time' : 2456535.333333, 
    'end_time' :   2456535.708333,    # 9 hours total.
    'Ntimes' : 20,    		      #Integrations per file  (ten minutes total)
    'time_format' : 'julian',   # julian

	'channel_width' : 492610.8374384237,	 #Hz
#	'sfreq' : 182435000.0, 

	'sfreq' : int(1e8),   #Hz
	'Nfreqs' : 203,

	'Npols' : 4,
	'polarization_array' : 'lin', # lin, circ, stokes

	'vis_units' : 'Jy'
}


## Frequencies are good. Need to check times with catalog.

instr_prms = {
    'instrument' : 'PAPER',
    'xyz_telescope_frame' : "ITRF",
    #'xyz_telescope_frame' : [None, None, None],
    'Nspws' : 1,
    'spw_array' : [1],
    'telescope_name' : 'PAPER',
    'integration_time': 30.0    #sec
}
