
# Holds relevant information for uvfits_zeros_gen to setup a uvfits for a given simulation.
import numpy as n

sim_prms= {
    'start_time' : 2456528.2434529355, 
    'end_time' : 2456528.2644212991, 
    'Ntimes' : 56,    		      #Integrations per file
    'time_format' : 'julian',   # julian
    'file_gap' : 9.5,    #Designated time between the end of one file and the start of another, in seconds
	'channel_width' : 80000.0,	 #Hz
#	'sfreq' : 182435000.0, 
	'sfreq' : 167075008.0,   #Hz
	'Nfreqs' : 384,
	'Npols' : 4,
	'polarization_array' : 'lin', # lin, circ, stokes
	'vis_units' : 'Jy'
}


## Frequencies are good. Need to check times with catalog.

instr_prms = {
    'instrument' : 'MWA',
    'xyz_telescope_frame' : "ITRF",
    #'xyz_telescope_frame' : [None, None, None],
    'Nspws' : 1,
    'spw_array' : [1],
    'telescope_name' : 'MWA',
    'integration_time': 2.0    #sec
}
