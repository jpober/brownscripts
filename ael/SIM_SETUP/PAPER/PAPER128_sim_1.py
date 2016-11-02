
# Holds relevant information for uvfits_zeros_gen to setup a uvfits for a given simulation.
import numpy as n

sim_prms= {
#    'start_time' : 2456528.19962963, 
#    'end_time' : 2456528.2412962965,    # 1 hour total --- original MWA golden timing
    'start_time' : 2456527.50000000, 
    'end_time' : 2456527.54166666,    # 1 hour total -- shifted to midnight locally
    'Ntimes' : 4,    		      #Integrations per file  (two minutes total)
    'time_format' : 'julian',   # julian

	'channel_width' : 4.8828125e4,	 #Hz   203 channels, 100 to 200 MHz
#	'sfreq' : 182435000.0, 

	'sfreq' : int(100e6),   #Hz
	'Nfreqs' : 203,   		#Select 20 course channels, 95 to 115 of the 203 between 100 and 200 MHz

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
