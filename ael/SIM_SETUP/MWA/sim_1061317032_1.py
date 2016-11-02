
# Holds relevant information for uvfits_zeros_gen to setup a uvfits for a given simulation.
import numpy as n

sim_prms= {
    'start_time' : 2456528.2617592663,
    'end_time' : 2456528.2630324066,
    'Ntimes' : 56,
    'time_format' : 'julian',   # julian, lst

	'channel_width' : 8e4,	 #Hz
#	'sfreq' : 182435000.0, 
	#'sfreq' : 91217500.0, 
	'sfreq' : 1.82435E+08,
	'Nfreqs' : 193,

	'Npols' : 4,
	'polarization_array' : 'lin', # lin, circ, stokes

	'vis_units' : 'Jy'

}

instr_prms = {
    'instrument' : 'MWA',
    'xyz_telescope_frame' : "ITRF",
    'xyz_telescope' : [None, None, None],
    'Nspws' : 1,
    'spw_array' : [1],
    'telescope_name' : 'MWA',
#    'integration_time': 2.3156404495239258e-05   #!!!! Should be in seconds!!!
    'integration_time': 2.000   
}
