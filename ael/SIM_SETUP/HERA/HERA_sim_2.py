
# Holds relevant information for uvfits_zeros_gen to setup a uvfits for a given simulation.
import numpy as n

sim_prms= {
    'start_time' : 2456528.19962963,  
    'end_time' : 2456528.3321990669,      # Three hours. This will be divided into 94 obs ids, 56 time samples per obs 
    'Ntimes' : 56,      		   # Time steps per Obs ID
    'time_format' : 'julian',   # julian

	'channel_width' : 5e4,	 #Hz
#	'sfreq' : 182435000.0, 

	'sfreq' : int(100e6),   #Hz
	'Nfreqs' : 200,   #Arbitrary.

	'Npols' : 4,
	'polarization_array' : 'lin', # lin, circ, stokes

	'vis_units' : 'Jy'

}

instr_prms = {
    'instrument' : 'HERA37',
    'xyz_telescope_frame' : "ITRF",
    #'xyz_telescope_frame' : [None, None, None],
    'Nspws' : 1,
    'spw_array' : [1],
    'telescope_name' : 'HERA',
    'integration_time': 8.0    #sec
}
