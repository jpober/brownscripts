
# Holds relevant information for uvfits_zeros_gen to setup a uvfits for a given simulation.
import numpy as n

sim_prms= {
    'altitude' : '',
    'phase_center_ra' : '',
    'phase_center_dec' : '',
    'object_name' : '',
    'phase_center_epoch' : '',
    'start_time' : 2456528.1996296346,
    'end_time' : 2456528.3321990669,
    'Ntimes' : 56,
    'time_format' : 'julian',   # julian, lst

	'channel_width' : 8e4,	 #Hz
#	'sfreq' : 182435000.0, 

	'sfreq' : 166995000.0+8e4,
	'Nfreqs' : 384,

	'Npols' : 4,
	'polarization_array' : 'lin', # lin, circ, stokes

	'vis_units' : 'Jy'

}

instr_prms = {
    'instrument' : 'MWA_Hex',
    'xyz_telescope_frame' : "ITRF",
#    'xyz_telescope' : [None, None, None],
    'Nspws' : 1,
    'spw_array' : [1],
    'telescope_name' : 'MWA',
    'integration_time': 2.3156404495239258e-05
}
