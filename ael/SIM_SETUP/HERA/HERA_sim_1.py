
# Holds relevant information for uvfits_zeros_gen to setup a uvfits for a given simulation.
import numpy as n

sim_prms= {
#    'altitude' : '',
#    'phase_center_ra' : '',
#    'phase_center_dec' : '',
#    'object_name' : '',
#    'phase_center_epoch' : '',
  #  'start_time' : 2456528.19962963 + 0.259,   # These times were for MWA when the source catalog was overhead
  #  'end_time' : 2456528.20090279 + 0.259,   # 		Shift by six hours to move the catalog over the HERA site
    'start_time' : 2456528.19962963,  
    'end_time' : 2456528.20090279,
    'Ntimes' : 56,
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
    'instrument' : 'HERA128',
    'xyz_telescope_frame' : "ITRF",
    #'xyz_telescope_frame' : [None, None, None],
    'Nspws' : 1,
    'spw_array' : [1],
    'telescope_name' : 'HERA',
    'integration_time': 8.0    #sec
}
