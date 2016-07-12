pro eor_firstpass_versions
except=!except
!except=0
heap_gc 

; wrapper to contain all the parameters for various runs we might do
; using firstpass.

; parse command line args
compile_opt strictarr
args = Command_Line_Args(count=nargs)
;obs_id = args[0]
;obs_id = '1061316296'
output_directory = args[0]
;output_directory = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/'
version = args[1]
;version = 'nb_autocal'
cmd_args={version:version}

;Options for using the IDL profiler
profile=0
profile_path='/dev/null'

; Set default values for everything
calibrate_visibilities=1
deconvolve=1
;;just added
snapshot_healpix_export=1
model_visibilities=1
calibration_visibilities_subtract=0
return_cal_visibilities=1
;;
min_cal_baseline=50.0
max_cal_baselin=172.937
time_avg=1
min_solns=5
ref_antenna=1
ref_antenna_name=2
conv_thresh=0.001
polyfit=0.0
bandpass=0
n_cal_src=518
n_pol=1
beam_threshold=0.05
max_iter=500.0
max_sources=15000.0
check_iter=7
gain_factor=0.15
add_threshold=0.8
pol_use=0
sigma_cut=2
local_max_radius=3
filter_background=1
decon_filter='filter_uv_natural'
decon_mode='Single snapshot'
catalog_file_path=filepath('MRC_full_radio_catalog.fits',root=rootdir('FHD'),subdir='catalog_data')
calibration_catalog_file_path=filepath('mwa_calibration_source_list.sav',root=rootdir('FHD'),subdir='catalog_data')
diffuse_calibrate=filepath('EoR0_diffuse_model_94.sav',root=rootdir('FHD'),subdir='catalog_data')

case version of
   'jrk': begin
      dimension=1024
      elements=1024
      nbaselines=496
      kpix=0.358099
      degpix=0.156250
      n_pol=1
      n_tile=64
      n_freq=150
      n_time=10
      pol_names='xx'
   end
endcase

;old version here -- for MIT
;SPAWN, 'read_uvfits_loc.py -v ' + STRING(uvfits_version) + ' -s ' + $
;  STRING(uvfits_subversion) + ' -o ' + STRING(obs_id), vis_file_list

; On Oscar
;SPAWN, 'locate_uvfits_oscar.py -o ' + STRING(obs_id), vis_file_list


IF (profile eq 1) THEN BEGIN
   RESOLVE_ROUTINE, 'general_obs',/QUIET	; Profiler only looks at compiled modules...
   RESOLVE_ROUTINE, 'slurm_ps_job', /QUIET
   RESOLVE_ALL,/CONTINUE_ON_ERROR,/QUIET
   PROFILER
;   PROFILER, /SYSTEM
ENDIF
;print,vis_file_list
;vis_file_list=vis_file_list ; this is silly, but it's so var_bundle sees it.
;data_directory= '/users/alanman/data/alanman'
;vis_file_list=[data_directory+'/1061311664.uvfits']
;print, vis_file_list[0]
vis_file_list='/users/jkerriga/zen.2456242.29909.uvcRREcACOcPM.uvfits'
undefine,uvfits_version ; don't need these passed further
undefine,uvfits_subversion
undefine,obs_id

fhd_file_list=fhd_path_setup(vis_file_list,version=version,output_directory=output_directory)
healpix_path=fhd_path_setup(output_dir=output_directory,subdir='Healpix',output_filename='Combined_obs',version=version)

extra=var_bundle() ; bundle all the variables into a structure

print,""
print,"Keywords set in wrapper:"
print,structure_to_text(extra)
print,""
general_obs,_Extra=extra

IF (profile eq 1) THEN BEGIN
   PROFILER, FILENAME=STRING(profile_path), /REPORT, /CODE_COVERAGE
ENDIF

end
