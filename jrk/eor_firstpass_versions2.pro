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
vis_file_list='/users/jkerriga/zen.2456242.29909.uvcRREcACOcPM.uvfits'
fhd_file_list=fhd_path_setup(vis_file_list,version=version)
data_directory=output_directory
healpix_path=fhd_path_setup(output_dir=data_directory,subdir='Healpix',output_filename='Combined_obs',version=version)
catalog_file_path=filepath('MRC_full_radio_catalog.fits',root=rootdir('FHD'),subdir='catalog_data')
calibration_catalog_file_path=filepath('mwa_commissioning_source_list_add_FHDaug23deconvolve_fornax_and_VLA_pic.sav',root=rootdir('FHD'),subdir='catalog_data')
model_catalog_file_path=filepath('mwa_commissioning_source_list_add_FHDaug23deconvolve_fornax_and_VLA_pic.sav',root=rootdir('FHD'),subdir='catalog_data')

;Options for using the IDL profiler
profile=0
profile_path='/dev/null'
image_filter_fn='filter_uv_uniform'
; Set default values for everything
calibrate_visibilities=1
deconvolve=0
firstpass=1
;;just added
snapshot_healpix_export=1
model_visibilities=1
calibration_visibilities_subtract=0
return_cal_visibilities=1
;delays='Null pointer'
FoV=160.
reorder_visibilities=1
precess=1 ;At least the new uvfits files from Jonnie need to be precessed to J2000
instrument='paper'
lat=Ten(-30,42,17.5)
lon=Ten(21,25,41)
time_offset=5.*60.
min_baseline=1.
min_cal_baseline=20.
bandpass_calibrate=0
calibration_polyfit=0.
no_restrict_cal_sources=1
no_rephase=1
freq_start=124.
freq_end=174.
beam_model_version=2
beam_cal_threshold=0.5
cal_gain_init=1E-02
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
export_images=1
;decon_filter='filter_uv_natural'
;decon_mode='Single snapshot'


case version of
   'jrk': begin
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
