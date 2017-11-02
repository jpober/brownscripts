PRO hera_hsa19,recalculate_all=recalculate_all,export_images=export_images,version=version,_Extra=extra
except=!except
!except=0 
heap_gc
;version='test'
;IF N_Elements(recalculate_all) EQ 0 THEN recalculate_all=1
;IF N_Elements(export_images) EQ 0 THEN export_images=1
;IF N_Elements(cleanup) EQ 0 THEN cleanup=0
;IF N_Elements(ps_export) EQ 0 THEN ps_export=0
;IF N_Elements(version) EQ 0 THEN version=0
image_filter_fn='filter_uv_uniform' ;applied ONLY to output images
compile_opt strictarr
args = Command_Line_Args(count=nargs)
obs_id = args[0]
;obs_id = 'Pzen.2456265.36872.uvcRREcACOTU.uvfits'
output_directory = args[1]
print,output_directory
;output_directory = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/'
version = args[2]
;version = 'nb_decon_Feb2016'
cmd_args={version:version}

data_directory='/users/jkerriga/data/jkerriga/FHDOutput/'
vis_file_list=file_search(data_directory,obs_id,count=n_files)
;vis_file_list=file_search(data_directory,'Pzen.24563*.*.*COM.uvfits',count=n_files)
print,vis_file_list
;vis_file_list='/users/jkerriga/zen.2456242.29909.uvcRREcACOM.uvfits'
print,version
fhd_file_list=fhd_path_setup(vis_file_list,version=version)
print,fhd_file_list

healpix_path=fhd_path_setup(output_dir=data_directory,subdir='Healpix',output_filename='Combined_obs',version=version)
catalog_file_path=filepath('mwa_calibration_source_list_gleam_kgs_fhd_fornax_pic.sav',root=rootdir('FHD'),subdir='catalog_data')
calibration_catalog_file_path=filepath('mwa_calibration_source_list_gleam_kgs_fhd_fornax_pic.sav',root=rootdir('FHD'),subdir='catalog_data')
model_catalog_file_path=filepath('mwa_calibration_source_list_gleam_kgs_fhd_fornax_pic.sav',root=rootdir('FHD'),subdir='catalog_data')
;catalog_file_path=filepath('mwa_commissioning_source_list_add_FHDaug23deconvolve_fornax_and_VLA_pic.sav',root=rootdir('FHD'),subdir='catalog_data')
;calibration_catalog_file_path=filepath('mwa_commissioning_source_list_add_FHDaug23deconvolve_fornax_and_VLA_pic.sav',root=rootdir('FHD'),subdir='catalog_data')
;model_catalog_file_path=filepath('mwa_commissioning_source_list_add_FHDaug23deconvolve_fornax_and_VLA_pic.sav',root=rootdir('FHD'),subdir='catalog_data')


;ADDED SEP10
;model_visibilities=1
;max_model_sources=20000
;kbinsize=0.1
allow_sidelobe_cal_sources=1
allow_sidelobe_model_sources=1
return_cal_visibilites=1
allow_sidelobe_sources=1
phase_degree=1



;TESTING PARAMETERS
reorder_visibilities=1 ;only use to check calibration scaling
;preserve_weights=1
vis_auto_model=0
;no_phase_calibration=1
;no_complex_beam=1
nfreq_avg=10.
;no_extend=0
;recalculate_all=0
;smooth_width=11.
;absolute_value=1


calibration_auto_fit=0
;FoV=160.
dimension=2048
;preserve_visibilities=1
pad_uv_image=1
beam_offset_time=279.17
residual_flag=0
;cal_gain_init=1E-2
;cal_amp_degree_fit=0
;cal_phase_degree_fit=0
n_pol=2
n_tile=64
n_tile_cut=0
export_images=1
flag_calibration=0
flag_visibilities=0
bandpass_calibrate=0
calibration_flag_iterate=0
snapshot_healpix_export=1
split_ps_export=1
;save_imagecube=1
;save_uvf=1
combine_obs=0
ps_export=0
precess=1 ;At least the new uvfits files from Jonnie need to be precessed to J2000
calibration_auto_initialize=0
instrument='hera'
no_ps=1
;gain_factor=0.15
lat=Ten(-30,42,17.5)
lon=Ten(21,25,41)
time_offset=279.1727364063263;5.*60. ;time offset of phase center from start time. PAPER data are phased to 5 minutes after the start time. 
min_baseline=0.
max_baseline=10000.
min_cal_baseline=1.
calibrate_visibilities=1
max_calibration_sources=1
calibration_polyfit=1 ;says order 1, forced in code to be 0
;no_restrict_cal_sources=1
no_rephase=1
freq_start=100 ;112
freq_end=200 ;182
psf_resolution=64.
beam_residual_threshold=0.001


output_residual_histogram=0
show_beam_contour=0
beam_model_version=2
firstpass=1
mark_zenith=1
;max_cal_iter=10000L
;cal_convergence_threshold=1E-4
no_restrict_cal_sources=1
;antenna_size=2.34
IF N_Elements(extra) GT 0 THEN cmd_args=extra

extra=var_bundle()
general_obs,_Extra=extra

!except=except
END
