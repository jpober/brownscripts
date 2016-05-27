PRO paper_psa64,recalculate_all=recalculate_all,export_images=export_images,version=version,_Extra=extra
except=!except
!except=0 
heap_gc

;IF N_Elements(recalculate_all) EQ 0 THEN recalculate_all=1
;IF N_Elements(export_images) EQ 0 THEN export_images=1
;IF N_Elements(cleanup) EQ 0 THEN cleanup=0
;IF N_Elements(ps_export) EQ 0 THEN ps_export=0
;IF N_Elements(version) EQ 0 THEN version=0
image_filter_fn='filter_uv_uniform' ;applied ONLY to output images

;data_directory=rootdir('mwa')+filepath('',root='PAPER_DATA',subdir=['psa32'])
;data_directory='/data2/PAPER/psa32/'
;IF StrLowCase(!version.os_family) EQ 'unix' THEN data_directory='/data2/PAPER/psa32/' $
;    ELSE data_directory=rootdir('mwa')+filepath('',root='PAPER_DATA',subdir=['psa32_redundant'])
data_directory='/users/jkerriga/'
vis_file_list=file_search(data_directory,'zen.2456242.29909.uvcRREcACOM.uvfits',count=n_files)
;vis_file_list='/users/jkerriga/zen.2456242.29909.uvcRREcACOM.uvfits'
fhd_file_list=fhd_path_setup(vis_file_list,version=version)
print,fhd_file_list
healpix_path=fhd_path_setup(output_dir=data_directory,subdir='Healpix',output_filename='Combined_obs',version=version)
catalog_file_path=filepath('master_sgal_cat.sav',root=rootdir('FHD'),subdir='catalog_data')
calibration_catalog_file_path=filepath('master_sgal_cat.sav',root=rootdir('FHD'),subdir='catalog_data')
model_catalog_file_path=filepath('master_sgal_cat.sav',root=rootdir('FHD'),subdir='catalog_data')

FoV=160.
dimension=1024
beam_offset_time=0
n_pol=2
n_tile=64
n_tile_cut=0
export_images=1
flag_calibration=1
flag_visibilities=0
phase_degree=0
calibration_flag_iterate=0
snapshot_healpix_export=1
split_ps_export=1
combine_obs=0
ps_export=0
precess=1 ;At least the new uvfits files from Jonnie need to be precessed to J2000
calibration_auto_initialize=0
instrument='paper'
no_ps=1
gain_factor=0.15
lat=Ten(-30,42,17.5)
lon=Ten(21,25,41)
time_offset=5.*60. ;time offset of phase center from start time. PAPER data are phased to 5 minutes after the start time. 
min_baseline=1.
min_cal_baseline=30.
calibrate_visibilities=1
max_calibration_sources=10000
calibration_polyfit=0
no_restrict_cal_sources=1
no_rephase=1
freq_start=124
freq_end=174
;beam_cal_threshold=0.5
beam_model_version=2
firstpass=1
antenna_size=2.34
IF N_Elements(extra) GT 0 THEN cmd_args=extra

extra=var_bundle()
general_obs,_Extra=extra

!except=except
END
