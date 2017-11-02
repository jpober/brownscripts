PRO test_psa64,cleanup=cleanup,ps_export=ps_export,recalculate_all=recalculate_all,export_images=export_images,version=version,$
    beam_recalculate=beam_recalculate,healpix_recalculate=healpix_recalculate,mapfn_recalculate=mapfn_recalculate,$
    grid=grid,deconvolve=deconvolve,channel=channel,_Extra=extra
except=!except
!except=0 
heap_gc

;IF N_Elements(recalculate_all) EQ 0 THEN recalculate_all=1
;IF N_Elements(export_images) EQ 0 THEN export_images=1
;IF N_Elements(cleanup) EQ 0 THEN cleanup=0
;IF N_Elements(ps_export) EQ 0 THEN ps_export=0
;IF N_Elements(version) EQ 0 THEN version=0
;IF N_Elements(channel) EQ 0 THEN channel=97
image_filter_fn='filter_uv_radial' ;applied ONLY to output images
args = Command_Line_Args(count=nargs)
obs_id = args[0]
export_images=0
;obs_id = 'Pzen.2456265.36872.uvcRREcACOTU.uvfits'
output_directory = args[1]
print,output_directory
;output_directory = '/nfs/mwa-09/r1/djc/EoR2013/Aug23/'
version = args[2]
;version = 'nb_decon_Feb2016'
cmd_args={version:version}
data_directory='/users/jkerriga/data/jkerriga/2DayOutput/'
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
;noise_calibrate=0
exclude_flagged_sources=0
no_extend=1
deconvolve=1
deconvolution_filter='filter_uv_uniform'
return_sidelobe_catalog=1
lat=Ten(-30,42,17.5)
lon=Ten(21,25,41)
n_tile=64
sigma_threshold=2.
combine_healpix=0
combine_obs=0
FoV=160.
time_offset=279.1727364063263
;align=1
instrument='paper'
flag=0
flag_visibilities=1
no_condense_sources=0
;export_images=1
return_decon_visibilities=1
dimension=2048.
gain_factor=0.15
max_sources=10.
pad_uv_image=2.
precess=1 ;set to 1 ONLY for X16 PXX scans (i.e. Drift_X16.pro)
FoV=160.
no_ps=1 ;don't save postscript copy of images
psf_dim=32

general_obs,cleanup=cleanup,ps_export=ps_export,recalculate_all=recalculate_all,export_images=export_images,version=version,$
    beam_recalculate=beam_recalculate,healpix_recalculate=healpix_recalculate,mapfn_recalculate=mapfn_recalculate,$
    grid=grid,deconvolve=deconvolve,image_filter_fn=image_filter_fn,data_directory=data_directory,$
    vis_file_list=vis_file_list,fhd_file_list=fhd_file_list,healpix_path=healpix_path,catalog_file_path=catalog_file_path,$
    dimension=dimension,max_sources=max_sources,pad_uv_image=pad_uv_image,precess=precess,psf_dim=psf_dim,$
    complex_beam=complex_beam,double_precison_beam=double_precison_beam,FoV=FoV,no_ps=no_ps,_Extra=extra
!except=except
END
