pro quick_n_dirty_image, obsID, freq=freq, pol=pol, fhd_dir=fhd_dir

;; Given an fhd_directory, an obsID, a freq, and a pol, draw a dirty image.

if n_elements(fhd_dir) eq 0 then fhd_dir='.'
if n_elements(freq) eq 0 then freq=50   ; Chan num
if n_elements(pol) eq 0 then pol="XX"

; Load visibility data from vis_data (includes vis_ptr, pol_i, and obs)
restore, fhd_dir+'/vis_data/'+obsID+"_vis_"+pol+".sav", /verbose

; Load psf
restore, fhd_dir+'/beams/'+obsID+"_beams.sav",/verbose

; Load params
restore, fhd_dir+"/metadata/"+obsID+"_params.sav",/verbose

n_param = n_elements(params.uu)

vis_weights = Ptr_new(Replicate(1.,obs.n_freq,n_param))
print, 'Begin gridding'

vis_ptr=ptr_new((*vis_ptr)[freq,*])
;vis_ptr=ptr_new((*vis_ptr)[*,*])

; Below is what's needed to get the dirty_uv_image
dirty_UV=visibility_grid(vis_ptr,vis_weights,obs,status_str,psf,params,file_path_fhd=fhd_dir,polarization=pol_i,timing=timing)
print, "Gridding time: ", timing

; Then apply dirty_image_generate to the uv image.

beam_arr = beam_image(psf, obs,pol_i=pol_i,square=0)

dimension=obs.dimension
elements = obs.elements
horizon_mask=fltarr(dimension,elements)+1.
apply_astrometry, obs, x_arr=meshgrid(dimension,elements,1), y_arr=meshgrid(dimension,elements,2), ra_arr=ra_arr, dec_arr=dec_arr, /xy2ad, /ignore_refraction
horizon_test=where(Finite(ra_arr,/nan),n_horizon_mask)
IF n_horizon_mask GT 0 THEN horizon_mask[horizon_test]=0

beam_base_out = Rebin(beam_arr,dimension,elements)*horizon_mask
beam_i=where(beam_base_out GE 0.01) 
beam_mask=fltarr(dimension,elements) & beam_mask[beam_i]=1.

x_inc=beam_i mod dimension
y_inc=Floor(beam_i/dimension)
zoom_low=min(x_inc)<min(y_inc)
zoom_high=max(x_inc)>max(y_inc)


obs.astr.crpix-=zoom_low
obs.astr.naxis=[zoom_high-zoom_low+1,zoom_high-zoom_low+1]

weights_arr = fltarr(dimension, elements) + 1.

instr_dirty_arr=Ptr_new(dirty_image_generate(dirty_UV,degpix=obs.degpix,weights=weights_arr,/antialias, file_path_fhd=fhd_dir,beam_ptr=beam_base_out))

; Renormalize
normalization_arr=1./(dirty_image_generate(weights_arr,degpix=obs.degpix,obs=obs,psf=psf,params=params,$
    weights=weights_arr,/antialias,beam_ptr=ptr_new(beam_base_out)))[dimension/2.,elements/2.]
normalization_arr*=((beam_base_out)[obs.obsx,obs.obsy])^2.

*instr_dirty_arr *= normalization_arr * 10000

Imagefast,(*instr_dirty_arr)[zoom_low:zoom_high,zoom_low:zoom_high],file_path=fhd_dir+'Dirty_'+obs.pol_names[pol_i]+"_chan-"+strtrim(freq,2),reverse_image=0,/right,sig=2,color_table=0,back='white',title=obsID,/show_grid,astr=obs.astr


END
