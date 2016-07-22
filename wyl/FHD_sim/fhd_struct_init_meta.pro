FUNCTION fhd_struct_init_meta,file_path_vis,hdr,params,lon=lon,lat=lat,alt=alt,n_tile=n_tile,$
    zenra_in=zenra_in,zendec_in=zendec_in,obsra_in=obsra_in,obsdec_in=obsdec_in,phasera_in=phasera_in,phasedec_in=phasedec_in,$
    rephase_to_zenith=rephase_to_zenith,precess=precess,degpix=degpix,dimension=dimension,elements=elements,$
    obsx=obsx,obsy=obsy,instrument=instrument,mirror_X=mirror_X,mirror_Y=mirror_Y,no_rephase=no_rephase,$
    meta_data=meta_data,meta_hdr=meta_hdr,time_offset=time_offset,$
    cotter_precess_fix=cotter_precess_fix,force_rephase_to_zenith=force_rephase_to_zenith,$
    override_target_phasera=override_target_phasera,override_target_phasedec=override_target_phasedec,_Extra=extra

IF N_Elements(instrument) EQ 0 THEN instrument=''
metafits_ext='.metafits'
metafits_dir=file_dirname(file_path_vis)
metafits_name=file_basename(file_path_vis,'.sav',/fold_case)
metafits_name=file_basename(metafits_name,'.uvfits',/fold_case)
metafits_name=file_basename(metafits_name,'_cal',/fold_case) ;sometimes "_cal" is present, sometimes not.
metafits_path=metafits_dir+path_sep()+metafits_name+metafits_ext
Imadethisup_path=metafits_dir+path_sep()+metafits_name+'.uvfits'
print, '###CHECK PATH###'
print, Imadethisup_path
time=params.time
b0i=Uniq(time)
jdate=double(hdr.jd0)+time[b0i]
IF Tag_exist(hdr,'lat') THEN BEGIN
    lat=hdr.lat
    lon=hdr.lon
    alt=hdr.alt
ENDIF ELSE BEGIN
    IF N_Elements(lon) EQ 0 THEN lon=116.67081524 & lon=Float(lon);degrees (MWA, from Tingay et al. 2013)
    IF N_Elements(lat) EQ 0 THEN lat=-26.7033194 & lat=Float(lat);degrees (MWA, from Tingay et al. 2013)
    IF N_Elements(alt) EQ 0 THEN alt=377.827 & alt=Float(alt);altitude (meters) (MWA, from Tingay et al. 2013)
ENDELSE

IF N_Elements(dimension) EQ 0 THEN dimension=1024.
IF N_Elements(elements) EQ 0 THEN elements=dimension
IF N_Elements(obsx) EQ 0 THEN obsx=dimension/2.
IF N_Elements(obsy) EQ 0 THEN obsy=elements/2.

degpix2=[degpix,degpix]
IF Keyword_Set(mirror_X) THEN degpix2[0]*=-1
IF Keyword_Set(mirror_Y) THEN degpix2[1]*=-1
n_pol=hdr.n_pol

IF file_test(metafits_path) THEN BEGIN
    meta_hdr=headfits(metafits_path,exten=0,/silent)
    
    meta_data=mrdfits(metafits_path,1,hdr1,/silent)
    tile_nums=meta_data.antenna
    tile_nums=radix_sort(tile_nums,index=tile_order)
    meta_data=meta_data[tile_order]
    pol_names=meta_data.pol
    single_i=where(pol_names EQ pol_names[0],n_single)
    tile_names=meta_data.tile
    tile_names=tile_names[single_i]
    tile_height=meta_data.height
    tile_height=tile_height[single_i]-alt
    tile_flag=Ptrarr(n_pol) & FOR pol_i=0,n_pol-1 DO tile_flag[pol_i]=Ptr_new(meta_data(single_i+pol_i).flag)
    
    obsra=sxpar(meta_hdr,'RA')
    obsdec=sxpar(meta_hdr,'Dec')
    phasera=sxpar(meta_hdr,'RAPHASE')
    phasedec=sxpar(meta_hdr,'DECPHASE')
    print, '###check ra and dec###'
    print, obsra,obsdec,phasera,phasedec
;    LST=sxpar(meta_hdr,'LST')
;    HA=sxpar(meta_hdr,'HA')
;    HA=ten([Fix(Strmid(HA,0,2)),Fix(Strmid(HA,3,2)),Fix(Strmid(HA,6,2))])*15.
    date_obs=sxpar(meta_hdr,'DATE-OBS')
    IF ~Keyword_Set(time_offset) THEN time_offset=0d
    time_offset/=(24.*3600.)
    JD0=Min(Jdate)+time_offset
    
    epoch=date_conv(JD0,'REAL')/1000.
    epoch_year=Floor(epoch)
    epoch_fraction=(epoch-epoch_year)*1000./365.24218967
    epoch=epoch_year+epoch_fraction    
    
    hor2eq,90.,0.,jd0,zenra,zendec,ha_out,lat=lat,lon=lon,/precess,/nutate    
    
    beamformer_delays=sxpar(meta_hdr,'DELAYS')
    beamformer_delays=Ptr_new(Float(Strsplit(beamformer_delays,',',/extract)))
ENDIF ELSE BEGIN
    ;use hdr and params to guess metadata
    print,'### NOTE ###'
    print,'Metafits file not found! Calculating obs settings from the uvfits header instead'
    
;    print,metafits_path+' not found. Calculating obs settings from the uvfits header instead'
    ;256 tile upper limit is hard-coded in CASA format
    ;these tile numbers have been verified to be correct
    name_mod=2.^((Ceil(Alog(Sqrt(hdr.nbaselines*2.-n_tile))/Alog(2.)))>Floor(Alog(Min(params.baseline_arr))/Alog(2.)))
    tile_A1=Long(Floor(params.baseline_arr/name_mod)) ;tile numbers start from 1
    tile_B1=Long(Fix(params.baseline_arr mod name_mod))
    hist_A1=histogram(tile_A1,min=1,max=n_tile,/binsize,reverse_ind=ria)
    hist_B1=histogram(tile_B1,min=1,max=n_tile,/binsize,reverse_ind=rib)
    hist_AB=hist_A1+hist_B1
    tile_names=indgen(n_tile)+1
    tile_use=where(hist_AB,n_tile_exist,complement=missing_i,ncomplement=missing_n)+1
    tile_height=Fltarr(n_tile)
    tile_flag0=intarr(n_tile)
    IF missing_n GT 0 THEN tile_flag0[missing_i]=1
    tile_flag=Ptrarr(n_pol) & FOR pol_i=0,n_pol-1 DO tile_flag[pol_i]=Ptr_new(tile_flag0)
    date_obs=hdr.date_obs
    
    IF ~Keyword_Set(time_offset) THEN time_offset=0d
    time_offset/=(24.*3600.)
    JD0=Min(Jdate)+time_offset
    epoch=date_conv(JD0,'REAL')/1000.
    epoch_year=Floor(epoch)
    epoch_fraction=(epoch-epoch_year)*1000./365.24218967
    epoch=epoch_year+epoch_fraction   
   
    hdr0=headfits(Imadethisup_path,exten=0,/silent) 
    obsra=hdr.obsra
    obsdec=hdr.obsdec
    print, '### check ra and dec ###'
    print, obsra,obsdec
    IF Keyword_Set(precess) THEN Precess,obsra,obsdec,epoch,2000.
    IF N_Elements(phasera) EQ 0 THEN BEGIN
        phasera=sxpar(hdr0,'RAPHASE',count=cnt)
        IF cnt>0 THEN phasera=phasera ELSE phasera=obsra
    IF N_Elements(phasedec) EQ 0 THEN BEGIN
        phasedec=sxpar(hdr0,'DECPHASE',count=cnt)        
        IF cnt>0 THEN phasedec=phasedec ELSE phasedec=obsdec
    print, phasera, phasedec

    hor2eq,90.,0.,jd0,zenra,zendec,ha_out,lat=lat,lon=lon,/precess,/nutate
    ;beamformer_delays=Ptr_new()
    beamformer_delays=sxpar(hdr0,'DELAYS',count=cnt)
    print, 'cnt=',cnt
    print, beamformer_delays
    IF cnt>0 THEN BEGIN
        beamformer_delays=Ptr_new(Float(Strsplit(beamformer_delays,',',/extract)))
    ENDIF ELSE BEGIN
        beamformer_delays=Ptr_new()
    ENDELSE
ENDELSE

IF N_Elements(zenra_in) EQ 1 THEN zenra=zenra_in
IF N_Elements(zendec_in) EQ 1 THEN zendec=zendec_in
IF N_Elements(obsra_in) EQ 1 THEN obsra=obsra_in
IF N_Elements(obsdec_in) EQ 1 THEN obsdec=obsdec_in
IF N_Elements(phasera_in) EQ 1 THEN phasera=phasera_in
IF N_Elements(phasedec_in) EQ 1 THEN phasedec=phasedec_in
IF N_Elements(override_target_phasera) EQ 0 THEN orig_phasera=phasera ELSE orig_phasera=override_target_phasera
IF N_Elements(override_target_phasedec) EQ 0 THEN orig_phasedec=phasedec ELSE orig_phasedec=override_target_phasedec
IF Keyword_Set(no_rephase) THEN BEGIN
    phasera=obsra
    phasedec=obsdec
ENDIF

;IF Abs(obsra-zenra) LT degpix THEN zenra=obsra
;IF Abs(obsdec-zendec) LT degpix THEN zendec=obsdec

IF Keyword_Set(force_rephase_to_zenith) THEN rephase_to_zenith=1
IF Keyword_Set(rephase_to_zenith) THEN BEGIN
    IF Keyword_Set(force_rephase_to_zenith) THEN BEGIN
        phasera=zenra
        phasedec=zendec
    ENDIF ELSE BEGIN
        phasera=obsra
        phasedec=obsdec
    ENDELSE
    obsra=zenra
    obsdec=zendec
ENDIF
projection_slant_orthographic,astr=astr,degpix=degpix2,obsra=obsra,obsdec=obsdec,zenra=zenra,zendec=zendec,$
    dimension=dimension,elements=elements,obsx=obsx,obsy=obsy,zenx=zenx,zeny=zeny,phasera=phasera,phasedec=phasedec,$
    epoch=2000.,JDate=JD0

Eq2Hor,obsra,obsdec,JD0,obsalt,obsaz,lat=lat,lon=lon,alt=Mean(alt)
meta={obsra:Float(obsra),obsdec:Float(obsdec),zenra:Float(zenra),zendec:Float(zendec),phasera:Float(phasera),phasedec:Float(phasedec),$
    epoch:Float(epoch),tile_names:tile_names,lon:Float(lon),lat:Float(lat),alt:Float(alt),JD0:Double(JD0),Jdate:Double(Jdate),$
    astr:astr,obsx:Float(obsx),obsy:Float(obsy),zenx:Float(zenx),zeny:Float(zeny),obsaz:Float(obsaz),obsalt:Float(obsalt),$
    delays:beamformer_delays,tile_height:Float(tile_height),tile_flag:tile_flag,orig_phasera:Float(orig_phasera),orig_phasedec:Float(orig_phasedec)}

RETURN,meta
END
