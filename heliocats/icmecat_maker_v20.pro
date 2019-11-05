;make_icmecat

;new ICME events - add in ICME master list only the times. 
;edit ICMECAT in master file, save as txt mit tabstopp getrennt
;this program makes the MVA and gets all parameters and saves them as the ICMECAT structure


;-----------------------------------------------------------

common apexdata3hours, apex_r, apex_t, apex_l, apex_s, apex_pa, apex_pa_thick


FUNCTION write_header_comments, asciifilename, icme_obs, icme_master, icmecat, ctype

headerfilename= strmid(asciifilename,0,55)+'header.txt'


openw,lun,headerfilename,/get_lun,width=900
printf,lun,'HEADER FILE FOR: ',asciifilename
printf,lun,' '
printf,lun,'AUTHORS: '
printf,lun,'Christian Moestl, Peter Boakes, University of Graz, Austria and SRI, Austrian Academy of Sciences, Graz, Austria.'
printf,lun,'Alexey Isavnin, Emilia Kilpua, University of Helsinki, Finland.'
printf,lun,'Benoit Lavraud, IRAP, Toulouse, France.'
printf,lun,'We sincerely thank the providers of individual ICME lists: Simon Good, Reka Winslow, Lan Jian and Teresa Nieves-Chinchilla.'
printf,lun,' '
printf,lun,'FILE CREATION DATE: ',systime()
printf,lun,' '
printf,lun,'INPUT FILE:', icme_master
printf,lun,' '
printf,lun,'Number of events in ICMECAT: ', num2str(n_elements(icmecat))
printf, lun, ' ' 
printf, lun,'ICME observatories:', icme_obs
printf, lun, '           '
printf, lun, 'Time range: January 2007 - December 2017.'
printf, lun, '           '
printf, lun, 'Coordinate system for all spacecraft except Ulysses (RTN) and MAVEN (MSO):  ', ctype
printf, lun, '           '
if ctype eq 'HEEQ' then begin 
  printf, lun, 'Heliocentric Earth Equatorial Coordinates (HEEQ): '
  printf, lun, '	Z is the solar rotation axis.'
  printf, lun, '	X points from the Sun to the Earth, at the intersection of the solar equator and solar central meridian as seen from Earth.'
  printf, lun, '	Y completes the right handed triad and points to solar west.'
  printf, lun, '	This system is always centered on Earth, regardless of the spacecraft position.'
endif
  
if ctype eq 'SCEQ' then begin
  printf, lun, 'SpaceCraft Equatorial Coordinates (SCEQ):'
  printf, lun, '	Z is the solar rotation axis.'
  printf, lun, '	X points from the Sun to the spacecraft, projected in the solar equatorial plane. '
  printf, lun, '	Y completes the right handed triad and points to solar west.'
  printf, lun, ' This system is thus centered on the respective in situ spacecraft. '
  printf, lun, ' The solar equatorial plane as the reference plane is similar for all spacecraft.'
endif

printf, lun, '           '
printf, lun,'Notes: '
printf, lun, '(1) For all spacecraft: If there is no sheath region, so the ICME starts immediately with a magnetic obstacle, the ICME_START_TIME is similar to MO_START_TIME.'
printf, lun, '(2) For all spacecraft the times that are given are: the ICME_START_TIME (defined by shock or elevated magnetic field start), the MO_START_TIME,'
printf, lun, '    and the MO_END_TIME (MO for magnetic obstacle, defined by elevated field and smooth magnetic field rotations, but also includes (few) complex ejecta).'
printf, lun, '    Only for Wind the ICME_END_TIME is given extra, which contains interval of disturbed solar wind by the ICME until the values return to slow solar wind.'
printf, lun, '(3) MESSENGER: For ICME events catalogued by Winslow et al. the end time of the ICME is used for MO_END_TIME (not ICME_END_TIME). '
printf, lun, '    Caution: after 2011 March 18 when MESSENGER is in orbit around Mercury the times are sometimes not exact because '
printf, lun, '(4) MESSENGER and VEX: For events cataloged by Simon Good ICME_START_TIME has been added by V. Krupar (Imperial College) and C. Moestl (Graz).'
printf, lun, '(5) For the calculation of the parameters at MESSENGER during the orbit around Mercury, all data points inside the bowshock of Mercury have been removed.'
printf, lun, '    (according to a list thankfully provided to us by by R. Winslow, UNH, B. Anderson, APL, Lydia Philpott, UBC).'
printf, lun, '(6) For MVA to be applied to the MO interval, up to 25% of the MO data points may be NaN. Otherwise MVA results are set to NaN.'
printf, lun, '    This is relevant to MESSENGER in orbit around Mercury and VEX at Venus, where the ICME observation in the solar wind sometimes contains too few data points.'
printf, lun, '(7) Calculation of the MO parameters at VEX is done with an approximate removal of the induced magnetosphere, with a modified equation'
printf, lun, '    according to the one in the discussion of Zhang et al. 2008 (doi: 10.1016/j.pss.2007.09.012), with a constant of 3.5 instead of 2.14/2.364,'
printf, lun, '    in order to account for a larger bowshock distance during solar maximum than studied in this paper. '
printf, lun, '(8) Similarly, for MAVEN data, all data inside the bow shock were removed with the model from Edberg et al. 2008 (JGR), including one standard deviation.'
printf, lun, '    For the remaining data the median for each orbit is taken as 1 data point. 4.5 hour time resolution.'
printf, lun, ''
printf, lun, 'updates in this version 20 compared to version 10: added MAVEN events, extended the Wind catalog to end of 2017 and added dynamic pressure parameters.'
printf,lun,'VARIABLES: ' 

;Magnetic obstacles are defined by elevated magnetic fields over a long duration and can contain ordered or disordered magnetic fields. Categorizations of magnetic obstacles are: flux ropes (FR, elevated, smooth rotating fields), flux-rope-like (FRL, elevated, smooth constant fields), and ejecta (elevated but disordered fields).

;write variable
comments = strarr(34)
comments[0]  = ' 	1: ICMECAT_ID: The unique identifier for the observed ICME. unit: string.'
comments[1]  = ' 	2: SC_INSITU: The name of the in situ observatory. unit: string.'
comments[2]  = ' 	3: ICME_START_TIME: The shock arrival or density enhancement time, can be similar to MO_START_TIME. unit: UTC.'
comments[3]  = ' 	4: MO_START_TIME: The start time of the magnetic obstacle (MO), including flux ropes, flux-rope-like, and ejecta signatures. unit: UTC.'
comments[4]  = ' 	5: MO_END_TIME: The end time of the magnetic obstacle. unit: UTC.'
comments[5]  = ' 	6: ICME_END_TIME: The end time of the ICME, can be similar to MO_END_TIME. unit: UTC.'
comments[6]  = ' 	7: MO_BMAX: The maximum total magnetic field in the magnetic obstacle. unit: nT.'
comments[7]  = ' 	8: MO_BMEAN: The mean total magnetic field of the magnetic obstacle. unit: nT.'
comments[8]  = ' 	9: MO_BSTD: The standard deviation of the total magnetic field of the magnetic obstacle. unit: nT.'
comments[9]  = '	10: MO_BZMEAN: The mean magnetic field Bz component in the magnetic obstacle. unit: nT.'
comments[10] = '	11: MO_BZMIN: The minimum magnetic field Bz component of the magnetic obstacle. unit: nT.'
comments[11] = '	12: MO_DURATION: Duration of interval between MO_START_TIME and MO_END_TIME. unit: hours.'
comments[12] = '	13: SC_HELIODISTANCE: Average heliocentric distance of the spacecraft during the MO. unit: AU.'
comments[13] = '	14: SC_LONG_HEEQ: Average heliospheric longitude of the spacecraft during the MO, range [-180,180]. unit: degree (HEEQ).'
comments[14] = '	15: SC_LAT_HEEQ: Average heliospheric latitude of the spacecraft during the MO, range [-90,90]. unit: degree (HEEQ).'
comments[15] = '	16: MO_MVA_AXIS_LONG: Longitude of axis from minimum variance analysis with magnetic field unit vectors (MVA): X=0 deg, Y=90 deg, range [0,360]. unit: degree ('+ctype+').'
comments[16] = '	17: MO_MVA_AXIS_LAT: Latitude of axis from MVA, +Z=-90 deg, -Z=-90, range [-90,90]. unit: degree ('+ctype+').'
comments[17] = '	18: MO_MVA_RATIO: Eigenvalue 2 over 3 ratio as indicator of reliability of MVA, must be > 2, otherwise NaN. unit: number.'
comments[18] = '	19: SHEATH_SPEED: For STEREO-A/B, Wind, MAVEN: average proton speed from ICME_START_TIME to MO_START_TIME, NaN if these times are similar. unit: km/s.'
comments[19] = '	20: SHEATH_SPEED_STD: For STEREO-A/B, Wind, MAVEN: standard deviation of proton speed from ICME_START_TIME to MO_START_TIME, NaN if these times are similar. unit: km/s.'
comments[20] = '	21: MO_SPEED: For STEREO-A/B, Wind, MAVEN: average proton speed from MO_START_TIME to MO_END_TIME. unit: km/s.'
comments[21] = '	22: MO_SPEED_STD: For STEREO-A/B, Wind, MAVEN: standard deviation of proton speed from MO_START_TIME to MO_END_TIME. unit: km/s.'
comments[22] = '	23: SHEATH_DENSITY: For STEREO-A/B, Wind, MAVEN: average proton density from ICME_START_TIME to MO_START_TIME, NaN if these times are similar. unit: ccm^-3.'
comments[23] = '	24: SHEATH_DENSITY_STD: For STEREO-A/B, Wind, MAVEN: standard deviation of proton density from ICME_START_TIME to MO_START_TIME, NaN if these times are similar. unit: cm^-3.'
comments[24] = '	25: MO_DENSITY: For STEREO-A/B, Wind, MAVEN: average proton density from MO_START_TIME to MO_END_TIME. unit: cm^-3.'
comments[25] = '	26: MO_DENSITY_STD: For STEREO-A/B, Wind, MAVEN: standard deviation of proton density from MO_START_TIME to MO_END_TIME. unit: cm^-3.'
comments[26] = '	27: SHEATH_TEMPERATURE: For STEREO-A/B, Wind, MAVEN:average proton temperature from ICME_START_TIME to MO_START_TIME, NaN if these times are similar. unit: K.'
comments[27] = '	28: SHEATH_TEMPERATURE_STD: For STEREO-A/B, Wind, MAVEN: standard deviation of proton temperature from ICME_START_TIME to MO_START_TIME, NaN if these times are similar. unit: K.'
comments[28] = '	29: MO_TEMPERATURE: For STEREO-A/B, Wind, MAVEN: average proton temperature from MO_START_TIME to MO_END_TIME. unit: K.'
comments[29] = '	30: SHEATH_PDYN_MEAN: For STEREO-A/B, Wind, MAVEN: mean of dynamic pressure assuming only protons contribute, from ICME_START_TIME to MO_START_TIME. unit: nPa.'
comments[30] = '	31: SHEATH_PDYN_STD: For STEREO-A/B, Wind, MAVEN: standard deviation of dynamic pressure assuming only protons contribute, from ICME_START_TIME to MO_START_TIME. unit: nPa.'
comments[31] = '	32: MO_PDYN_MEAN: For STEREO-A/B, Wind, MAVEN: mean of dynamic pressure assuming only protons contribute, from MO_START_TIME to MO_END_TIME. unit: nPa.'
comments[32] = '	33: MO_PDYN_STD: For STEREO-A/B, Wind, MAVEN: standard deviation of dynamic pressure assuming only protons contribute, from MO_START_TIME to MO_END_TIME. unit: nPa.'


for i=0,n_elements(comments)-1 do begin
	printf,lun,comments[i]
endfor


free_lun, lun

return, comments
	

END



FUNCTION dynamic_pressure, pdensity, pspeed

   ;pdyn=np.zeros(len([density])) #in nano Pascals
   protonmass=1.6726219*1e-27  ;kg
   ;assume pdyn is only due to protons
   dynamicpressure=((pspeed*1e3)^2)*pdensity*1e6*protonmass*1e9
   return, dynamicpressure ;in nanoPascal

END



;_--------------------------------------- MINIMUM VARIANCE ANALYSIS

FUNCTION mvub, bx, by, bz, time

print, 'START MVUB analysis'


;linearly interpolate NaNs
gooddata = where(Finite(bx), ngooddata, comp=baddata, ncomp=nbaddata)
ratiox=float(nbaddata)/n_elements(bx)
if nbaddata gt 0 then bx[baddata] = interpol(bx[gooddata], gooddata, baddata) 

gooddata = where(Finite(by), ngooddata, comp=baddata, ncomp=nbaddata)
ratioy=float(nbaddata)/n_elements(by)
if nbaddata gt 0 then by[baddata] = interpol(by[gooddata], gooddata, baddata) 

gooddata = where(Finite(bz), ngooddata, comp=baddata, ncomp=nbaddata)
ratioz=float(nbaddata)/n_elements(bz)
if nbaddata gt 0 then bz[baddata] = interpol(bz[gooddata], gooddata, baddata) 

ratiobadall=(ratiox+ratioy+ratioz)/3

print, 'Ratio bad data to all data	 in percent',ratiobadall*100



B=sqrt(bx^2+by^2+bz^2);

Bhat=dblarr(n_elements(bx),3)
Bhat(*,0)=bx/B
Bhat(*,1)=by/B
Bhat(*,2)=bz/B

;covariance Matrix
M=dblarr(3,3)
for i=0,2 do begin
 for j=0,2 do begin 
  M(i,j)=mean((Bhat(*,i)*Bhat(*,j)))-mean(Bhat(*,i))*mean(Bhat(*,j)); 
 endfor
endfor 

;in M no NaNs are allowed
if total(finite(M)) gt 0 then begin 


 eval = EIGENQL(M, EIGENVECTORS = evec, RESIDUAL = residual, /double)
 ;Print the eigenvalues and eigenvectors:
 PRINT, 'Eigenvalues: '
 PRINT, eval
 PRINT, 'Eigenvectors: ' ;in rows
 PRINT, evec
 ; The ith row of the returned array contains the eigenvector corresponding to the ith eigenvalue. 

 pol_vec=evec(*,0) ;max variance
 axis_vec=evec(*, 1);Zeile 2 ist die intermediate var direction
 rad_vec=evec(*,2) ;min variance


 ;now axis should point in the correct direction, not in the opposite

 ;project components to MVA system

 bxstar=dblarr(n_elements(bx))
 bystar=dblarr(n_elements(bx))
 bzstar=dblarr(n_elements(bx))

 ;find components in MVA system by dot products
 for i=0, n_elements(bx)-1 do begin
  bxstar(i)= Bhat(i,*) # pol_vec  ;bxstar is now the poloidal cloud component
  bystar(i)= Bhat(i,*) # axis_vec ;bystar is the axial component
  bzstar(i)= Bhat(i,*) # rad_vec; bzstar is the radial component
 endfor

;plot hodogram
;screen_prefs
;window, 5, xsize=500,ysize=500, retain=2, xpos=500, ypos=00      
;plot,  bystar
;oplot,  bystar
;oplot, bzstar
;window, 2, xsize=1500,ysize=750, retain=2, xpos=100, ypos=500      

 ;the axis has to point in the positive direction in the MVA system, so bystar needs to be positive
 if mean(bystar) lt 0 then begin 
   axis_vec=-axis_vec
   pol_vec=-pol_vec
   rad_vec=-rad_vec
 endif  


 ;criterion after Lepping and Behannon
 eigenvalue_ratio_2_3=eval(1)/eval(2)
 if eigenvalue_ratio_2_3 lt 2 then begin 
   axis_longitude=!VALUES.F_NAN
   axis_latitude=!VALUES.F_NAN
 endif

 if eigenvalue_ratio_2_3 ge 2 then begin 
 ;convert intermediate variance direction to angles 
  axis_latitude=round(asin(axis_vec(2))*180/3.14159265)
  axis_longitude=round(acos(axis_vec(0)/sqrt(axis_vec(0)^2+axis_vec(1)^2))*180/3.14159265)
  if axis_vec(1) lt 0 then axis_longitude=360-axis_longitude;
 endif

endif else begin 

   axis_longitude=!VALUES.F_NAN
   axis_latitude=!VALUES.F_NAN
   eigenvalue_ratio_2_3=!VALUES.F_NAN
   
endelse


;criterium: if more 25 percent bad data, set results to NaN
if ratiobadall gt 0.25 then begin

   axis_longitude=!VALUES.F_NAN
   axis_latitude=!VALUES.F_NAN
   eigenvalue_ratio_2_3=!VALUES.F_NAN
   
end   


print, 'MVUB results'
print, axis_longitude
print, axis_latitude
print, eigenvalue_ratio_2_3

print, 'END MVUB analysis'


	return, [axis_longitude, axis_latitude, eigenvalue_ratio_2_3]

END







FUNCTION make_plot, tp, bp, bxp, byp, bzp, range, start, ende, spacecraft, id, lid, mvalong, mvalat, mvaratio, ctype, icme_start, icme_ende


    common apexdata3hours
  
     
     	;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;left ------------------------- plot spacecraft positions and HI circles
     	
      !P.POSITION=[0.05, 0.1, 0.45, 0.9]	 ;
	     
	     ;als time einen frame aus dem circle move nehmen der am nächsten der icme_start time ist
	     ;find closest time in all_apex_t
	     closetime=min(abs(apex_t-icme_start),closetimeindex)
	     
     
      ;spacecraft is taken by ssc_plot_where_chris to mark the active spacecraft with a square
     	if apex_t(closetimeindex) lt anytim('18-Mar-2011') then begin
  		   ;ssc_plot_where_chris, anytim(apex_t(closetimeindex), /vms), /WHITE_BG, /inner, /MESSENGER, /ULYSSES, /MAVEN, /Rosetta, charsize=5, marksc=spacecraft, scchar=1.5, xthick=2,ythick=2
  		   ssc_plot_where_chris, anytim(icme_start, /vms), /WHITE_BG, /inner, /MESSENGER, /ULYSSES, /MAVEN, /Rosetta, charsize=5, marksc=spacecraft, scchar=1.5, xthick=2,ythick=2    
    		 endif else begin 
      	 if spacecraft eq 'MESSENGER' then spacecraft='Mercury' ;show Mercury position after orbit insertion - no /MESSENGER because it would overplot the text
     		  ;ssc_plot_where_chris, anytim(apex_t(closetimeindex),/vms), /WHITE_BG, /inner, /ULYSSES, /MAVEN, /Rosetta, charsize=3, marksc=spacecraft, scchar=1.5, xthick=2,ythick=2
  	  		  ssc_plot_where_chris, anytim(icme_start,/vms), /WHITE_BG, /inner, /ULYSSES, /MAVEN, /Rosetta, charsize=3, marksc=spacecraft, scchar=1.5, xthick=2,ythick=2
    	 endelse
    	 
    	 
      check_on_linkcat=where(lid eq id)   
      
      ;plot circles for event in active frame
            
        lambda=30
    	 ;time
      timeindicator=anytim(apex_t(closetimeindex), /ccsds)
      timeindicator_plot=anytim(icme_start, /ccsds)
   
   
		     													
									sun=[0,0]
			
									r_earth=get_stereo_lonlat(anytim(timeindicator), 'Earth', system='HEEQ')
									earth=[0,r_earth/149597870]
									PLOTS, [sun(0), earth(0)], [sun(1), earth(1)], color=0, linestyle=0, /DATA, thick=1 
			
									;get indices of CMEs active at the current time
									active_cme=where(apex_t eq anytim(timeindicator))
			
									;when CME apex is found for this time
									if active_cme(0) ne -1 then begin 
					
										;draw all active CMEs as circle segments
											for p=0, n_elements(active_cme)-1 do begin
												;print, p
												;print, 'plotting apex and longitude:', apex_r(active_cme(p)), apex_l(active_cme(p))
												dir=[sin(apex_l(active_cme(p))*!dtor),cos(apex_l(active_cme(p))*!dtor)]*apex_r(active_cme(p)) 	    
												;print, dir
												;plot symbol to apex
												;PLOTS, [sun(0)+dir(0)], [sun(1)+dir(1)], psym=8,color=255, /DATA, symsize=1
							
												;plot circles      
												;make points of circles for half circle -110 to +110 deg away of apex direction
												;angle for circle is in radians
												circ_ang = (findgen(111)*2-20)*!dtor-apex_l(active_cme(p))*!dtor  
												;use this line for -120 to 120 around apex
												circ_ang = (findgen(121)*2-30)*!dtor-apex_l(active_cme(p))*!dtor      
												
												;if apex_r(active_cme(p)) lt 2.5 then begin
  										xc = sun(0)+dir(0)/(1+sin(lambda*!dtor)) + (apex_r(active_cme(p))*sin(lambda*!dtor)/(1+sin(lambda*!dtor)))*cos(circ_ang)
		  								yc = sun(1)+dir(1)/(1+sin(lambda*!dtor)) + (apex_r(active_cme(p))*sin(lambda*!dtor)/(1+sin(lambda*!dtor)))*sin(circ_ang)
								
												
											 loadct, 13, /silent	
		
												if apex_s(active_cme(p)) eq 'A' then begin
																PLOTS, xc,yc, /data, linestyle=0, thick=apex_pa_thick(active_cme(p)), color=255, noclip=0      

																endif
												if apex_s(active_cme(p)) eq 'B' then begin 
																PLOTS, xc,yc, /data, linestyle=0, thick=apex_pa_thick(active_cme(p)), color=50, noclip=0    
																endif
												;endif				

										endfor   	    	  	   			  
									endif
       ;endif 
    	  
 	    
    	 
	 	
	     loadct, 0, /silent
	      
      xyouts, 0.35, 0.035, strmid(timeindicator_plot,0,16), /normal, charsize=2.5, charthick=2, font=1, alignment=0.5, color=0
	    
	 	
	 	 
	 
      ;right ---------------------------------------------------plot in situ data 
      
      !P.POSITION=[0.55, 0.1, 0.95, 0.9]	 
     
     
      coord_indicator=ctype
  
      loadct, 13, /silent          
      ;plotten MO 
      utplot, tp,bp, color=0, yrange=[-range-30,range+30], timerange=[start-24*60*60.,ende+24*60*60.], $
 	     ytitle='B [nT] '+coord_indicator, xstyle=1, charsize=2, charthick=2, font=1, minors=2, thick=2,$
	      ystyle=1, yticklen=0.03, xticklen=0.03,  /NOERASE, ANYTIM_LABEL_EXTRA = {ccsds:1}
	     ;this is the ICME boundary
	    	outplot,[icme_start, icme_start], [-range-30,range+30], color=0, linestyle=2,thick=2, /DATA
	    	outplot,[icme_ende, icme_ende], [-range-30,range+30], color=0, linestyle=2,thick=2, /DATA
	       	     
	     
	     ;this is the MO boundary
	    	outplot,[start, start], [-range-30,range+30], color=0, thick=2, /DATA
	    	outplot,[ende, ende], [-range-30,range+30], color=0, thick=2, /DATA
	    	
	    	
	    	;0 line	    	
   	 	outplot,[start, ende], [0,0], color=0, thick=2, /DATA, linestyle=2
	      
      ;title   	     
      xyouts, 0.5, 0.95, id, /normal, charsize=5, charthick=3, font=1, alignment=0.5, color=0
      
      ;plot components
	     outplot, tp,bxp, color=250, thick=2
	     outplot, tp,byp, color=200, thick=2
	     outplot, tp,bzp, color=50, thick=2
	     
	    
     
      ;legend
      xyouts,0.8,0.92, /normal, 'Bx',charsize=2, charthick=2, /data, font=1, color=255
      xyouts,0.8+0.05,0.92, /normal, 'By',charsize=2, charthick=2, /data, font=1, color=200
      xyouts,0.8+0.1,0.92, /normal, 'Bz',charsize=2, charthick=2, /data, font=1, color=50
      
         ;legend
      xyouts,0.89,0.85, /normal, 'north',charsize=2, charthick=2, /data, font=1, color=50
      xyouts,0.89,0.82, /normal, 'west',charsize=2, charthick=2, /data, font=1, color=200
      
      xyouts,0.89,0.18, /normal, 'south',charsize=2, charthick=2, /data, font=1, color=50
      xyouts,0.89,0.15, /normal, 'east',charsize=2, charthick=2, /data, font=1, color=200
      
           
      ;MVA results
      if finite(mvalong) then longstring=['longitude: '+ string(mvalong, format='(I5)')] else longstring=['longitude: -']
      if finite(mvalat) then latstring=['latitude: '+ string(mvalat, format='(I5)')] else latstring=['latitude: -']
      if finite(mvaratio) then ratiostring=['ratio: '+ string(mvaratio, format='(F5.1)')] else ratiostring=['ratio: -']
      xyouts,0.58+0.0,0.30, /normal, 'MVA',charsize=2, charthick=2, /data, font=1, color=0
      xyouts,0.58+0.0,0.25, /normal, longstring,charsize=2, charthick=2, /data, font=1, color=0
      xyouts,0.58+0.0,0.20, /normal,latstring,charsize=2, charthick=2, /data, font=1, color=0
      xyouts,0.58+0.0,0.15, /normal, ratiostring,charsize=2, charthick=2, /data, font=1, color=0

	          	     
	          	     
	          	     
	          	     
      filejpg='/home/cmoestl/helcats/products/ICMECAT/plots/events_jpg_'+ctype+'/'+id+'_'+ctype+'.jpg'
      x2jpeg, filejpg
      
      ;checken ob das event im linkcat ist - extra unter /plots/linkcat/events/ speichern
      if  check_on_linkcat(0) gt 0 then begin
       print, 'event is in linkcat'
       filejpg='/home/cmoestl/helcats/products/LINKCAT/plots/events_jpg_'+ctype+'/'+id+'_'+ctype+'.jpg'
       x2jpeg, filejpg
      endif 
















     ;--------------------------------------------------------------do the same as eps
   
   

     ;plot this twice if the event is in linkcat
     if  check_on_linkcat(0) gt 0 then  fileeps=['/home/cmoestl/helcats/products/ICMECAT/plots/events_eps_'+ctype+'/'+id+'_'+ctype+'.eps',$
                                                 '/home/cmoestl/helcats/products/LINKCAT/plots/events_eps_'+ctype+'/'+id+'_'+ctype+'.eps' ]
     if  check_on_linkcat(0) lt 0 then   fileeps='/home/cmoestl/helcats/products/ICMECAT/plots/events_eps_'+ctype+'/'+id+'_'+ctype+'.eps'


     for epsplots=0, n_elements(fileeps)-1 do begin 
     
   
      set_plot, 'ps'

 
      device, /encapsulated, filename=fileeps[epsplots], xsize=20, ysize=10, yoffset=15, /inches, /color, bits_per_pixel=8, portrait=1


     	;left ------------------------- plot spacecraft positions and HI circles
     	
      !P.POSITION=[0.05, 0.1, 0.45, 0.9]	 ;
	     
	     ;als time einen frame aus dem circle move nehmen der am nächsten ist
	     ;find closest time in all_apex_t
	     closetime=min(abs(apex_t-icme_start),closetimeindex)
	     
     
       ;spacecraft is taken by ssc_plot_where_chris to mark the active spacecraft with a square
     	if apex_t(closetimeindex) lt anytim('18-Mar-2011') then begin
  		   ;ssc_plot_where_chris, anytim(apex_t(closetimeindex), /vms), /WHITE_BG, /inner, /MESSENGER, /ULYSSES, /MAVEN, /Rosetta, charsize=5, marksc=spacecraft, scchar=1.5, xthick=2,ythick=2
  		   ssc_plot_where_chris, anytim(icme_start, /vms), /WHITE_BG, /inner, /MESSENGER, /ULYSSES, /MAVEN, /Rosetta, charsize=5, marksc=spacecraft, scchar=1.5, xthick=2,ythick=2    
    		 endif else begin 
      	 if spacecraft eq 'MESSENGER' then spacecraft='Mercury' ;show Mercury position after orbit insertion - no /MESSENGER because it would overplot the text
     		  ;ssc_plot_where_chris, anytim(apex_t(closetimeindex),/vms), /WHITE_BG, /inner, /ULYSSES, /MAVEN, /Rosetta, charsize=3, marksc=spacecraft, scchar=1.5, xthick=2,ythick=2
  	  		  ssc_plot_where_chris, anytim(icme_start,/vms), /WHITE_BG, /inner, /ULYSSES, /MAVEN, /Rosetta, charsize=3, marksc=spacecraft, scchar=1.5, xthick=2,ythick=2
    	 endelse
    	 
      ;time
      xyouts, 0.35, 0.035, strmid(timeindicator_plot,0,16), /normal, charsize=2.5, charthick=2, font=1, alignment=0.5
	         	 
      check_on_linkcat=where(lid eq id)   
      
      ;plot circles for event in active frame
      
      ;if  check_on_linkcat(0) gt 0 then begin
      
        lambda=30
  
		 		
													
									sun=[0,0]
			
									r_earth=get_stereo_lonlat(anytim(timeindicator), 'Earth', system='HEEQ')
									earth=[0,r_earth/149597870]
									PLOTS, [sun(0), earth(0)], [sun(1), earth(1)], color=0, linestyle=0, /DATA, thick=1 
			
									;get indices of CMEs active at the current time
									active_cme=where(apex_t eq anytim(timeindicator))
			
									;when CME apex is found for this time
									if active_cme(0) ne -1 then begin 
					
										;draw all active CMEs as circle segments
											for p=0, n_elements(active_cme)-1 do begin
												;print, p
												;print, 'plotting apex and longitude:', apex_r(active_cme(p)), apex_l(active_cme(p))
												dir=[sin(apex_l(active_cme(p))*!dtor),cos(apex_l(active_cme(p))*!dtor)]*apex_r(active_cme(p)) 	    
												;print, dir
												;plot symbol to apex
												;PLOTS, [sun(0)+dir(0)], [sun(1)+dir(1)], psym=8,color=255, /DATA, symsize=1
							
												;plot circles      
												;make points of circles for half circle -110 to +110 deg away of apex direction
												;angle for circle is in radians
												circ_ang = (findgen(111)*2-20)*!dtor-apex_l(active_cme(p))*!dtor  
												;use this line for -120 to 120 around apex
												circ_ang = (findgen(121)*2-30)*!dtor-apex_l(active_cme(p))*!dtor      
												
												;if apex_r(active_cme(p)) lt 2.5 then begin
  										xc = sun(0)+dir(0)/(1+sin(lambda*!dtor)) + (apex_r(active_cme(p))*sin(lambda*!dtor)/(1+sin(lambda*!dtor)))*cos(circ_ang)
		  								yc = sun(1)+dir(1)/(1+sin(lambda*!dtor)) + (apex_r(active_cme(p))*sin(lambda*!dtor)/(1+sin(lambda*!dtor)))*sin(circ_ang)
								
												loadct, 13, /silent	
												if apex_s(active_cme(p)) eq 'A' then begin
																PLOTS, xc,yc, /data, linestyle=0, thick=apex_pa_thick(active_cme(p)), color=255, noclip=0      

																endif
												if apex_s(active_cme(p)) eq 'B' then begin 
																PLOTS, xc,yc, /data, linestyle=0, thick=apex_pa_thick(active_cme(p)), color=50, noclip=0    
																endif
												;endif				

										endfor   	    	  	   			  
									endif
       ;endif 
    	  
	 	 
	 
      ;right ---------------------------------------------------plot in situ data 
      
      !P.POSITION=[0.55, 0.1, 0.95, 0.9]	 
     
     
      coord_indicator=ctype


      loadct, 13, /silent
      ;plotten MO 
      utplot, tp,bp, color=0, yrange=[-range-30,range+30], timerange=[start-24*60*60.,ende+24*60*60.], $
 	     ytitle='B [nT] '+coord_indicator, xstyle=1, charsize=2, charthick=2, font=1, minors=2, thick=2,$
	      ystyle=1, yticklen=0.03, xticklen=0.03,  /NOERASE,  ANYTIM_LABEL_EXTRA = {ccsds:1},xthick=4,ythick=4
	     
	     
	     ;this is the ICME boundary
	    	outplot,[icme_start, icme_start], [-range-30,range+30], color=0, linestyle=2,thick=2, /DATA
	    	outplot,[icme_ende, icme_ende], [-range-30,range+30], color=0, linestyle=2,thick=2, /DATA
	    	
	     	     
	     ;this is the MO boundary
	    	outplot,[start, start], [-range-30,range+30], color=0, thick=2, /DATA
	    	outplot,[ende, ende], [-range-30,range+30], color=0, thick=2, /DATA
	    	
	    	
	    	;0 line	    	
   	 	outplot,[start, ende], [0,0], color=0, thick=2, /DATA, linestyle=2
	      
      ;title   	     
      xyouts, 0.5, 0.95, id, /normal, charsize=5, charthick=3, font=1, alignment=0.5
      
      ;plot components
	     outplot, tp,bxp, color=250, thick=2
	     outplot, tp,byp, color=200, thick=2
	     outplot, tp,bzp, color=50, thick=2
	     
	    
     
      ;legend
      xyouts,0.8,0.92, /normal, 'Bx',charsize=2, charthick=2, /data, font=1, color=255
      xyouts,0.8+0.05,0.92, /normal, 'By',charsize=2, charthick=2, /data, font=1, color=200
      xyouts,0.8+0.1,0.92, /normal, 'Bz',charsize=2, charthick=2, /data, font=1, color=50
      
        ;legend
      xyouts,0.89,0.85, /normal, 'north',charsize=2, charthick=2, /data, font=1, color=50
      xyouts,0.89,0.82, /normal, 'west',charsize=2, charthick=2, /data, font=1, color=200
      
      xyouts,0.89,0.18, /normal, 'south',charsize=2, charthick=2, /data, font=1, color=50
      xyouts,0.89,0.15, /normal, 'east',charsize=2, charthick=2, /data, font=1, color=200
      
      
      
      ;MVA results
      if finite(mvalong) then longstring=['longitude: '+ string(mvalong, format='(I5)')] else longstring=['longitude: -']
      if finite(mvalat) then latstring=['latitude: '+ string(mvalat, format='(I5)')] else latstring=['latitude: -']
      if finite(mvaratio) then ratiostring=['ratio: '+ string(mvaratio, format='(F5.1)')] else ratiostring=['ratio: -']
      xyouts,0.58+0.0,0.30, /normal, 'MVA',charsize=2, charthick=2, /data, font=1, color=0
      xyouts,0.58+0.0,0.25, /normal, longstring,charsize=2, charthick=2, /data, font=1, color=0
      xyouts,0.58+0.0,0.20, /normal,latstring,charsize=2, charthick=2, /data, font=1, color=0
      xyouts,0.58+0.0,0.15, /normal, ratiostring,charsize=2, charthick=2, /data, font=1, color=0
           
      device, /close
      
      set_plot, 'x'     
    endfor     
 
 return, 0
 
 end
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 





;----------------------------------- MAIN


;make plots if gt 0
make_plots_on=0

;choose HEEQ or SCEQ

;ctype='HEEQ' 
ctype='SCEQ' 



;get HIKINCAT for plots (circles and spacecraft positions) 
restore, '/home/cmoestl/helcats/products/ARRCAT/HIKINCAT/apex_from_hicat_final_3hours.sav'


;---------------- read shock lists first 


;restore, filename= '/home/cmoestl/helcats/products/ICMECAT/shocklists/ips_mes_vex.sav', /verbose

;print, anytim(mercury_ips_anytime, /ccsd)
;print, anytim(venus_ips_anytime, /ccsd)


;restore, filename= '/home/cmoestl/helcats/products/ICMECAT/shocklists/krupar_shocks.sav', /verbose




set_plot, 'x'



;------------------------------------------------------ read in MASTER TABLE

;make new icmecat structure from EXCEL master table (manually editable, on mac)

;write EXCEL master file into txt "tabstoppgetrennter text", read this in as ascii and make icmecat structure
masterfile='/home/cmoestl/helcats/products/ICMECAT/HELCATS_ICMECAT_v20_master.txt'
;read ascii
print, 'reading in ICMECAT MASTER'

;define template first time
;template_icmemaster=ascii_template(masterfile)
;save, template_icmemaster, filename= '/home/cmoestl/helcats/data/icmelists.sav'

;restore template for reading in ICMECAT master
restore, filename= '/home/cmoestl/helcats/data/icmelists.sav'

icmecat=read_ascii(masterfile, template=template_icmemaster, /verbose)

;convert to array of structures
icmecat = convert_structure_of_arrays_format(icmecat)

help, icmecat, /st


;Which observatories they contain - needed for header files
icme_observatories=['Wind','STEREO-A', 'STEREO-B', 'VEX', 'MESSENGER','ULYSSES', 'MAVEN' ]



; ----------------------------------- check for LINKCAT overlaps

;get linkcat
restore, filename='/home/cmoestl/helcats/products/LINKCAT/HELCATS_LINKCAT_v10.sav', /verbose

;check which ICME events overlap with linkcat

;to be used by
lid=linkcat.icmecat_id;
    


;------------------------------------------------ add parameters from data





AU=149598000


;define master parameters

;magnetic field in MO
mbmax=dblarr(n_elements(icmecat))
mbmean=dblarr(n_elements(icmecat))
mbstd=dblarr(n_elements(icmecat))
mbzmin=dblarr(n_elements(icmecat))
mbzmean=dblarr(n_elements(icmecat))
mbzstd=dblarr(n_elements(icmecat))




;magnetic field in sheath
;mbmax_s=dblarr(n_elements(icmecat))
;mbmean_s=dblarr(n_elements(icmecat))
;mbstd_s=dblarr(n_elements(icmecat))
;mbzmin_s=dblarr(n_elements(icmecat))
;mbzmean_s=dblarr(n_elements(icmecat))
;mbzstd_s=dblarr(n_elements(icmecat))


;duration and position
mmo_duration=dblarr(n_elements(icmecat))
mheliodistance=dblarr(n_elements(icmecat))
mlongitude=dblarr(n_elements(icmecat))
mlatitude=dblarr(n_elements(icmecat))

;MVA results
maxis_longitude=dblarr(n_elements(icmecat))
maxis_latitude=dblarr(n_elements(icmecat))
meigenvalue_ratio_2_3=dblarr(n_elements(icmecat))

;speed results
mmospeed=dblarr(n_elements(icmecat))
mmospeedstd=dblarr(n_elements(icmecat))
msheathspeed=dblarr(n_elements(icmecat))
msheathspeedstd=dblarr(n_elements(icmecat))

;density, temperature
mmodensity=dblarr(n_elements(icmecat))
mmodensitystd=dblarr(n_elements(icmecat))
mmotemperature=dblarr(n_elements(icmecat))
mmotemperaturestd=dblarr(n_elements(icmecat))

msheathdensity=dblarr(n_elements(icmecat))
msheathdensitystd=dblarr(n_elements(icmecat))
msheathtemperature=dblarr(n_elements(icmecat))
msheathtemperaturestd=dblarr(n_elements(icmecat))



msheathpdyn=dblarr(n_elements(icmecat))
msheathpdynstd=dblarr(n_elements(icmecat))
mmopdyn=dblarr(n_elements(icmecat))
mmopdynstd=dblarr(n_elements(icmecat))


;------------------------------------------------------------------ ULYSSES
    
;get ICME parameters (or rather MO)  
ulyicmes=where(icmecat.sc_insitu eq 'ULYSSES')
start_uly=anytim(icmecat(ulyicmes).mo_start_time)
end_uly=anytim(icmecat(ulyicmes).mo_end_time)

;----------------for ICME start and end times
icme_start_uly=anytim(icmecat(ulyicmes).icme_start_time)
icme_end_uly=anytim(icmecat(ulyicmes).icme_end_time)


;define parameters
bmax=dblarr(n_elements(ulyicmes))
bmean=dblarr(n_elements(ulyicmes))
bstd=dblarr(n_elements(ulyicmes))
bzmean=dblarr(n_elements(ulyicmes))
bzmin=dblarr(n_elements(ulyicmes))



duration=dblarr(n_elements(ulyicmes))
heliodistance=dblarr(n_elements(ulyicmes))
longitude=dblarr(n_elements(ulyicmes))
latitude=dblarr(n_elements(ulyicmes))

axis_longitude=dblarr(n_elements(ulyicmes))
axis_latitude=dblarr(n_elements(ulyicmes))
eigenvalue_ratio_2_3=dblarr(n_elements(ulyicmes))

mospeed=dblarr(n_elements(ulyicmes))
mospeedstd=dblarr(n_elements(ulyicmes))
sheathspeed=dblarr(n_elements(ulyicmes))
sheathspeedstd=dblarr(n_elements(ulyicmes))

modensity=dblarr(n_elements(ulyicmes))
modensitystd=dblarr(n_elements(ulyicmes))
motemperature=dblarr(n_elements(ulyicmes))
motemperaturestd=dblarr(n_elements(ulyicmes))

sheathdensity=dblarr(n_elements(ulyicmes))
sheathdensitystd=dblarr(n_elements(ulyicmes))
sheathtemperature=dblarr(n_elements(ulyicmes))
sheathtemperaturestd=dblarr(n_elements(ulyicmes))


       

savefilename='/home/cmoestl/helcats/products/DATACAT/finaldata/ULY_2006to2009.sav'
restore, savefilename, /verbose
	
  ;get data only for 2007 because the other years are not near ecliptic
  helior=uly_2007.ULY_RADIUS_IN_KM_HEEQ
  longit=uly_2007.ULY_LONGITUDE_IN_RADIANS_HEEQ
  latit=uly_2007.ULY_LATITUDE_IN_RADIANS_HEEQ
  bx=uly_2007.br		
  by=uly_2007.bt
  bz=uly_2007.bn
  btot=uly_2007.btot
  timeb=uly_2007.time
  vtot=uly_2007.v
  den=uly_2007.den
  temp=uly_2007.tp1
    
    
   for p=0,n_elements(ulyicmes)-1 do begin
   
     ;speeds in sheath
     time_inside_sheath=where(timeb gt icme_start_uly(p) and timeb lt start_uly(p))       
     if time_inside_sheath(0) gt -1 then begin     
       sheathspeed(p)=mean(vtot(time_inside_sheath), /NAN)
       sheathspeedstd(p)=stddev(vtot(time_inside_sheath), /NAN)
       sheathdensity(p)=mean(den(time_inside_sheath), /NAN)
       sheathdensitystd(p)=stddev(den(time_inside_sheath), /NAN)
       sheathtemperature(p)=mean(temp(time_inside_sheath), /NAN)
       sheathtemperaturestd(p)=stddev(temp(time_inside_sheath), /NAN)       
     endif else begin
      sheathspeed(p)=!VALUES.F_NAN
      sheathspeedstd(p)=!VALUES.F_NAN
      sheathdensity(p)=!VALUES.F_NAN
      sheathdensitystd(p)=!VALUES.F_NAN
      sheathtemperature(p)=!VALUES.F_NAN
      sheathtemperaturestd(p)=!VALUES.F_NAN       
     endelse
     
   
     ;all parameters inside MO
     time_inside_mo=where(timeb gt start_uly(p) and timeb lt end_uly(p))       
     if time_inside_mo(0) gt -1 then begin          
      bmax(p)=max(btot(time_inside_mo), /NAN)
      bmean(p)=mean(btot(time_inside_mo), /NAN)
      bstd(p)=stddev(btot(time_inside_mo), /NAN)
      bzmean(p)=mean(bz(time_inside_mo), /NAN)
      bzmin(p)=min(bz(time_inside_mo), /NAN)
      
      mospeed(p)=mean(vtot(time_inside_mo), /NAN)
      mospeedstd(p)=stddev(vtot(time_inside_mo), /NAN)
      
      modensity(p)=mean(den(time_inside_mo), /NAN)
      modensitystd(p)=stddev(den(time_inside_mo), /NAN)
      motemperature(p)=mean(temp(time_inside_mo), /NAN)
      motemperaturestd(p)=stddev(temp(time_inside_mo), /NAN)   


      duration(p)=(end_uly(p)-start_uly(p))/3600.
      heliodistance(p)=mean(helior(time_inside_mo), /NAN) /AU
      longitude(p)=mean(longit(time_inside_mo), /NAN)*!RADEG   
      latitude(p)=mean(latit(time_inside_mo), /NAN)*!RADEG 

      print, anytim(start_uly(p), /vms)
      print, 'uly/ICME max/mean/std btot:', bmax(p), bmean(p), bstd(p)
      print, 'duration', duration(p)
      
      ;get MVA results for magnetic obstacle
      mvub_results=mvub(bx(time_inside_mo),by(time_inside_mo),bz(time_inside_mo),timeb(time_inside_mo))     
  
      axis_longitude(p)=mvub_results(0)
      axis_latitude(p)=mvub_results(1)
      eigenvalue_ratio_2_3(p)=mvub_results(2)
      
      ;plot MO interval +/ 1 day      
      time_for_plot=where(timeb gt start_uly(p)-24*60*60. and timeb lt end_uly(p)+24*60*60.) 

      ;make in situ plots
      if make_plots_on gt 0 then r=make_plot(timeb(time_for_plot),btot(time_for_plot), $
             bx(time_for_plot), by(time_for_plot), bz(time_for_plot), $
             bmax(p),start_uly(p), end_uly(p),'Wind',icmecat(ulyicmes(p)).id,$
             lid,axis_longitude(p), axis_latitude(p),eigenvalue_ratio_2_3(p), ctype, icme_start_uly(p),icme_end_uly(p))

      
     endif           
    endfor     

  
 
; write results to master parameter variables
mbmax(ulyicmes)=bmax
mbmean(ulyicmes)=bmean
mbstd(ulyicmes)=bstd
mbzmean(ulyicmes)=bzmean
mbzmin(ulyicmes)=bzmin

mmo_duration(ulyicmes)=duration
mheliodistance(ulyicmes)=heliodistance
mlongitude(ulyicmes)=longitude
mlatitude(ulyicmes)=latitude
;MVA results
maxis_longitude(ulyicmes)=axis_longitude
maxis_latitude(ulyicmes)=axis_latitude
meigenvalue_ratio_2_3(ulyicmes)=eigenvalue_ratio_2_3
    
;speed results
mmospeed(ulyicmes)=mospeed
mmospeedstd(ulyicmes)=mospeedstd
msheathspeed(ulyicmes)=sheathspeed
msheathspeedstd(ulyicmes)=sheathspeedstd

;plasma results

mmodensity(ulyicmes)=modensity
mmodensitystd(ulyicmes)=modensitystd
mmotemperature(ulyicmes)=motemperature
mmotemperaturestd(ulyicmes)=motemperaturestd

msheathdensity(ulyicmes)=sheathdensity
msheathdensitystd(ulyicmes)=sheathdensitystd
msheathtemperature(ulyicmes)=sheathtemperature
msheathtemperaturestd(ulyicmes)=sheathtemperaturestd


    
msheathpdyn(ulyicmes)=!VALUES.F_NAN
msheathpdynstd(ulyicmes)=!VALUES.F_NAN
mmopdyn(ulyicmes)=!VALUES.F_NAN
mmopdynstd(ulyicmes)=!VALUES.F_NAN




;------------------------------------------------------ MESSENGER




;get ICME parameters (or rather MO) from VEX, for each year
;indices of all messenger icmes
mesicmes=where(icmecat.sc_insitu eq 'MESSENGER')


;get all start times - this is always mo_start_time
start_mes=anytim(icmecat(mesicmes).mo_start_time)

;first set all end_mes times to mo_end_time
end_mes=anytim(icmecat(mesicmes).mo_end_time)

;----------------for ICME start and end times
icme_start_mes=anytim(icmecat(mesicmes).icme_start_time)
icme_end_mes=anytim(icmecat(mesicmes).icme_end_time)



;define parameters 
bmax=dblarr(n_elements(mesicmes))
bmean=dblarr(n_elements(mesicmes))
bstd=dblarr(n_elements(mesicmes))
bzmean=dblarr(n_elements(mesicmes))
bzmin=dblarr(n_elements(mesicmes))

duration=dblarr(n_elements(mesicmes))
heliodistance=dblarr(n_elements(mesicmes))
longitude=dblarr(n_elements(mesicmes))
latitude=dblarr(n_elements(mesicmes))

;MVA
axis_longitude=dblarr(n_elements(mesicmes))
axis_latitude=dblarr(n_elements(mesicmes))
eigenvalue_ratio_2_3=dblarr(n_elements(mesicmes))











!p.background = 255
!p.color = 0
;!P.MULTI=[0,1,1]
;fuer single panel
;window, 2, xsize=1000,ysize=600, retain=2, xpos=100, ypos=500    
window, 2, xsize=1500,ysize=750, retain=2, xpos=100, ypos=500      


;for only solar wind

if ctype eq 'HEEQ' then savefilename=['/home/cmoestl/helcats/products/DATACAT/finaldata/MES_2007to2015_HEEQ_removed.sav']
if ctype eq 'SCEQ' then savefilename=['/home/cmoestl/helcats/products/DATACAT/finaldata/MES_2007to2015_SCEQ_removed.sav']

;for solar wind + bowshock data
;savefilename=['/home/cmoestl/helcats/products/DATACAT/finaldata/MES_2007to2014_RTN_1min_only_bowshock_single_array.sav']

restore, savefilename, /verbose
  
  ;get data for current year
  helior=mes.MES_RADIUS_IN_KM_HEEQ
  longit=mes.MES_LONGITUDE_IN_RADIANS_HEEQ
  latit=mes.MES_LATITUDE_IN_RADIANS_HEEQ
  bx=mes.bx
  by=mes.by
  bz=mes.bz
  btot=mes.btot
  timeb=mes.time
    
  ;go through all MESSENGER ICMEs and extract data inside mo_start to mo_end_time, and get parameters
   for p=0,n_elements(mesicmes)-1 do begin
     time_inside_mo=where(timeb gt start_mes(p) and timeb lt end_mes(p)) 
      
     if time_inside_mo(0) gt -1 then begin          
      bmax(p)=max(btot(time_inside_mo), /NAN)
      bmean(p)=mean(btot(time_inside_mo), /NAN)
      bstd(p)=stddev(btot(time_inside_mo), /NAN)
      bzmean(p)=mean(bz(time_inside_mo), /NAN)
      bzmin(p)=min(bz(time_inside_mo), /NAN)

      heliodistance(p)=mean(helior(time_inside_mo), /NAN) /AU
      longitude(p)=mean(longit(time_inside_mo), /NAN)*!RADEG 
      latitude(p)=mean(latit(time_inside_mo), /NAN)*!RADEG 

      duration(p)=(end_mes(p)-start_mes(p))/3600. ;hours
      print, anytim(start_mes(p), /vms)
      print, 'MES/ICME max/mean/std btot:', bmax(p), bmean(p), bstd(p)
      print, 'duration', duration(p)
      
      ;get MVA results for magnetic obstacle
      mvub_results=mvub(bx(time_inside_mo),by(time_inside_mo),bz(time_inside_mo),timeb(time_inside_mo))       
      axis_longitude(p)=mvub_results(0)
      axis_latitude(p)=mvub_results(1)
      eigenvalue_ratio_2_3(p)=mvub_results(2)  
       
      
      ;plot MO interval +/ 1 day
      time_for_plot=where(timeb gt start_mes(p)-24*60*60. and timeb lt end_mes(p)+24*60*60.) 

      ;make in situ plots
      
      if make_plots_on gt 0 then $
            r=make_plot(timeb(time_for_plot),btot(time_for_plot), $
             bx(time_for_plot), by(time_for_plot), bz(time_for_plot), $
             bmax(p),start_mes(p), end_mes(p),'MESSENGER',icmecat(mesicmes(p)).id,lid, $
             axis_longitude(p), axis_latitude(p),eigenvalue_ratio_2_3(p), ctype,icme_start_mes(p),icme_end_mes(p))
    
     endif      
     
    endfor    
 

 
; write results to master parameter variables
mbmax(mesicmes)=bmax
mbmean(mesicmes)=bmean
mbstd(mesicmes)=bstd
mbzmean(mesicmes)=bzmean
mbzmin(mesicmes)=bzmin


mmo_duration(mesicmes)=duration
mheliodistance(mesicmes)=heliodistance
mlongitude(mesicmes)=longitude
mlatitude(mesicmes)=latitude
;MVA results
maxis_longitude(mesicmes)=axis_longitude
maxis_latitude(mesicmes)=axis_latitude
meigenvalue_ratio_2_3(mesicmes)=eigenvalue_ratio_2_3

;speed
mmospeed(mesicmes)=!VALUES.F_NAN
mmospeedstd(mesicmes)=!VALUES.F_NAN
msheathspeed(mesicmes)=!VALUES.F_NAN
msheathspeedstd(mesicmes)=!VALUES.F_NAN


mmodensity(mesicmes)=!VALUES.F_NAN
mmodensitystd(mesicmes)=!VALUES.F_NAN
mmotemperature(mesicmes)=!VALUES.F_NAN
mmotemperaturestd(mesicmes)=!VALUES.F_NAN

msheathdensity(mesicmes)=!VALUES.F_NAN
msheathdensitystd(mesicmes)=!VALUES.F_NAN
msheathtemperature(mesicmes)=!VALUES.F_NAN
msheathtemperaturestd(mesicmes)=!VALUES.F_NAN


       
    
msheathpdyn(mesicmes)=!VALUES.F_NAN
msheathpdynstd(mesicmes)=!VALUES.F_NAN
mmopdyn(mesicmes)=!VALUES.F_NAN
mmopdynstd(mesicmes)=!VALUES.F_NAN


    
    
    
    
    
    
    
    
    

;------------------------------------------------------ VEX


;get ICME parameters (or rather MO) from VEX, for each year
;indices of all messenger icmes
vexicmes=where(icmecat.sc_insitu eq 'VEX')
;get all start times -this is mo_start_time
start_vex=anytim(icmecat(vexicmes).mo_start_time)
end_vex=anytim(icmecat(vexicmes).mo_end_time)

;----------------for ICME start and end times

;there is no ICME end time for VEX, and only few ICME start times

icme_start_vex=anytim(icmecat(vexicmes).icme_start_time)
icme_end_vex=anytim(icmecat(vexicmes).icme_end_time)










;define parameters
bmax=dblarr(n_elements(vexicmes))
bmean=dblarr(n_elements(vexicmes))
bstd=dblarr(n_elements(vexicmes))
bzmean=dblarr(n_elements(vexicmes))
bzmin=dblarr(n_elements(vexicmes))


duration=dblarr(n_elements(vexicmes))
heliodistance=dblarr(n_elements(vexicmes))
longitude=dblarr(n_elements(vexicmes))
latitude=dblarr(n_elements(vexicmes))
;MVA
axis_longitude=dblarr(n_elements(vexicmes))
axis_latitude=dblarr(n_elements(vexicmes))
eigenvalue_ratio_2_3=dblarr(n_elements(vexicmes))




if ctype eq 'HEEQ' then savefilename=['/home/cmoestl/helcats/products/DATACAT/finaldata/VEX_2007to2014_HEEQ_removed.sav']
if ctype eq 'SCEQ' then savefilename=['/home/cmoestl/helcats/products/DATACAT/finaldata/VEX_2007to2014_SCEQ_removed.sav']

restore, savefilename, /verbose
	
  ;get data for current year
  helior=vex.VEX_RADIUS_IN_KM_HEEQ
  longit=vex.VEX_LONGITUDE_IN_RADIANS_HEEQ
  bx=vex.bx
  by=vex.by
  bz=vex.bz
  btot=vex.btot
  timeb=vex.time
 
  ;go through all ICMEs and extract data inside mo_start to mo_end_time, and get parameters
   for p=0,n_elements(vexicmes)-1 do begin
     time_inside_mo=where(timeb gt start_vex(p) and timeb lt end_vex(p)) 
      
     if time_inside_mo(0) gt -1 then begin          
      bmax(p)=max(btot(time_inside_mo), /NAN)
      bmean(p)=mean(btot(time_inside_mo), /NAN)
      bstd(p)=stddev(btot(time_inside_mo), /NAN)
      bzmean(p)=mean(bz(time_inside_mo), /NAN)
      bzmin(p)=min(bz(time_inside_mo), /NAN)

      heliodistance(p)=mean(helior(time_inside_mo), /NAN) /AU
      longitude(p)=mean(longit(time_inside_mo), /NAN)*!RADEG   
      latitude(p)=mean(latit(time_inside_mo), /NAN)*!RADEG 

      duration(p)=(end_vex(p)-start_vex(p))/3600.
      print, anytim(start_vex(p), /vms)
      print, 'VEX/ICME max/mean/std btot:', bmax(p), bmean(p), bstd(p)
      print, 'duration', duration(p)
      
      ;get MVA results for magnetic obstacle
      mvub_results=mvub(bx(time_inside_mo),by(time_inside_mo),bz(time_inside_mo),timeb(time_inside_mo))     
      axis_longitude(p)=mvub_results(0)
      axis_latitude(p)=mvub_results(1)
      eigenvalue_ratio_2_3(p)=mvub_results(2)
      
      ;plot MO interval +/ 1 day      
      time_for_plot=where(timeb gt start_vex(p)-24*60*60. and timeb lt end_vex(p)+24*60*60.) 


      ;make in situ plots
      if make_plots_on gt 0 then $
        r=make_plot(timeb(time_for_plot),btot(time_for_plot), $
             bx(time_for_plot), by(time_for_plot), bz(time_for_plot), $
             bmax(p),start_vex(p), end_vex(p),'VEX',icmecat(vexicmes(p)).id,$
             lid,axis_longitude(p), axis_latitude(p),eigenvalue_ratio_2_3(p), ctype,$
             icme_start_vex(p),icme_end_vex(p))

     
      
      
     endif      
     
    endfor    
 
  
 
; write results to master parameter variables
mbmax(vexicmes)=bmax
mbmean(vexicmes)=bmean
mbstd(vexicmes)=bstd
mbzmean(vexicmes)=bzmean
mbzmin(vexicmes)=bzmin


mmo_duration(vexicmes)=duration
mheliodistance(vexicmes)=heliodistance
mlongitude(vexicmes)=longitude
mlatitude(vexicmes)=latitude

;MVA results
maxis_longitude(vexicmes)=axis_longitude
maxis_latitude(vexicmes)=axis_latitude
meigenvalue_ratio_2_3(vexicmes)=eigenvalue_ratio_2_3


;set all VEX speeds to NaN
mmospeed(vexicmes)=!VALUES.F_NAN
msheathspeed(vexicmes)=!VALUES.F_NAN
mmospeedstd(vexicmes)=!VALUES.F_NAN
msheathspeedstd(vexicmes)=!VALUES.F_NAN
  
mmodensity(vexicmes)=!VALUES.F_NAN
mmodensitystd(vexicmes)=!VALUES.F_NAN
mmotemperature(vexicmes)=!VALUES.F_NAN
mmotemperaturestd(vexicmes)=!VALUES.F_NAN

msheathdensity(vexicmes)=!VALUES.F_NAN
msheathdensitystd(vexicmes)=!VALUES.F_NAN
msheathtemperature(vexicmes)=!VALUES.F_NAN
msheathtemperaturestd(vexicmes)=!VALUES.F_NAN

  
  
  
    
msheathpdyn(vexicmes)=!VALUES.F_NAN
msheathpdynstd(vexicmes)=!VALUES.F_NAN
mmopdyn(vexicmes)=!VALUES.F_NAN
mmopdynstd(vexicmes)=!VALUES.F_NAN

  
  
  
  
  
  
  
  
  
  
  
  
  
;------------------------------------------------------ STEREO-A



;get ICME parameters (or rather MO) from VEX, for each year
staicmes=where(icmecat.sc_insitu eq 'STEREO-A')

start_sta=anytim(icmecat(staicmes).mo_start_time)
end_sta=anytim(icmecat(staicmes).mo_end_time)

;----------------for ICME start and end times
icme_start_sta=anytim(icmecat(staicmes).icme_start_time)
icme_end_sta=anytim(icmecat(staicmes).icme_end_time)


;define parameters
bmax=dblarr(n_elements(staicmes))
bmean=dblarr(n_elements(staicmes))
bstd=dblarr(n_elements(staicmes))
bzmean=dblarr(n_elements(staicmes))
bzmin=dblarr(n_elements(staicmes))


duration=dblarr(n_elements(staicmes))
heliodistance=dblarr(n_elements(staicmes))
longitude=dblarr(n_elements(staicmes))
latitude=dblarr(n_elements(staicmes))
;MVA
axis_longitude=dblarr(n_elements(staicmes))
axis_latitude=dblarr(n_elements(staicmes))
eigenvalue_ratio_2_3=dblarr(n_elements(staicmes))

mospeed=dblarr(n_elements(staicmes))
mospeedstd=dblarr(n_elements(staicmes))
sheathspeed=dblarr(n_elements(staicmes))
sheathspeedstd=dblarr(n_elements(staicmes))


modensity=dblarr(n_elements(staicmes))
modensitystd=dblarr(n_elements(staicmes))
motemperature=dblarr(n_elements(staicmes))
motemperaturestd=dblarr(n_elements(staicmes))

sheathdensity=dblarr(n_elements(staicmes))
sheathdensitystd=dblarr(n_elements(staicmes))
sheathtemperature=dblarr(n_elements(staicmes))
sheathtemperaturestd=dblarr(n_elements(staicmes))

mopdyn=dblarr(n_elements(staicmes))
mopdynstd=dblarr(n_elements(staicmes))
sheathpdyn=dblarr(n_elements(staicmes))
sheathpdynstd=dblarr(n_elements(staicmes))
       


if ctype eq 'HEEQ' then savefilename=['/home/cmoestl/helcats/products/DATACAT/finaldata/STA_2007to2015_HEEQ.sav']
if ctype eq 'SCEQ' then savefilename=['/home/cmoestl/helcats/products/DATACAT/finaldata/STA_2007to2015_SCEQ.sav']
restore, savefilename, /verbose
	
  ;get data for current year
  helior=sta.STA_RADIUS_IN_KM_HEEQ
  longit=sta.STA_LONGITUDE_IN_RADIANS_HEEQ
  bx=sta.bx
  by=sta.by
  bz=sta.bz
  btot=sta.btot
  timeb=sta.time
  vtot=sta.vtot
  den=sta.density
  temp=sta.temperature
  
    
   for p=0,n_elements(staicmes)-1 do begin
  
  
    ;speeds in sheath
      time_inside_sheath=where(timeb gt icme_start_sta(p) and timeb lt start_sta(p))       
      if time_inside_sheath(0) gt -1 then begin     
       sheathspeed(p)=mean(vtot(time_inside_sheath), /NAN)
       sheathspeedstd(p)=stddev(vtot(time_inside_sheath), /NAN)
       sheathdensity(p)=mean(den(time_inside_sheath), /NAN)
       sheathdensitystd(p)=stddev(den(time_inside_sheath), /NAN)
       sheathtemperature(p)=mean(temp(time_inside_sheath), /NAN)
       sheathtemperaturestd(p)=stddev(temp(time_inside_sheath), /NAN)      
       sheathpdyn(p)=mean(dynamic_pressure(den(time_inside_sheath),vtot(time_inside_sheath)), /NAN)
       sheathpdynstd(p)=stddev(dynamic_pressure(den(time_inside_sheath),vtot(time_inside_sheath)), /NAN)
      endif else begin
       sheathspeed(p)=!VALUES.F_NAN
       sheathspeedstd(p)=!VALUES.F_NAN
       sheathdensity(p)=!VALUES.F_NAN
       sheathdensitystd(p)=!VALUES.F_NAN
       sheathtemperature(p)=!VALUES.F_NAN
       sheathtemperaturestd(p)=!VALUES.F_NAN
       sheathpdyn(p)=!VALUES.F_NAN
       sheathpdynstd(p)=!VALUES.F_NAN  
     
      endelse
   
     time_inside_mo=where(timeb gt start_sta(p) and timeb lt end_sta(p))       
   
     if time_inside_mo(0) gt -1 then begin          
      bmax(p)=max(btot(time_inside_mo), /NAN)
      bmean(p)=mean(btot(time_inside_mo), /NAN)
      bstd(p)=stddev(btot(time_inside_mo), /NAN)
      bzmean(p)=mean(bz(time_inside_mo), /NAN)
      bzmin(p)=min(bz(time_inside_mo), /NAN)
      
      mospeed(p)=mean(vtot(time_inside_mo), /NAN)
      mospeedstd(p)=stddev(vtot(time_inside_mo), /NAN)
            
      modensity(p)=mean(den(time_inside_mo), /NAN)
      modensitystd(p)=stddev(den(time_inside_mo), /NAN)
      motemperature(p)=mean(temp(time_inside_mo), /NAN)
      motemperaturestd(p)=stddev(temp(time_inside_mo), /NAN)   
      
      
      mopdyn(p)=mean(dynamic_pressure(den(time_inside_mo),vtot(time_inside_mo)), /NAN)
      mopdynstd(p)=stddev(dynamic_pressure(den(time_inside_mo),vtot(time_inside_mo)), /NAN)


      duration(p)=(end_sta(p)-start_sta(p))/3600.
      heliodistance(p)=mean(helior(time_inside_mo), /NAN) /AU
      longitude(p)=mean(longit(time_inside_mo), /NAN)*!RADEG   
      latitude(p)=mean(latit(time_inside_mo), /NAN)*!RADEG 

      print, anytim(start_sta(p), /vms)
      print, 'STA/ICME max/mean/std btot:', bmax(p), bmean(p), bstd(p)
      print, 'duration', duration(p)
      
      ;get MVA results for magnetic obstacle
      mvub_results=mvub(bx(time_inside_mo),by(time_inside_mo),bz(time_inside_mo),timeb(time_inside_mo))     
 
      axis_longitude(p)=mvub_results(0)
      axis_latitude(p)=mvub_results(1)
      eigenvalue_ratio_2_3(p)=mvub_results(2)
      
      ;plot MO interval +/ 1 day      
      time_for_plot=where(timeb gt start_sta(p)-24*60*60. and timeb lt end_sta(p)+24*60*60.) 

      ;make in situ plots
      if make_plots_on gt 0 then r=make_plot(timeb(time_for_plot),btot(time_for_plot), $
                bx(time_for_plot), by(time_for_plot), bz(time_for_plot), $
                bmax(p),start_sta(p), end_sta(p),'STEREO-A',icmecat(staicmes(p)).id,$
                lid,axis_longitude(p), axis_latitude(p),eigenvalue_ratio_2_3(p), ctype, icme_start_sta(p),icme_end_sta(p))
 
      
     endif           
    endfor    
    
  
 
; write results to master parameter variables
mbmax(staicmes)=bmax
mbmean(staicmes)=bmean
mbstd(staicmes)=bstd
mbzmean(staicmes)=bzmean
mbzmin(staicmes)=bzmin


mmo_duration(staicmes)=duration
mheliodistance(staicmes)=heliodistance
mlongitude(staicmes)=longitude
mlatitude(staicmes)=latitude

;MVA results
maxis_longitude(staicmes)=axis_longitude
maxis_latitude(staicmes)=axis_latitude
meigenvalue_ratio_2_3(staicmes)=eigenvalue_ratio_2_3


;speed results
mmospeed(staicmes)=mospeed
mmospeedstd(staicmes)=mospeedstd
msheathspeed(staicmes)=sheathspeed
msheathspeedstd(staicmes)=sheathspeedstd


;plasma results


mmodensity(staicmes)=modensity
mmodensitystd(staicmes)=modensitystd
mmotemperature(staicmes)=motemperature
mmotemperaturestd(staicmes)=motemperaturestd

msheathdensity(staicmes)=sheathdensity
msheathdensitystd(staicmes)=sheathdensitystd
msheathtemperature(staicmes)=sheathtemperature
msheathtemperaturestd(staicmes)=sheathtemperaturestd



msheathpdyn(staicmes)=sheathpdyn
msheathpdynstd(staicmes)=sheathpdynstd
mmopdyn(staicmes)=mopdyn
mmopdynstd(staicmes)=mopdynstd










  
  
;------------------------------------------------------ STEREO-B



;get ICME parameters (or rather MO) from VEX, for each year
stbicmes=where(icmecat.sc_insitu eq 'STEREO-B')
start_stb=anytim(icmecat(stbicmes).mo_start_time)
end_stb=anytim(icmecat(stbicmes).mo_end_time)

;----------------for ICME start and end times
icme_start_stb=anytim(icmecat(stbicmes).icme_start_time)
icme_end_stb=anytim(icmecat(stbicmes).icme_end_time)



;define parameters
bmax=dblarr(n_elements(stbicmes))

bmean=dblarr(n_elements(stbicmes))
bstd=dblarr(n_elements(stbicmes))
bzmean=dblarr(n_elements(stbicmes))
bzmin=dblarr(n_elements(stbicmes))

duration=dblarr(n_elements(stbicmes))
heliodistance=dblarr(n_elements(stbicmes))
longitude=dblarr(n_elements(stbicmes))
latitude=dblarr(n_elements(stbicmes))


;MVA
axis_longitude=dblarr(n_elements(stbicmes))
axis_latitude=dblarr(n_elements(stbicmes))
eigenvalue_ratio_2_3=dblarr(n_elements(stbicmes))


mospeed=dblarr(n_elements(stbicmes))
mospeedstd=dblarr(n_elements(stbicmes))
sheathspeed=dblarr(n_elements(stbicmes))
sheathspeedstd=dblarr(n_elements(stbicmes))

sheathdensity=dblarr(n_elements(stbicmes))
sheathdensitystd=dblarr(n_elements(stbicmes))
sheathtemperature=dblarr(n_elements(stbicmes))
sheathtemperaturestd=dblarr(n_elements(stbicmes))

modensity=dblarr(n_elements(stbicmes))
modensitystd=dblarr(n_elements(stbicmes))
motemperature=dblarr(n_elements(stbicmes))
motemperaturestd=dblarr(n_elements(stbicmes))
       

mopdyn=dblarr(n_elements(stbicmes))
mopdynstd=dblarr(n_elements(stbicmes))
sheathpdyn=dblarr(n_elements(stbicmes))
sheathpdynstd=dblarr(n_elements(stbicmes))
       


if ctype eq 'HEEQ' then savefilename=['/home/cmoestl/helcats/products/DATACAT/finaldata/STB_2007to2014_HEEQ.sav']
if ctype eq 'SCEQ' then savefilename=['/home/cmoestl/helcats/products/DATACAT/finaldata/STB_2007to2014_SCEQ.sav']
restore, savefilename, /verbose
	
  ;get data for current year
  helior=stb.STB_RADIUS_IN_KM_HEEQ
  longit=stb.STB_LONGITUDE_IN_RADIANS_HEEQ
  latit=stb.STB_LATITUDE_IN_RADIANS_HEEQ
  bx=stb.bx
  by=stb.by
  bz=stb.bz
  btot=stb.btot
  timeb=stb.time
  vtot=stb.vtot
  den=stb.density
  temp=stb.temperature
    
    
   for p=0,n_elements(stbicmes)-1 do begin
   
    ;speeds in sheath
     time_inside_sheath=where(timeb gt icme_start_stb(p) and timeb lt start_stb(p))       
     if time_inside_sheath(0) gt -1 then begin     
       sheathspeed(p)=mean(vtot(time_inside_sheath), /NAN)
       sheathspeedstd(p)=stddev(vtot(time_inside_sheath), /NAN)
       sheathdensity(p)=mean(den(time_inside_sheath), /NAN)
       sheathdensitystd(p)=stddev(den(time_inside_sheath), /NAN)
       sheathtemperature(p)=mean(temp(time_inside_sheath), /NAN)
       sheathtemperaturestd(p)=stddev(temp(time_inside_sheath), /NAN)       
       sheathpdyn(p)=mean(dynamic_pressure(den(time_inside_sheath),vtot(time_inside_sheath)), /NAN)
       sheathpdynstd(p)=stddev(dynamic_pressure(den(time_inside_sheath),vtot(time_inside_sheath)), /NAN)
     endif else begin ;if there is no sheath
      sheathspeed(p)=!VALUES.F_NAN
      sheathspeedstd(p)=!VALUES.F_NAN
      sheathdensity(p)=!VALUES.F_NAN
      sheathdensitystd(p)=!VALUES.F_NAN
      sheathtemperature(p)=!VALUES.F_NAN
      sheathtemperaturestd(p)=!VALUES.F_NAN     
      sheathpdyn(p)=!VALUES.F_NAN
      sheathpdynstd(p)=!VALUES.F_NAN         
     endelse
   
   
     time_inside_mo=where(timeb gt start_stb(p) and timeb lt end_stb(p))       
     if time_inside_mo(0) gt -1 then begin          
      bmax(p)=max(btot(time_inside_mo), /NAN)
      bmean(p)=mean(btot(time_inside_mo), /NAN)
      bstd(p)=stddev(btot(time_inside_mo), /NAN)
      bzmean(p)=mean(bz(time_inside_mo), /NAN)
      bzmin(p)=min(bz(time_inside_mo), /NAN)
      
      mospeed(p)=mean(vtot(time_inside_mo), /NAN)
      mospeedstd(p)=stddev(vtot(time_inside_mo), /NAN)

      modensity(p)=mean(den(time_inside_mo), /NAN)
      modensitystd(p)=stddev(den(time_inside_mo), /NAN)
      motemperature(p)=mean(temp(time_inside_mo), /NAN)
      motemperaturestd(p)=stddev(temp(time_inside_mo), /NAN)   
      
      mopdyn(p)=mean(dynamic_pressure(den(time_inside_mo),vtot(time_inside_mo)), /NAN)
      mopdynstd(p)=stddev(dynamic_pressure(den(time_inside_mo),vtot(time_inside_mo)), /NAN)



      duration(p)=(end_stb(p)-start_stb(p))/3600.
      heliodistance(p)=mean(helior(time_inside_mo), /NAN) /AU
      longitude(p)=mean(longit(time_inside_mo), /NAN)*!RADEG   
      latitude(p)=mean(latit(time_inside_mo), /NAN)*!RADEG 

      print, anytim(start_stb(p), /vms)
      print, 'stb/ICME max/mean/std btot:', bmax(p), bmean(p), bstd(p)
      print, 'duration', duration(p)
      
      ;get MVA results for magnetic obstacle
      mvub_results=mvub(bx(time_inside_mo),by(time_inside_mo),bz(time_inside_mo),timeb(time_inside_mo))       
      axis_longitude(p)=mvub_results(0)
      axis_latitude(p)=mvub_results(1)
      eigenvalue_ratio_2_3(p)=mvub_results(2)
       
      ;plot MO interval +/ 1 day      
      time_for_plot=where(timeb gt start_stb(p)-24*60*60. and timeb lt end_stb(p)+24*60*60.) 

      ;make in situ plots
      if make_plots_on gt 0 then r=make_plot(timeb(time_for_plot),btot(time_for_plot), $
                    bx(time_for_plot), by(time_for_plot), bz(time_for_plot), $
                     bmax(p),start_stb(p), end_stb(p),'STEREO-B',icmecat(stbicmes(p)).id,$
                     lid,axis_longitude(p), axis_latitude(p),eigenvalue_ratio_2_3(p), ctype, icme_start_stb(p),icme_end_stb(p))
      
      
     endif           
    endfor     
  
 
; write results to master parameter variables
mbmax(stbicmes)=bmax
mbmean(stbicmes)=bmean
mbstd(stbicmes)=bstd
mbzmean(stbicmes)=bzmean
mbzmin(stbicmes)=bzmin


mmo_duration(stbicmes)=duration
mheliodistance(stbicmes)=heliodistance
mlongitude(stbicmes)=longitude
mlatitude(stbicmes)=latitude

;MVA results
maxis_longitude(stbicmes)=axis_longitude
maxis_latitude(stbicmes)=axis_latitude
meigenvalue_ratio_2_3(stbicmes)=eigenvalue_ratio_2_3


;speed results
mmospeed(stbicmes)=mospeed
mmospeedstd(stbicmes)=mospeedstd
msheathspeed(stbicmes)=sheathspeed
msheathspeedstd(stbicmes)=sheathspeedstd


mmodensity(stbicmes)=modensity
mmodensitystd(stbicmes)=modensitystd
mmotemperature(stbicmes)=motemperature
mmotemperaturestd(stbicmes)=motemperaturestd

msheathdensity(stbicmes)=sheathdensity
msheathdensitystd(stbicmes)=sheathdensitystd
msheathtemperature(stbicmes)=sheathtemperature
msheathtemperaturestd(stbicmes)=sheathtemperaturestd



msheathpdyn(stbicmes)=sheathpdyn
msheathpdynstd(stbicmes)=sheathpdynstd
mmopdyn(stbicmes)=mopdyn
mmopdynstd(stbicmes)=mopdynstd





  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
;------------------------------------------------------ Wind


;get ICME parameters (or rather MO)  for each year
winicmes=where(icmecat.sc_insitu eq 'Wind')
start_win=anytim(icmecat(winicmes).mo_start_time)
end_win=anytim(icmecat(winicmes).mo_end_time)

;----------------for ICME start and end times
icme_start_win=anytim(icmecat(winicmes).icme_start_time)
icme_end_win=anytim(icmecat(winicmes).icme_end_time)


;define parameters
bmax=dblarr(n_elements(winicmes))
bmean=dblarr(n_elements(winicmes))
bstd=dblarr(n_elements(winicmes))
bzmean=dblarr(n_elements(winicmes))
bzmin=dblarr(n_elements(winicmes))



duration=dblarr(n_elements(winicmes))
heliodistance=dblarr(n_elements(winicmes))
longitude=dblarr(n_elements(winicmes))
latitude=dblarr(n_elements(winicmes))

axis_longitude=dblarr(n_elements(winicmes))
axis_latitude=dblarr(n_elements(winicmes))
eigenvalue_ratio_2_3=dblarr(n_elements(winicmes))

mospeed=dblarr(n_elements(winicmes))
mospeedstd=dblarr(n_elements(winicmes))
sheathspeed=dblarr(n_elements(winicmes))
sheathspeedstd=dblarr(n_elements(winicmes))

modensity=dblarr(n_elements(winicmes))
modensitystd=dblarr(n_elements(winicmes))
motemperature=dblarr(n_elements(winicmes))
motemperaturestd=dblarr(n_elements(winicmes))

sheathdensity=dblarr(n_elements(winicmes))
sheathdensitystd=dblarr(n_elements(winicmes))
sheathtemperature=dblarr(n_elements(winicmes))
sheathtemperaturestd=dblarr(n_elements(winicmes))

mopdyn=dblarr(n_elements(winicmes))
mopdynstd=dblarr(n_elements(winicmes))
sheathpdyn=dblarr(n_elements(winicmes))
sheathpdynstd=dblarr(n_elements(winicmes))
       


;same file for both HEEQ and SCEQ because no difference
;savefilename=['/home/cmoestl/helcats/products/DATACAT/finaldata/WIND_2007to2018_HEEQ.sav']
savefile=['/home/cmoestl/helcats/products/DATACAT/finaldata/WIND_2007to2018_HEEQ_plasma_median21.sav']

restore, savefilename, /verbose
	
  ;get data 
  helior=wind.WIND_RADIUS_IN_KM_HEEQ
  longit=wind.WIND_LONGITUDE_IN_RADIANS_HEEQ
  latit=wind.WIND_LATITUDE_IN_RADIANS_HEEQ
  bx=wind.bx
  by=wind.by
  bz=wind.bz
  btot=wind.btot
  timeb=wind.time
  vtot=wind.vtot
  den=wind.density
  temp=wind.temperature
    
    
   for p=0,n_elements(winicmes)-1 do begin
   
     ;speeds in sheath
     time_inside_sheath=where(timeb gt icme_start_win(p) and timeb lt start_win(p))       
     if time_inside_sheath(0) gt -1 then begin     
       sheathspeed(p)=mean(vtot(time_inside_sheath), /NAN)
       sheathspeedstd(p)=stddev(vtot(time_inside_sheath), /NAN)
       sheathdensity(p)=mean(den(time_inside_sheath), /NAN)
       sheathdensitystd(p)=stddev(den(time_inside_sheath), /NAN)
       sheathtemperature(p)=mean(temp(time_inside_sheath), /NAN)
       sheathtemperaturestd(p)=stddev(temp(time_inside_sheath), /NAN)       
       sheathpdyn(p)=mean(dynamic_pressure(den(time_inside_sheath),vtot(time_inside_sheath)), /NAN)
       sheathpdynstd(p)=stddev(dynamic_pressure(den(time_inside_sheath),vtot(time_inside_sheath)), /NAN)
       
     endif else begin
      sheathspeed(p)=!VALUES.F_NAN
      sheathspeedstd(p)=!VALUES.F_NAN
      sheathdensity(p)=!VALUES.F_NAN
      sheathdensitystd(p)=!VALUES.F_NAN
      sheathtemperature(p)=!VALUES.F_NAN
      sheathtemperaturestd(p)=!VALUES.F_NAN     
      sheathpdyn(p)=!VALUES.F_NAN
      sheathpdynstd(p)=!VALUES.F_NAN  
     endelse
     
   
     ;all parameters inside MO
     time_inside_mo=where(timeb gt start_win(p) and timeb lt end_win(p))       
     if time_inside_mo(0) gt -1 then begin          
      bmax(p)=max(btot(time_inside_mo), /NAN)
      bmean(p)=mean(btot(time_inside_mo), /NAN)
      bstd(p)=stddev(btot(time_inside_mo), /NAN)
      bzmean(p)=mean(bz(time_inside_mo), /NAN)
      bzmin(p)=min(bz(time_inside_mo), /NAN)
      
      mospeed(p)=mean(vtot(time_inside_mo), /NAN)
      mospeedstd(p)=stddev(vtot(time_inside_mo), /NAN)
      
      modensity(p)=mean(den(time_inside_mo), /NAN)
      modensitystd(p)=stddev(den(time_inside_mo), /NAN)
      motemperature(p)=mean(temp(time_inside_mo), /NAN)
      motemperaturestd(p)=stddev(temp(time_inside_mo), /NAN)   

      mopdyn(p)=mean(dynamic_pressure(den(time_inside_mo),vtot(time_inside_mo)), /NAN)
      mopdynstd(p)=stddev(dynamic_pressure(den(time_inside_mo),vtot(time_inside_mo)), /NAN)


      duration(p)=(end_win(p)-start_win(p))/3600.
      heliodistance(p)=mean(helior(time_inside_mo), /NAN) /AU
      longitude(p)=mean(longit(time_inside_mo), /NAN)*!RADEG   
      latitude(p)=mean(latit(time_inside_mo), /NAN)*!RADEG 

      print, anytim(start_win(p), /vms)
      print, 'win/ICME max/mean/std btot:', bmax(p), bmean(p), bstd(p)
      print, 'duration', duration(p)
      
      ;get MVA results for magnetic obstacle
      mvub_results=mvub(bx(time_inside_mo),by(time_inside_mo),bz(time_inside_mo),timeb(time_inside_mo))     
  
      axis_longitude(p)=mvub_results(0)
      axis_latitude(p)=mvub_results(1)
      eigenvalue_ratio_2_3(p)=mvub_results(2)
      
      ;plot MO interval +/ 1 day      
      time_for_plot=where(timeb gt start_win(p)-24*60*60. and timeb lt end_win(p)+24*60*60.) 

      ;make in situ plots
      if make_plots_on gt 0 then r=make_plot(timeb(time_for_plot),btot(time_for_plot), $
             bx(time_for_plot), by(time_for_plot), bz(time_for_plot), $
             bmax(p),start_win(p), end_win(p),'Wind',icmecat(winicmes(p)).id,$
             lid,axis_longitude(p), axis_latitude(p),eigenvalue_ratio_2_3(p), ctype, icme_start_win(p),icme_end_win(p))

      
     endif           
    endfor     

  
 
; write results to master parameter variables
mbmax(winicmes)=bmax
mbmean(winicmes)=bmean
mbstd(winicmes)=bstd
mbzmean(winicmes)=bzmean
mbzmin(winicmes)=bzmin

mmo_duration(winicmes)=duration
mheliodistance(winicmes)=heliodistance
mlongitude(winicmes)=longitude
mlatitude(winicmes)=latitude
;MVA results
maxis_longitude(winicmes)=axis_longitude
maxis_latitude(winicmes)=axis_latitude
meigenvalue_ratio_2_3(winicmes)=eigenvalue_ratio_2_3
    
;speed results
mmospeed(winicmes)=mospeed
mmospeedstd(winicmes)=mospeedstd
msheathspeed(winicmes)=sheathspeed
msheathspeedstd(winicmes)=sheathspeedstd

;plasma results

mmodensity(winicmes)=modensity
mmodensitystd(winicmes)=modensitystd
mmotemperature(winicmes)=motemperature
mmotemperaturestd(winicmes)=motemperaturestd

msheathdensity(winicmes)=sheathdensity
msheathdensitystd(winicmes)=sheathdensitystd
msheathtemperature(winicmes)=sheathtemperature
msheathtemperaturestd(winicmes)=sheathtemperaturestd





msheathpdyn(winicmes)=sheathpdyn
msheathpdynstd(winicmes)=sheathpdynstd
mmopdyn(winicmes)=mopdyn
mmopdynstd(winicmes)=mopdynstd
















  
;------------------------------------------------------ MAVEN


;get ICME parameters (or rather MO)  for each year
mavicmes=where(icmecat.sc_insitu eq 'MAVEN')
start_mav=anytim(icmecat(mavicmes).mo_start_time)
end_mav=anytim(icmecat(mavicmes).mo_end_time)

;----------------for ICME start and end times
icme_start_mav=anytim(icmecat(mavicmes).icme_start_time)
icme_end_mav=anytim(icmecat(mavicmes).icme_end_time)


;define parameters
bmax=dblarr(n_elements(mavicmes))
bmean=dblarr(n_elements(mavicmes))
bstd=dblarr(n_elements(mavicmes))
bzmean=dblarr(n_elements(mavicmes))
bzmin=dblarr(n_elements(mavicmes))



duration=dblarr(n_elements(mavicmes))
heliodistance=dblarr(n_elements(mavicmes))
longitude=dblarr(n_elements(mavicmes))
latitude=dblarr(n_elements(mavicmes))

axis_longitude=dblarr(n_elements(mavicmes))
axis_latitude=dblarr(n_elements(mavicmes))
eigenvalue_ratio_2_3=dblarr(n_elements(mavicmes))

mospeed=dblarr(n_elements(mavicmes))
mospeedstd=dblarr(n_elements(mavicmes))
sheathspeed=dblarr(n_elements(mavicmes))
sheathspeedstd=dblarr(n_elements(mavicmes))

modensity=dblarr(n_elements(mavicmes))
modensitystd=dblarr(n_elements(mavicmes))
motemperature=dblarr(n_elements(mavicmes))
motemperaturestd=dblarr(n_elements(mavicmes))

sheathdensity=dblarr(n_elements(mavicmes))
sheathdensitystd=dblarr(n_elements(mavicmes))
sheathtemperature=dblarr(n_elements(mavicmes))
sheathtemperaturestd=dblarr(n_elements(mavicmes))

mopdyn=dblarr(n_elements(mavicmes))
mopdynstd=dblarr(n_elements(mavicmes))
sheathpdyn=dblarr(n_elements(mavicmes))
sheathpdynstd=dblarr(n_elements(mavicmes))
       
       


;same file for both HEEQ and SCEQ because no difference
savefilename=['/home/cmoestl/helcats/products/DATACAT/finaldata/MAVEN_2014to2018_MSO_removed_orbit.sav']
restore, savefilename, /verbose
	
  ;get data 
  helior=mav.MAVEN_RADIUS_IN_KM_HEEQ
  longit=mav.MAVEN_LONGITUDE_IN_RADIANS_HEEQ
  latit=mav.MAVEN_LATITUDE_IN_RADIANS_HEEQ
  bx=mav.bx
  by=mav.by
  bz=mav.bz
  btot=mav.btot
  timeb=mav.time
  vtot=mav.vtot
  den=mav.n
  temp=mav.t
    
    
   for p=0,n_elements(mavicmes)-1 do begin
   
     ;speeds in sheath
     time_inside_sheath=where(timeb gt icme_start_mav(p) and timeb lt start_mav(p))       
     if time_inside_sheath(0) gt -1 then begin     
       sheathspeed(p)=mean(vtot(time_inside_sheath), /NAN)
       sheathspeedstd(p)=stddev(vtot(time_inside_sheath), /NAN)
       sheathdensity(p)=mean(den(time_inside_sheath), /NAN)
       sheathdensitystd(p)=stddev(den(time_inside_sheath), /NAN)
       sheathtemperature(p)=mean(temp(time_inside_sheath), /NAN)
       sheathtemperaturestd(p)=stddev(temp(time_inside_sheath), /NAN)   
       
       sheathpdyn(p)=mean(dynamic_pressure(den(time_inside_sheath),vtot(time_inside_sheath)), /NAN)
       sheathpdynstd(p)=stddev(dynamic_pressure(den(time_inside_sheath),vtot(time_inside_sheath)), /NAN)
       
           
     endif else begin
      sheathspeed(p)=!VALUES.F_NAN
      sheathspeedstd(p)=!VALUES.F_NAN
      sheathdensity(p)=!VALUES.F_NAN
      sheathdensitystd(p)=!VALUES.F_NAN
      sheathtemperature(p)=!VALUES.F_NAN
      sheathtemperaturestd(p)=!VALUES.F_NAN       
      sheathpdyn(p)=!VALUES.F_NAN
      sheathpdynstd(p)=!VALUES.F_NAN  
     
     endelse
     
   
     ;all parameters inside MO
     time_inside_mo=where(timeb gt start_mav(p) and timeb lt end_mav(p))       
     if time_inside_mo(0) gt -1 then begin          
      bmax(p)=max(btot(time_inside_mo), /NAN)
      bmean(p)=mean(btot(time_inside_mo), /NAN)
      bstd(p)=stddev(btot(time_inside_mo), /NAN)
      bzmean(p)=mean(bz(time_inside_mo), /NAN)
      bzmin(p)=min(bz(time_inside_mo), /NAN)
      
      mospeed(p)=mean(vtot(time_inside_mo), /NAN)
      mospeedstd(p)=stddev(vtot(time_inside_mo), /NAN)
      
      modensity(p)=mean(den(time_inside_mo), /NAN)
      modensitystd(p)=stddev(den(time_inside_mo), /NAN)
      motemperature(p)=mean(temp(time_inside_mo), /NAN)
      motemperaturestd(p)=stddev(temp(time_inside_mo), /NAN)   

      mopdyn(p)=mean(dynamic_pressure(den(time_inside_mo),vtot(time_inside_mo)), /NAN)
      mopdynstd(p)=stddev(dynamic_pressure(den(time_inside_mo),vtot(time_inside_mo)), /NAN)



      duration(p)=(end_mav(p)-start_mav(p))/3600.
      heliodistance(p)=mean(helior(time_inside_mo), /NAN) /AU
      longitude(p)=mean(longit(time_inside_mo), /NAN)*!RADEG   
      latitude(p)=mean(latit(time_inside_mo), /NAN)*!RADEG 

      print, anytim(start_mav(p), /vms)
      print, 'mav/ICME max/mean/std btot:', bmax(p), bmean(p), bstd(p)
      print, 'duration', duration(p)
      
      ;get MVA results for magnetic obstacle
      mvub_results=mvub(bx(time_inside_mo),by(time_inside_mo),bz(time_inside_mo),timeb(time_inside_mo))     
  
      axis_longitude(p)=mvub_results(0)
      axis_latitude(p)=mvub_results(1)
      eigenvalue_ratio_2_3(p)=mvub_results(2)
      
      ;plot MO interval +/ 1 day      
      time_for_plot=where(timeb gt start_mav(p)-24*60*60. and timeb lt end_mav(p)+24*60*60.) 

      ;make in situ plots
      if make_plots_on gt 0 then r=make_plot(timeb(time_for_plot),btot(time_for_plot), $
             bx(time_for_plot), by(time_for_plot), bz(time_for_plot), $
             bmax(p),start_mav(p), end_mav(p),'mav',icmecat(mavicmes(p)).id,$
             lid,axis_longitude(p), axis_latitude(p),eigenvalue_ratio_2_3(p), ctype, icme_start_mav(p),icme_end_mav(p))

      
     endif           
    endfor     

  
 
; write results to master parameter variables
mbmax(mavicmes)=bmax
mbmean(mavicmes)=bmean
mbstd(mavicmes)=bstd
mbzmean(mavicmes)=bzmean
mbzmin(mavicmes)=bzmin

mmo_duration(mavicmes)=duration
mheliodistance(mavicmes)=heliodistance
mlongitude(mavicmes)=longitude
mlatitude(mavicmes)=latitude
;MVA results skipped for MAVEN
;maxis_longitude(mavicmes)=axis_longitude
;maxis_latitude(mavicmes)=axis_latitude
;meigenvalue_ratio_2_3(mavicmes)=eigenvalue_ratio_2_3

maxis_longitude(mavicmes)=!VALUES.F_NAN
maxis_latitude(mavicmes)=!VALUES.F_NAN
meigenvalue_ratio_2_3(mavicmes)=!VALUES.F_NAN

    
;speed results
mmospeed(mavicmes)=mospeed
mmospeedstd(mavicmes)=mospeedstd
msheathspeed(mavicmes)=sheathspeed
msheathspeedstd(mavicmes)=sheathspeedstd

;plasma results

mmodensity(mavicmes)=modensity
mmodensitystd(mavicmes)=modensitystd
mmotemperature(mavicmes)=motemperature
mmotemperaturestd(mavicmes)=motemperaturestd

msheathdensity(mavicmes)=sheathdensity
msheathdensitystd(mavicmes)=sheathdensitystd
msheathtemperature(mavicmes)=sheathtemperature
msheathtemperaturestd(mavicmes)=sheathtemperaturestd



msheathpdyn(mavicmes)=sheathpdyn
msheathpdynstd(mavicmes)=sheathpdynstd
mmopdyn(mavicmes)=mopdyn
mmopdynstd(mavicmes)=mopdynstd










    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
  
   
  
  

;add parameters to structure
struct_add_field, icmecat, 'sc_heliodistance', mheliodistance
struct_add_field, icmecat, 'sc_long_heeq', mlongitude
struct_add_field, icmecat, 'sc_lat_heeq', mlatitude

struct_add_field, icmecat, 'mo_bmax', mbmax
struct_add_field, icmecat, 'mo_bmean', mbmean
struct_add_field, icmecat, 'mo_bstd', mbstd
struct_add_field, icmecat, 'mo_bzmean', mbzmean
struct_add_field, icmecat, 'mo_bzmin', mbzmin
struct_add_field, icmecat, 'mo_duration', mmo_duration
struct_add_field, icmecat, 'mo_mva_axis_long', maxis_longitude
struct_add_field, icmecat, 'mo_mva_axis_lat', maxis_latitude
struct_add_field, icmecat, 'mo_mva_ratio', meigenvalue_ratio_2_3



struct_add_field, icmecat, 'sheath_speed', msheathspeed
struct_add_field, icmecat, 'sheath_speed_std', msheathspeedstd
struct_add_field, icmecat, 'mo_speed', mmospeed
struct_add_field, icmecat, 'mo_speed_std', mmospeedstd

struct_add_field, icmecat, 'sheath_density', msheathdensity
struct_add_field, icmecat, 'sheath_density_std', msheathdensitystd
struct_add_field, icmecat, 'mo_density', mmodensity
struct_add_field, icmecat, 'mo_density_std', mmodensitystd

struct_add_field, icmecat, 'sheath_temperature', msheathtemperature
struct_add_field, icmecat, 'sheath_temperature_std', msheathtemperaturestd
struct_add_field, icmecat, 'mo_temperature', mmotemperature
struct_add_field, icmecat, 'mo_temperature_std', mmotemperaturestd

struct_add_field, icmecat, 'sheath_pdyn', msheathpdyn
struct_add_field, icmecat, 'sheath_pdyn_std', msheathpdynstd
struct_add_field, icmecat, 'mo_pdyn', mmopdyn
struct_add_field, icmecat, 'mo_pdyn_std', mmopdynstd



;------------------ write icmecat to ASCII with header and comments


full_ascii_icmecat_filename='/home/cmoestl/helcats/products/ICMECAT/HELCATS_ICMECAT_v20_'+ctype+'.txt'
openw,lun,full_ascii_icmecat_filename, /get_lun
		for i=0,n_elements(icmecat)-1 do begin
					printf, lun, icmecat[i].id, icmecat[i].sc_insitu,	icmecat[i].icme_start_time, $
					             icmecat[i].mo_start_time, icmecat[i].mo_end_time, $
					             icmecat[i].icme_end_time, icmecat[i].mo_bmax, icmecat[i].mo_bmean, icmecat[i].mo_bstd, icmecat[i].mo_bzmean, icmecat[i].mo_bzmin, icmecat[i].mo_duration,$
					             icmecat[i].sc_heliodistance, icmecat[i].sc_long_heeq, icmecat[i].sc_lat_heeq, icmecat[i].mo_mva_axis_long,icmecat[i].mo_mva_axis_lat,icmecat[i].mo_mva_ratio, $
					             icmecat[i].sheath_speed, icmecat[i].sheath_speed_std, icmecat[i].mo_speed, icmecat[i].mo_speed_std, $
					             icmecat[i].sheath_density, icmecat[i].sheath_density_std, icmecat[i].mo_density, icmecat[i].mo_density_std, $
					             icmecat[i].sheath_temperature, icmecat[i].sheath_temperature_std, icmecat[i].mo_temperature, icmecat[i].mo_temperature_std,icmecat[i].sheath_pdyn, icmecat[i].sheath_pdyn_std,icmecat[i].mo_pdyn,  icmecat[i].mo_pdyn_std,$
					             FORMAT='(A32, A16, A20, A20, A20, A20, F10.1, F10.1, F10.1, F10.1,F10.1, F10.1, F10.4, F10.2, F10.2, F10.1, F10.1, F10.1,F10.1, F10.1, F10.1,F10.1,F10.1, F10.1, F10.1,F10.1,F10.1, F10.1, F10.1,F10.1,F10.1, F10.1, F10.1,F10.1)'
		endfor
	free_lun,lun
print, 'ASCII file created:      ', full_ascii_icmecat_filename

comments=write_header_comments(full_ascii_icmecat_filename, icme_observatories,masterfile,icmecat, ctype)

icmecatfilename='/home/cmoestl/helcats/products/ICMECAT/HELCATS_ICMECAT_v20_'+ctype+'.sav'
;save to ICMECAT products
save, icmecat, comments, filename=icmecatfilename, /verbose
print, 'SAV file with comments written to products/ICMECAT folder'
;read in product:
restore, icmecatfilename, /verbose
help, icmecat, /st


print, 'total number of events:  ',n_elements(icmecat.sc_insitu)


print, 'number of events for Wind, STA, STB, VEX, MES, ULY, MAV'
print, n_elements(where(icmecat.sc_insitu eq 'Wind'))
print, n_elements(where(icmecat.sc_insitu eq 'STEREO-A'))
print, n_elements(where(icmecat.sc_insitu eq 'STEREO-B'))
print, n_elements(where(icmecat.sc_insitu eq 'VEX'))
print, n_elements(where(icmecat.sc_insitu eq 'MESSENGER'))
print, n_elements(where(icmecat.sc_insitu eq 'ULYSSES'))
print, n_elements(where(icmecat.sc_insitu eq 'MAVEN'))

end