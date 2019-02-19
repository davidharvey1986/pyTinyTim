pro tinytim_rerrg, n_stars,		    $
                   X=x, 		    $
                   Y=y, 		    $
                   FOCUS_RANGE=focus_range, $
                   PIXEL_SCALE=pixel_scale, $
                   OUTPUT_DIR=output_dir,   $
                   SAME_STAR=same_star,     $
                   OVER=over,		    $
                   PLOTIT=plotit,	    $
                   SIMAGE=simage,	    $
                   MOMS=moms,		    $
                   DECIMAL=decimal,	    $
                   RRG=rrg,		    $
                   SHAPELETS=shapelets,     $
                   RAW=raw,		    $
                   VERTEX=vertex,	    $
                   EXACT_POSITION=exact_position

;+
; NAME:
;      TINYTIM_RERRG
;
; CATEGORY:
;      ACS data reduction
;
; PURPOSE:
;      Rerun RRG (& shapelets) on a set of TinyTim model stars.
;
; INPUTS:
;      As tinytim_create.pro.
;
; OPTIONAL INPUTS:
;      As tinytim_create.pro.
;
; KEYWORD PARAMETERS:
;      As tinytim_create.pro.
;
; OUTPUTS:
;      RRG (& shapelet) catalogue(s).
;
; MODIFICATION HISTORY:
;      Sep 05 - Written by Richard Massey.
;-

; Parse input options
compile_opt idl2
acspix=[4096,4096]                                                   ; Number of pixels in an ACS CCD.
if not keyword_set(n_stars) then n_stars=60                          ; Number of stars along one side of the grid (the TinyTim images are 34x34 pixels).
n_stars=fix(n_stars) & n_stars=n_stars+(n_stars mod 2)               ; Make sure this is an even integer, so it can be split easily along the chip boundary.
if not keyword_set(over) then over=5                                 ; Oversampling factor for Tiny Tim output - used internally only.
over=fix(over) & if over lt 2 then message,"Not oversampling enough!"; Must oversample by an integer!
if not keyword_set(focus_range) then focus_range=[0,0]               ; Range of focues values to be computed during this process.
if n_elements(focus_range) lt 2 then focus_range=[1,1]*focus_range   ; Cope with single values of focus.
if not keyword_set(output_dir) then output_dir="small_pixels"        ; Where to write the output file(s).
if not file_test(output_dir,/directory) then file_mkdir,output_dir   ; Make sure there's somewhere for the data to go
if not file_test(output_dir,/directory) then message,"Output directory does not exist and cannot be created!"
if not keyword_set(pixel_scale) then pixel_scale=0.03                                             ; Size of final pixels in arcseconds (only used if drizzling)
if keyword_set(raw) then pixel_scale=0.05                            ; Doesn't make sense otherwise 
if keyword_set(raw) then simagepix=acspix else simagepix=fix([4600,4600]*0.05/pixel_scale) ; Number of pixels in the simulated stargrid.
;if keyword_set(same_star) then message,"SAME_STAR doesn't work!"     ; Warning that I'd never quite got this to do what I wanted
if pixel_scale eq 0.05 then rrg_min_rad=3.5 $
  else if pixel_scale eq 0.03 then rrg_min_rad=6.0 $
  else rrg_min_rad=3.5*0.05/pixel_scale
if pixel_scale eq 0.05 then beta=2.0 $
  else if pixel_scale eq 0.03 then beta=3.0 $
  else beta=2.0*0.05/pixel_scale
n_max=20
spawn,"hostname",hostname                                            ; Obtain name of processor.
spawn,"echo $$",unique & unique=".temp_"+unique                      ; Obtain a unique number from process ID, with which to prepend all temporary files.
spawn,"setenv TINYTIM /usr/local/optical/tinytim" ; (This needs to be run from outside IDL - kept here for reference only)
if not keyword_set(decimal) then begin
  decimal=0
  decimal_string=strtrim(string(0),2)
endif else if decimal lt 1 then begin
  decimal_string=strtrim(string(fix(decimal*1000)),2)
endif else if alog10(decimal) le 1 then begin
  decimal_string=strtrim(string(fix(decimal)),2)
  decimal=decimal/10.
endif else begin
  decimal_string=strtrim(string(fix(decimal)),2)
  decimal=float(decimal)/(10^(fix(alog10(decimal)+1)))
endelse



; Read in previously-created image 
for focus=focus_range[0],focus_range[1] do begin

  ; Read in previously-created image array
  filename=output_dir+path_sep()+"TinyTim_f"+strtrim(string(focus),2)
  if not keyword_set(simage) then fits_read,filename+".fits",simage
  n_pixels_simage=size(simage,/DIMENSIONS)

  ; Find positions of stars
  if not file_test(filename+"_noisy.fits") then begin
    noisy_simage=simage
    noisy_simage+=randomn(seed,n_pixels_simage[0],n_pixels_simage[1])/(max(simage)*1e5)
    writefits,filename+"_noisy.fits",noisy_simage
    delvarx,noisy_simage
  endif
  if not file_test(filename+"_noisy.sexcat") then begin
    sex,filename+"_noisy.fits",/FULL_PATH,telescope="TinyTim"
  endif
  shapelets_read_sexcat,sexcat,filename+"_noisy.fits",/FULL_PATH
  xdist=reform(sexcat.x[*,0])
  ydist=reform(sexcat.x[*,1])
  n_stars=sexcat.n


  ; Create some empty variables to extend into the RRG catalogue
  if keyword_set(rrg) then begin
    xpos=0. & ypos=0. &   xx=0. &   yy=0. &   xy=0. 
    xxxx=0. & xxxy=0. & xxyy=0. & xyyy=0. & yyyy=0. 
  endif


  ; Loop over stars
  for i=0,n_stars-1 do begin
      
    ; Progress report
    print,"Focus position "+strtrim(string(focus),2)+"."+decimal_string+"um, "+$
    	  "star #"+strtrim(string(i+1),2)+"/"+strtrim(string(n_stars),2)+"."
      
    ; Cut out a postage stamp
    pstamp_size=105*pixel_scale/0.03
    left   = round(xdist[i]-pstamp_size/2.-1)>0
    right  = round(xdist[i]+pstamp_size/2.  )<(n_pixels_simage[0]-1)
    bottom = round(ydist[i]-pstamp_size/2.-1)>0
    top    = round(ydist[i]+pstamp_size/2.  )<(n_pixels_simage[1]-1)
    pstamp = simage[left:right,bottom:top]
    ;window,0 & shapelets_plot_image,alog10(pstamp>1e-6),/colbar,/frame,cran=cran


    ; Measure moments using RRG
    if keyword_set(rrg) or keyword_set(shapelets) then begin
      ;print,'Measuring shape moments using RRG'
      rrg_pstamp=(strtrim(string(unique),2)+"_"+strtrim(string(i),2)+".fits")[0]
      rrg_cat	=(strtrim(string(unique),2)+"_"+strtrim(string(i),2)+".cat")[0]
      rrg_mom	=(strtrim(string(unique),2)+"_"+strtrim(string(i),2)+".mom")[0]
      rrg_log	=(strtrim(string(unique),2)+"_"+strtrim(string(i),2)+".mom.log")[0]
      writefits,rrg_pstamp,pstamp
      spawn,"/usr/local/optical/sextractor/bin/sex "+rrg_pstamp+$
	    " -c rrg.sex -CATALOG_NAME "+rrg_cat+$
	    " -PARAMETERS_NAME rrg.param",junk,junk_error
      
      rrg_measure_moms,rrg_pstamp,rrg_cat[0],rrg_mom,min_rad=rrg_min_rad,/stellar
      restore,rrg_mom;,/verbose
      file_delete,rrg_pstamp
      file_delete,rrg_cat
      file_delete,rrg_mom
      if moms.xx ne -99 then begin
        file_delete,rrg_log
        xpos=[xpos,xdist[i]]
        ypos=[ypos,ydist[i]]
        xx=[xx,moms.xx[0]];(pixel_scale/0.05)^2]
        xy=[xy,moms.xy[0]];(pixel_scale/0.05)^2]
        yy=[yy,moms.yy[0]];(pixel_scale/0.05)^2]
        xxxx=[xxxx,moms.xxxx[0]];(pixel_scale/0.05)^4]
        xxxy=[xxxy,moms.xxxy[0]];(pixel_scale/0.05)^4]
        xxyy=[xxyy,moms.xxyy[0]];(pixel_scale/0.05)^4]
        xyyy=[xyyy,moms.xyyy[0]];(pixel_scale/0.05)^4]
        yyyy=[yyyy,moms.yyyy[0]];(pixel_scale/0.05)^4]
      endif 
    endif
   
    
    ; Measure shapelet coefficients
    if keyword_set(shapelets) then begin
      decomp=shapelets_decomp(pstamp,beta,n_max,centre=[moms.x,moms.y])
      ;decomp=shapelets_focus(pstamp,beta_guess=beta,n_max_range=[2,8]);,centre=grid_centre)
      decomp.x=[xdist[i],ydist[i]]
      if i eq 0 then begin
        shapecat={name:"TinyTim model stars at HST focus="+strtrim(string(focus),2),$
        	  type:"shapecat",			       $
        	  n:1L, 				       $
		  maxn_max:n_max,			       $
		  n_coeffs:long(decomp.n_coeffs),	       $
		  polar:0B,				       $
		  x:[[xdist[i]],[ydist[i]]],                   $
		  beta:beta,				       $
		  n_max:n_max,  			       $
		  coeffs:transpose(decomp.coeffs),	       $
		  coeffs_error:transpose(decomp.coeffs_error), $
		  flag:intarr(1,2),			       $
		  flag_interpret:["OK","OK"],		       $
		  chisq:decomp.chisq[1],		       $
        	  sextractor:0B,			       $
        	  morphology:0B,			       $
        	  shear_estimates:0B}
      endif else shapelets_add,shapecat,decomp
    endif

  endfor


  ; Clear simulated image
  delvarx,simage


  ; Write simulated RRG catalogue
  if keyword_set(rrg) then begin
    print,"Writing RRG catalogue to file "+filename+".moms"
    moms={x:xpos[1:*],                            $
          y:ypos[1:*],                            $
          xx:xx[1:*],                             $
          yy:yy[1:*],                             $
          xy:xy[1:*],                             $
          xxxx:xxxx[1:*],                         $
          xxxy:xxxy[1:*],                         $
          xxyy:xxyy[1:*],                         $
          xyyy:xyyy[1:*],                         $
          yyyy:yyyy[1:*],                         $
          e1:((xx-yy)/(xx+yy))[1:*],              $
          e2:(2*xy/(xx+yy))[1:*]     }
    save,filename=filename+".moms",moms
  endif


  ; Write simulated shapelet catalogue
  if keyword_set(shapelets) then begin
    print,"Writing shapelet catalogue to file "+filename+".shapecat"
    shapelets_write_shapecat,shapecat,filename,/FULL_PATH
  endif


endfor



end
