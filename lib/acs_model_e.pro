function acs_model_e, x, y, focus, catalogues=catalogues, datadir=datadir, dither=dither

;+
; NAME:
;       ACS_MODEL_E
;
; PURPOSE:
;       Obtain the ellipticity of a TinyTim model PSF at one
;       position and any focus value(s) on ACS.
;  
; CATEGORY:
;       ACS data reduction.
;
; CALLING SEQUENCE:
;       result=acs_model_e(x, y, focus, catalogues=catalogues)
;
; INPUTS:
;       X          - Object x coordinate.
;       Y          - Object y coordinate.
;
; OPTIONAL INPUTS:
;       FOCUS      - Focus offset from optimal [um].
;                    Can be one number or an array.
;       CATALOGUES - Array of model ellipticities. If not set (e.g. on first
;                    use), will read in from disk. If already set, will recall
;                    from memory, for speed.
;       DITHER     - Dither pattern (PSF is averaged from exposures at all of
;                    these positions)
;
; KEYWORD PARAMETERS:
;       None.
;
; OUTPUTS:
;       Returns a structure containing model ellipticities and other parameters.
;
; OPTIONAL OUTPUTS:
;       None.
;
; COMMON BLOCKS:
;       None.
;
; SIDE EFFECTS:
;       None.
;
; RESTRICTIONS:
;       None.
;
; EXAMPLE:
;       See use in acs_determine_focus.pro.
;
; MODIFICATION HISTORY:
;       Feb 05 - All shape moments needed by RRG now returned. RM.
;       Feb 05 - Written by Richard Massey.
;-
; 

if keyword_set(dither) then begin

  size_dither=size(dither,/DIMENSIONS)
  if size_dither[0] ne 2 then message,"DITHERs must be entered in [2,n_dither] format!"
  n_dither=size_dither[1]

  for i=0,n_dither-1 do begin
    single_exposure=acs_model_e(x+dither[0,i],y+dither[1,i],focus,catalogues=catalogues,datadir=datadir)
    if i eq 0 then begin
      answer     =single_exposure
      answer.xx  =answer.xx/n_dither
      answer.yy  =answer.yy/n_dither
      answer.xxxx=answer.xxxx/n_dither
      answer.xxxy=answer.xxxy/n_dither
      answer.xxyy=answer.xxyy/n_dither
      answer.xyyy=answer.xyyy/n_dither
      answer.yyyy=answer.yyyy/n_dither
    endif else begin
      answer.xx  =answer.xx+single_exposure.xx/n_dither
      answer.yy  =answer.yy+single_exposure.yy/n_dither
      answer.xxxx=answer.xxxx+single_exposure.xxxx/n_dither
      answer.xxxy=answer.xxxy+single_exposure.xxxy/n_dither
      answer.xxyy=answer.xxyy+single_exposure.xxyy/n_dither
      answer.xyyy=answer.xyyy+single_exposure.xyyy/n_dither
      answer.yyyy=answer.yyyy+single_exposure.yyyy/n_dither
    endelse 
  endfor
  answer.x=x
  answer.y=y

;  print,[416,412,408,402]-413.253
;  print,[3788,3688,3590,3484]-3645.08
; [[2.74701,142.920],[-1.25299,42.9199],[-5.25299,-55.0801],[-11.2530,-161.080]]





endif else begin

  ; Initialise variables
  ;if not keyword_set(datadir) then datadir="/raid1/rjm/idl/TinyTim/small_pixels/"
  if not keyword_set(datadir) then datadir="/raid1/rjm/idl/TinyTim/small_pixels_swap_nominal_dithered/"
  filenames1=datadir+"TinyTim_f" & filenames2=".moms"
  catalogues_focus=strmid(["f-10","f-9","f-8","f-7","f-6","f-5","f-4","f-3","f-2","f-1","f0","f1","f2","f3","f4","f5"],1)
  n_coords=n_elements(x)
  
  ; Load data (required on first use only)
  if not keyword_set(catalogues) then begin
    catalogues={focus:float(catalogues_focus),focus_string:catalogues_focus}
    for i=0,n_elements(catalogues.focus)-1 do begin
      restore,filenames1+catalogues_focus[i]+filenames2;,/verbose
      moms=create_struct(moms,"focus",catalogues.focus[i],"focus_string",catalogues.focus_string[i])
      catalogues=create_struct(catalogues,"moms"+strtrim(string(i),2),moms)
    endfor
  endif
  
  ; Which focus values are we interested in? (default: all of them)
  if n_elements(focus) eq 0 then focus=catalogues.focus
  n_focus=n_elements(focus)
  
  ; Prepare empty arrays to contain the answer
  e1=fltarr(n_focus,n_coords)
  e2=fltarr(n_focus,n_coords)
  xx=fltarr(n_focus,n_coords)
  xy=fltarr(n_focus,n_coords)
  yy=fltarr(n_focus,n_coords)
  xxxx=fltarr(n_focus,n_coords)
  xxxy=fltarr(n_focus,n_coords)
  xxyy=fltarr(n_focus,n_coords)
  xyyy=fltarr(n_focus,n_coords)
  yyyy=fltarr(n_focus,n_coords)
  x_model=fltarr(n_focus,n_coords)
  y_model=fltarr(n_focus,n_coords)
  offset_model=fltarr(n_focus,n_coords)

  ; Loop over each focus value in turn
  for f=0L,n_focus-1 do begin

    ; Isolate the RRG catalogue for that focus value
    moms=catalogues.(where(catalogues.focus eq focus[f])+2)

    ; Loop over each position in which we're interested
    for i=0L,n_coords-1 do begin
      ; Find closest model star (i.e. don't inerpolate)
      rmin=sqrt(min((moms.x-x[i])^2+(moms.y-y[i])^2,match))
      ;print,x[i],y[i],r,moms.x[match],moms.y[match]
      
      ; Store its ellipticity and shape moments
      e1[f,i]=moms.e1[match]
      e2[f,i]=moms.e2[match]
      xx[f,i]=moms.xx[match]
      xy[f,i]=moms.xy[match]
      yy[f,i]=moms.yy[match]
      xxxx[f,i]=moms.xxxx[match]
      xxxy[f,i]=moms.xxxy[match]
      xxyy[f,i]=moms.xxyy[match]
      xyyy[f,i]=moms.xyyy[match]
      yyyy[f,i]=moms.yyyy[match]
      
      ; Store its other quantities
      x_model[f,i]=moms.x[match]
      y_model[f,i]=moms.y[match]
      offset_model[f,i]=rmin
    endfor
    
  endfor
  
  ; Tell the world what we've found
  answer={x:x,y:y,focus:focus,$
          e1:e1,e2:e2,$
          xx:xx,xy:xy,yy:yy,$
          xxxx:xxxx,xxxy:xxxy,xxyy:xxyy,xyyy:xyyy,yyyy:yyyy,$
          offset_model:offset_model,x_model:x_model,y_model:y_model}

endelse

return,answer

end
