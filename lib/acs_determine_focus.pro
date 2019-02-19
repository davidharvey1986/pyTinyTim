function acs_determine_focus_metric,true,model,A=A,B=B
if n_elements(A) eq 0 then A=0
if n_elements(B) eq 0 then B=1
goodness_of_fit=total( sqrt(A^2*(true.xx+true.yy-model.xx-model.yy)^2)+B^2*((true.e1-model.e1)^2+(true.e2-model.e2)^2) ) /float((n_elements(model.e1)-1)>1)
return,goodness_of_fit
end

; **************************************************************************

pro acs_determine_focus, results, CATALOGUES=catalogues

;+      
; NAME:
;      ACS_DETERMINE_FOCUS
;
; CATEGORY:
;      Reduction of ACS COSMOS data.
;
; PURPOSE:
;      Decides the actual focus value of HST during observations, due to 
;      thermal fluctuations that change the distance between the primary
;      and secondary mirrors. The offset from the nominal focus position
;      is returned, in microns.
;
; INPUTS:
;      None.
;
; OPTIONAL INPUTS:
;      CATALOGUES - Structure containing model TinyTim PSF catalogues. Can
;                   be set for speed. If not set, ;read from disk on first run.
;
; OUTPUTS:
;      RESULTS - Structure containing a list of field names and focus offsets.
;
; KEYWORDS:
;      None.
;
; EXAMPLE USE:
;      acs_determine_focus, results
;
; MODIFICATION HISTORY:
;      Feb 05 - All shape moments needed by RRG used, by RM.
;      Jan 05 - Written by Richard Massey.
;-

;goto,time_series

; Initialise default variables
datadir="/scr/rhodes/Cosmos/Momstar/"
modeldir="/disk/cosmos/user/rjm/idl/TinyTim/grid_30x30/"
allstars=0B
stars_mag_ran =[19,23];-28
stars_fwhm_ran=[1.5,2.2]
sg_class_cut=0.95
e_cut=0.2
sn_cut=1
x_cut=[400,4200]
y_cut=[400,4200]
r_match=60
loadct,4,/silent
colours=[0,255,150]
n_grid=40                   ; Number of points in the whisker grid
shapelets_make_xarr,[n_grid,n_grid],grid_x,grid_y 
dist=acs_map_xy(4096.*grid_x/n_grid+2048,4096.*grid_y/n_grid+2048,offset=[199,210])
grid_x=dist.x & grid_y=dist.y

; Look for input files, containing RRG catalogues of objects run with /stellar flag
files=file_search(datadir,"acs*.momstar",count=n_files)
files=strmid(files,strlen(datadir))                                 ; Strip path
for i=0,n_files-1 do files[i]=strmid(files[i],0,strlen(files[i])-8) ; Strip extension
readcol,"cosmos_dates.txt",field,mjd,year,month,day,hour,format="A,F,I,I,I,I",skip=1,/silent
message,"Found "+strtrim(string(n_files),2)+" fields.",/info

; Create empty global variables in which to later put the answers
global_n_stars=intarr(n_files)
global_focus=fltarr(n_files)
global_focus_error=fltarr(n_files)
global_mjd=fltarr(n_files)
global_date=intarr(n_files,4)

; Loop over all fields
for i=0,n_files-1 do begin


  ; Load measured catalogue
  restore,datadir+files[i]+".momstar"


; *************
;moms.x=moms.x+200
;moms.y=moms.y+200
  
  ; Find date of exposure
  date_match=where(field eq files[i])
  if date_match[0] ne -1 then begin
    global_mjd[i]=mjd[date_match]
    global_date[i,*]=[hour[date_match],day[date_match],month[date_match],year[date_match]]
  endif
  
  
  ; Pick the useful stars in that image
  stars=where(moms.mag   gt stars_mag_ran[0]  $
          and moms.mag   lt stars_mag_ran[1]  $
          and moms.fwhm  gt stars_fwhm_ran[0] $
          and moms.fwhm  lt stars_fwhm_ran[1] $
          and moms.x     gt x_cut[0]          $
          and moms.x     lt x_cut[1]          $
          and moms.y     gt y_cut[0]          $
          and moms.y     lt y_cut[1]          $
          ;and moms.sn    lt sn_cut            $
          and moms.class gt sg_class_cut      $
          and (moms.e1^2+moms.e2^2) lt e_cut^2)
  n_objects=n_elements(moms.x)
  n_stars=n_elements(stars)
  ;plot,moms.mag,moms.fwhm,psym=3,yran=[0,10]

  
  ; Calculate model ellipticities
  model_e=acs_model_e(moms.x[stars], moms.y[stars], focus, catalogues=catalogues, datadir=modeldir)
  n_focus=n_elements(focus)
  average_distance=total(model_e.offset_model,1)/n_focus


  ; Select only those stars with a suitably close model (since we're not interpolating)
  close_match=where(average_distance le r_match,n_stars)
  global_n_stars[i]=n_stars
  if n_stars lt 2 then begin
    message,"No stars found in field "+files[i],/info
    if i gt 0 then begin
      global_focus[i]=global_focus_error[i-1]
      global_focus_error[i]=global_focus_error[i-1]  
    endif else begin
      global_focus[i]=0.
      global_focus_error[i]=100.
    endelse
    continue
  endif
  
  ; Tabulate model ellipticities
  model_e1=model_e.e1[*,close_match]
  model_e2=model_e.e2[*,close_match]

  ; Tabulate measured ellipticities
  true_e1=moms.e1[stars[close_match]]
  true_e2=moms.e2[stars[close_match]]
  true={e1:moms.e1[stars[close_match]],$
        e2:moms.e2[stars[close_match]],$
        x:moms.x[stars[close_match]],$
        y:moms.y[stars[close_match]],$
        xy:moms.xy[stars[close_match]],$
        xx:moms.xx[stars[close_match]],$
        yy:moms.yy[stars[close_match]],$
        xxxx:moms.xxxx[stars[close_match]],$
        xxxy:moms.xxxy[stars[close_match]],$
        xxyy:moms.xxyy[stars[close_match]],$
        xyyy:moms.xyyy[stars[close_match]],$
        yyyy:moms.yyyy[stars[close_match]]}

  ; Find best-fit focus
  chisq=fltarr(n_focus) ; Absolute chi squared
  gof=fltarr(n_focus)   ; Goodness-of-fit
  if allstars then begin
    for f=0,n_focus-1 do begin
      ; Tabulate model ellipticities
      model={e1:reform(model_e.e1[f,close_match]),$
             e2:reform(model_e.e2[f,close_match]),$
             x:reform(model_e.x[close_match]),$
             y:reform(model_e.y[close_match]),$
             xx:reform(model_e.xx[f,close_match]),$
             xy:reform(model_e.xy[f,close_match]),$
             yy:reform(model_e.yy[f,close_match]),$
             xxxx:reform(model_e.xxxx[f,close_match]),$
             xxxy:reform(model_e.xxxy[f,close_match]),$
             xxyy:reform(model_e.xxyy[f,close_match]),$
             xyyy:reform(model_e.xyyy[f,close_match]),$
             yyyy:reform(model_e.yyyy[f,close_match])}
      gof[f]=acs_determine_focus_metric(true,model)
      chisq[f]=gof[f]*50
    endfor
    ; Store the best one
    junk=min(chisq,best_fit_focus) 
    model={e1:reform(model_e.e1[best_fit_focus,close_match]),$
           e2:reform(model_e.e2[best_fit_focus,close_match]),$
           x:reform(model_e.x[close_match]),$
           y:reform(model_e.y[close_match]),$
           xx:reform(model_e.xx[best_fit_focus,close_match]),$
           xy:reform(model_e.xy[best_fit_focus,close_match]),$
           yy:reform(model_e.yy[best_fit_focus,close_match]),$
           xxxx:reform(model_e.xxxx[best_fit_focus,close_match]),$
           xxxy:reform(model_e.xxxy[best_fit_focus,close_match]),$
           xxyy:reform(model_e.xxyy[best_fit_focus,close_match]),$
           xyyy:reform(model_e.xyyy[best_fit_focus,close_match]),$
           yyyy:reform(model_e.yyyy[best_fit_focus,close_match])}
    best_fit_focus=focus[best_fit_focus]
    error_focus=0.
    
    ; Report findings to screen
    !p.multi=0
    ;oplot,focus,chisq,thick=4;,title="!6"+files[i],$
         ;xtitle="!6HST Focus position [!7l!6m]";,yran=[0,1]
    print,"model:",mean(model.xx)," truth:",mean(true.xx)
    ;read,temp
    ;plot,[0,0],/nodata,xtitle="!6Data",ytitle="!6Tiny Tim model",$
    ;     /iso,psym=4,xran=[1.5,2],yran=[1.5,2],/xstyle,/ystyle
    ;;oplot,(true.xx+true.yy)/2,(model.xx+model.yy)/2,psym=4
    ;oplot,true.xx,model.xx,psym=5,symsize=0.5
    ;oplot,true.yy,model.yy,psym=6,symsize=0.5
    ;;oplot,abs(true.xy)*10+1,abs(model.xy)*10+1,psym=1
    ;;oplot,[0,100],[0,100]
    ;read,temp
    !p.multi=[0,1,2]
    binsize=0.005
    ;plothist,true.xx,bin=binsize,/fill,xran=[1.5,2.2],yran=[0,76],xtitle="!6!8I!6!Dxx!N"
    ;;plothist,model_e.xx,bin=binsize/20.,/over,linestyle=0,color=150
    ;plothist,model.xx,bin=binsize,/over,linestyle=0,color=250
    ;plothist,true.yy,bin=binsize,/fill,xran=[1.5,2.2],yran=[0,76],xtitle="!6!8I!6!Dyy!N"
    ;;plothist,model_e.yy,bin=binsize/20.,/over,linestyle=0,color=150
    ;plothist,model.yy,bin=binsize,/over,linestyle=0,color=250
    ;read,temp
    ;plothist,true.xx+true.yy,bin=binsize,/fill,xran=[3,4.4],yran=[0,76],xtitle="!6!8I!6!Dxx!N+!8I!6!Dyy!N"
    ;;plothist,model_e.xx+model_e.yy,bin=binsize/20.,/over,linestyle=0,color=150
    ;plothist,model.xx+model.yy,bin=binsize,/over,linestyle=0,color=250
    ;plothist,true.xy,bin=binsize/4,/fill,xran=[-0.2,0.2],yran=[0,76],xtitle="!6!8I!6!Dxy!N"
    ;plothist,model.xy,bin=binsize/4,/over,linestyle=0,color=250
    ;read,temp
    !p.multi=0
  endif else begin
    chisq_indiv=fltarr(n_focus,n_stars)
    best_fit_focus_indiv=fltarr(n_stars)
    for s=0,n_stars-1 do begin
      for f=0,n_focus-1 do begin
        chisq_indiv[f,s]= $
          total(sqrt((true_e1[s]-model_e1[f,s])^2+(true_e2[s]-model_e2[f,s])^2))
        ; Tabulate ellipticities in real data
        true_star={e1:true.e1[s],$
                   e2:true.e2[s],$
                   x:true.x[s],$
                   y:true.y[s],$
                   xx:true.xx[s],$
                   xy:true.xy[s],$
                   yy:true.yy[s],$
                   xxxx:true.xxxx[s],$
                   xxxy:true.xxxy[s],$
                   xxyy:true.xxyy[s],$
                   xyyy:true.xyyy[s],$
                   yyyy:true.yyyy[s]}
        ; Tabulate model ellipticities
        model={e1:model_e.e1[f,close_match[s]],$
               e2:model_e.e2[f,close_match[s]],$
               x:model_e.x[close_match[s]],$
               y:model_e.y[close_match[s]],$
               xx:model_e.xx[f,close_match[s]],$
               xy:model_e.xy[f,close_match[s]],$
               yy:model_e.yy[f,close_match[s]],$
               xxxx:model_e.xxxx[f,close_match[s]],$
               xxxy:model_e.xxxy[f,close_match[s]],$
               xxyy:model_e.xxyy[f,close_match[s]],$
               xyyy:model_e.xyyy[f,close_match[s]],$
               yyyy:model_e.yyyy[f,close_match[s]]}
        ; Caluculate goodness of fit
        chisq_indiv[f,s]=acs_determine_focus_metric(true_star,model,A=A,B=B)
      endfor
      junk=min(chisq_indiv[*,s],temp) & best_fit_focus_indiv[s]=focus[temp]
    endfor
    best_fit_focus=median(best_fit_focus_indiv)
    if n_stars lt 2 then begin
      chisq=chisq_indiv
      error_focus=15 
    endif else begin
      chisq=total(chisq_indiv,2)/float(n_stars-1)
      error_focus=stddev(best_fit_focus_indiv);/sqrt(n_stars)
    endelse
      plot,focus,chisq,thick=4,title="!6"+files[i],$
           xtitle="!6HST Focus position [!7l!6m]",yran=[0,0.2]
      oplot,focus,chisq,thick=4
      oplot,[1,1]*best_fit_focus,[0,1e5],thick=4
    for s=0,n_stars-1 do oplot,focus,chisq_indiv[*,s]
    ; Report findings to screen
  endelse
  ; Save results in variables with a global scope
  print,files[i]," (",strtrim(string(i),2),"/",strtrim(string(n_files),2),")",$
        n_stars,best_fit_focus,error_focus
  global_focus[i]=best_fit_focus
  global_focus_error[i]=error_focus  

  ; Print higher-order shape moments
  model_moms=acs_model_e(moms.x[stars], moms.y[stars], round(best_fit_focus[0]), catalogues=catalogues)
;  print,mean(model_moms.xx[close_match]),mean(moms.xx[stars[close_match]])
;  print,mean(model_moms.xx[close_match]-model_moms.yy[close_match]),mean(moms.xx[stars[close_match]]-moms.yy[stars[close_match]])
;  print,2*mean(model_moms.xy[close_match]),2*mean(moms.xy[stars[close_match]])
;  print,mean(model_moms.xx[close_match]+model_moms.yy[close_match]),mean(moms.xx[stars[close_match]]+moms.yy[stars[close_match]])

  ; Make whisker ;plot
  ;model_moms=catalogues.(where(focus eq round(best_fit_focus[0]))+2)
  model_moms=acs_model_e(grid_x, grid_y, round(best_fit_focus[0]), catalogues=catalogues, datadir=modeldir)
  scale=1000
  plot,[0,0],/nodata,xran=[0,4650],yran=[0,4650],/xstyle,/ystyle,/iso,$
       title="!6"+files[i]+": focus "+strtrim(string(round(best_fit_focus[0])),2)+"!7l!6m",$
       xtitle="!6x [pixels]",ytitle="!6y [pixels]"
  plt_evec,model_moms.x,model_moms.y,model_moms.e1,model_moms.e2,/e1e2,xscale=scale,yscale=scale
  plt_evec,moms.x[stars],moms.y[stars],moms.e1[stars],moms.e2[stars],/e1e2,xscale=scale,yscale=scale,thick=2
  plt_evec,200,4400,0.1,0,/e1e2,xscale=scale,yscale=scale
  xyouts,200,4450,"!610%",align=0.5
  ;read,temp
  
endfor

; Save results to disc
results={field:files,                    $
         mjd:global_mjd,                 $
         date:global_date,               $
         focus:round(global_focus),      $
         focus_best_fit:global_focus,    $
         focus_error:global_focus_error, $
         focus_n_stars:global_n_stars}
save,results,filename="cosmos_focus.sav"


time_series:

restore,"cosmos_focus.sav"

t=results.mjd-min(results.mjd)
f=results.focus[*,0]

;plot,t,f,psym=7,symsize=0.3,xrange=[140,210],yrange=[-10,5],/xstyle,/ystyle;,/nodata
;slice_fit,t,f,100,2,xfit,yfit,ysigma,/plot;,/bars
x=findgen(1200)/5.
;oplot,x,-2.*cos(0.42*(x-146))-4.5;,/linestyle


ops,file="cosmos_focus.ps"
for i=0,1 do begin

  if i eq 0 then xran=[140,210] else if i eq 1 then xran=[0,215]
  
  usersym, 0.7*cos(2*!pi*findgen(21)/20), 0.7*sin(2*!pi*findgen(21)/20), /fill ; Set psym=8 to be a circle
  plot,[0,0],xran=xran,yran=[-10,5],/xstyle,/ystyle,/nodata,$
     title="!6Focus time series",$
     xtitle="!6Days after start of COSMOS observations",ytitle="!6HST focus offset [!7l!6m]"
  oploterr,results.mjd-min(results.mjd),$
           results.focus[*,0],$
           results.focus_error[*,0]/sqrt(results.focus_n_stars[*,0]),8
  
endfor
cps

end
