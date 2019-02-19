PRO tinytim_rrg, FITS_IMAGE, SEX_CATALOG, OUTFILE,$
                 WIDTH=width,$
                 SATURATION=saturation,$
   	             BADVAL=badval,$
                 CUT_OFF=cut_off,$
                 FIELDBACK=fieldback,$
                 MULT=mult,$
                 MIN_RAD=min_rad,$
	               STARTX=startx,$
                 STARTY=starty,$
                 SILENT=silent,$
                 STELLAR=stellar,$
	               WEIGHT=weight,$
                 WT_EXT=wt_ext,$
                 NOCENTER=nocenter,$
	               NAME_WEIGHT_IM=name_weight_im,$
                 GAIN=gain,$
                 EXP_TIME=exp_time,$
	               SGM_IM=sgm_im,$
                 ASCII_SEXCAT=ascii_sexcat
;
; NAME:
;    TINYTIM_RRG (a dedicated version of RRG_MEASURE_MOMS)
;
; PURPOSE:
;    Uses the RRG method to measure the moments of objects in a fits image.
;
; INPUTS:
;   fits_image- fits image with single extension
;   sex_catalog-SExtractor catalog in fits format
;   outfile- output file to store IDL array of variables
;
; OUTPUTS:
;   stores moments and other info in the output file as an idl structure  
;
; KEYWORD PARAMETERS:
;   ascii_sexcat - Switch to assume that SExtractor catalogue is in ASCII format (default is FITS_LDAC)
;   width- gaussian window function width, if not input is calculated from area
;   saturation,bad_value are saturation level and value of pixels not to use (outside image)
;   badval
;   cut_off- cut_off for measuring objects default is 2.5x width
;   fieldback -1 means use a field backgroup as calculated by IDL SKY, 0(default) uses object back from SEx
;   mult- multiplier used to find gaussian width from area
;   min_rad- stellar radius
;   startx, starty- lower left corner of image in pixel values, default to zero
;   silent set to 1 to supress messages
;   stellar set to 1 to evaluate all objects with the same radius (min_rad)
;   weight=use a weight image to determine good/bad pixels (yes=1 0=no/default)
;   wt_ext= extension of weight image in fits file
;  nocenter= if 1 then uses Sextractor centers and doesnt centroid
;  name_weight_im=name of the weight image
;  sgm_im Sextractor segmentation image to get sky background 
;
;MODIFICATION HISTORY:
;    30 Dec 2002 jrhodes
;    11 March 2003 jrhodes
;    02 January 2005 Ability to read in ASCII SExtractor catalogues added by RJM

	      
IF N_PARAMS(0) LT 1 THEN BEGIN
   PRINT, 'tinytim_rrg,fits_image,sex_catalog,outfile,width=width,saturation=saturation,'
   PRINT, '    badval=badval,cut_off=cut_off,fieldback=fieldback,mult=mult,min_rad=min_rad,'
   PRINT, '    startx=startx,starty=starty,silent=silent,stellar=stellar,'
   PRINT, '    weight=weight,wt_ext=wt_ext,nocenter=nocenter,name_weight_im=name_weight_im'
   return
ENDIF

; Read in SExtractor catalogue
cat=mrdfits(sex_catalog,2,hcat,/silent)
if size(cat,/TYPE) eq 2 then if cat[0] eq 0 then begin
  print,"Object not found!"
  moms={xx:-99,yy:-99}
  save,filename=outfile,moms
  return
endif
;
openw,1,outfile+'.log'
fits_read,fits_image,img,imhead
;sky,img,skymed,skysd
mmm,img,skymed,skysd
printf,1,skymed,skysd,' skymed and skysd'
num=n_elements(cat.x_image)
sexnum=intarr(num)
xin=dblarr(num) & yin=dblarr(num)
x=dblarr(num) & y=dblarr(num)
xx=dblarr(num) & yy=dblarr(num) & xy=dblarr(num)
xxxx=dblarr(num) & xxxy=dblarr(num) & xxyy=dblarr(num) & xyyy=dblarr(num) & yyyy=dblarr(num)
class=dblarr(num) & flags=dblarr(num) & back=dblarr(num)
area=dblarr(num) & x=dblarr(num) & y=dblarr(num)
xxerr=dblarr(num) & xyerr=dblarr(num) & yyerr=dblarr(num)
c_xx_yy=dblarr(num) &  c_xx_xy=dblarr(num) & c_xy_yy=dblarr(num)
xcerr=dblarr(num) & ycerr=dblarr(num) & det_area=dblarr(num)
x_sex=dblarr(num) & y_sex=dblarr(num)
theta_world=dblarr(num) & theta_image=dblarr(num)
x_sex=cat.x_image & y_sex=cat.y_image
theta_world=cat.theta_world & theta_image=cat.theta_image
fwhm=dblarr(num) & fwhm=cat.fwhm_image
ysize=sxpar(imhead,'NAXIS2') & xsize=sxpar(imhead,'NAXIS1')
kron=dblarr(num)
test=dblarr(num)
mag_iso=dblarr(num)
magerr_iso=dblarr(num)
mag_auto=dblarr(num)
magerr_auto=dblarr(num)
mu_max=dblarr(num)
mu_max=cat.mu_max
mag_iso=cat.mag_iso
mag_auto=cat.mag_auto
magerr_auto=cat.magerr_auto
magerr_iso=cat.magerr_iso
xin=cat.x_image & yin=cat.y_image
prob=intarr(num) ;lets you know if there is a centerprob (2), badpixl(4),edgeprob (1)
det_area=cat.isoarea_image
;area_ab=!pi*((cat.a_image+cat.b_image)/2)^2
area_ab=!pi*cat.a_image*cat.b_image
area=area_ab
class=cat.class_star
back=cat.background
flags=cat.flags
sn=cat.flux_auto/cat.fluxerr_auto
kron=cat.kron_radius
;print,weight,' weight'
IF NOT keyword_set(sgm_im) THEN sgm_im = 'null'
IF NOT keyword_set(fieldback) THEN fieldback = 0
;if fieldback then back=replicate(skymed,num) 
if fieldback then back=replicate(skymed,num) 
IF NOT keyword_set(cut_off) THEN cut_off=2.5 
IF NOT keyword_set(saturation) THEN saturation=160000.
IF NOT keyword_set(bad_val) THEN bad_val=-99
IF NOT keyword_set(mult) THEN mult = 1
IF NOT keyword_set(min_rad) THEN min_rad = 1.5
IF NOT keyword_set(startx) THEN startx = 0.
IF NOT keyword_set(starty) THEN starty = 0.
IF NOT keyword_set(silent) THEN silent = 0.
IF NOT keyword_set(stellar) THEN stellar = 0.
IF NOT keyword_set(weight) THEN weight = 0.
IF NOT keyword_set(wt_ext) THEN wt_ext = 0.
IF NOT keyword_set(nocenter) THEN nocenter = 0 
IF NOT keyword_set(name_weight_im) THEN name_weight_im =fits_image 
IF NOT keyword_set(gain) THEN gain=1
IF NOT keyword_set(exp_time) THEN exp_time=1
if (weight eq 1) then begin
    ;print,'weight image ',name_weight_im
    fits_read,name_weight_im,wt_image,exten_no=wt_ext
endif else begin
    wt_image=replicate(1,xsize,ysize)
endelse    
if ((sgm_im ne 'null')and(weight eq 1)and(fieldback eq 1)) then begin
    fits_read,sgm_im,seg,shdr
    sel=where((seg eq 0)and (wt_image gt 0))
    mmm,img(sel),seg_sky,seg_skysd
    back=replicate(seg_sky,num) 
endif ;sgm
radius=mult*sqrt(area/!pi)
for i=0,(num-1) do begin
    if (radius(i) lt min_rad) then radius(i) = min_rad
    if (stellar eq 1) then radius(i) = min_rad
endfor
IF keyword_set(width) then radius=replicate(width,num)
cut_rad=radius*cut_off ; cut off radius of object
offedge=0;number of objects off the endge initially
centerprob=0; number of objects with centroiding problems
badpix_prob=0;number of objects with bad pixels
;

;openw,1,'tmp.txt'
;
;cycle through objects
;sky,img,sky_mode,sky_sigma
sky_mode=skymed
sky_sigma=skysd
dI=sky_sigma
gain=1
;print,center,' center'
for i=0,(num-1) do begin
    sexnum(i)=i
;following changed by jrhodes to account for different indexing in SExtractor and IDL
    xc=xin(i)-0.5
    yc=yin(i)-0.5
    ;print,xin(i),yin(i),radius(i),back(i)
    deltax=1.0
    deltay=1.0
    count=1
    sum_int=0;
    xx_int=0;
    yy_int=0;
    xy_int=0;
    x_int=0;
    y_int=0;
    blank=0;
    go_on=1;
    xy_test=0;
    badpix_centroid='no'

    if ( ((xc-cut_rad(i)-1) lt startx) or ((yc-cut_rad(i)-1) lt starty) or $
	((xc+cut_rad(i)+1) gt xsize)or  ((yc+cut_rad(i)+1) gt ysize) ) then  begin
		;if not silent then printf,1,i,' too close to edge at (iteration,xin,yin) ',(count-1),xin(i),yin(i)
		if not silent then printf,1,i,xin(i),yin(i),' too close to edge at iteration 1'
		offedge=offedge+1;
		go_on=0;
		prob(i)=prob(i)+1
    endif ;xc
    
;find the centroid
    if (nocenter eq 0) then begin
    while  ( ( (abs(deltax) gt 0.01) or (abs(deltay) gt 0.01)) and (count lt 500) and (go_on eq 1)) do begin 
    ;while  ( ( (abs(deltax) gt 0.1) or (abs(deltay) gt 0.1)) and (count lt 500) and (go_on eq 1)) do begin
     ;while  ( ( (abs(deltax) gt 1.0) or (abs(deltay) gt 1.0)) and (count lt 500) and (go_on eq 1)) do begin
     	checkx=xc
        checky=yc
         
    	for county=round((yc-cut_rad(i)-1)),round((yc+cut_rad(i)+1)) do begin
    	for countx=round((xc-cut_rad(i)-1)),round((xc+cut_rad(i)+1)) do begin
		;dist=(sqrt((countx-xc)*(countx-xc)+(county-yc)*(county-yc)))
		dist=sqrt( (fix(countx)+0.5-xc)^2 + (fix(county)+0.5-yc)^2)
		if (dist le cut_rad(i))  then begin	
		    	if (wt_image(countx,county) eq 0.) then begin 
			badpix_centroid='yes'
			;print,countx,county,' badpix'
			endif else begin	
			    g_f=exp(-(dist*dist)/(2*radius(i)*radius(i)))
			    sum_int=sum_int+g_f*(img(countx,county)-back(i))
			    ;x_int=x_int+g_f*countx*(img(countx,county)-back(i))
			    ;y_int=y_int+g_f*county*(img(countx,county)-back(i))
			    x_int=x_int+g_f*(fix(countx)+0.5)*(img(countx,county)-back(i))
    	    	    	    y_int=y_int+g_f*(fix(county)+0.5)*(img(countx,county)-back(i))
			endelse
		endif ;dist
        endfor
	endfor
	;print,i,count,xc,yc,g_f,sum_int,x_int
	xc=x_int/sum_int;
	yc=y_int/sum_int;
	sum_int=0;
	x_int=0;
	y_int=0;
	sum_int=0;	
	count=count+1
	deltax=xc-checkx
	deltay=yc-checky;
    	;print,xc,yc,deltax,deltay,count
    

	if ( ((xc-cut_rad(i)-1) lt startx) or ((yc-cut_rad(i)-1) lt starty) or $
	((xc+cut_rad(i)+1) gt xsize)or  ((yc+cut_rad(i)+1) gt ysize) ) then  begin
			if not silent then printf,1,i,xin(i),yin(i),' too close to edge at iteration ',(count-1)
                        go_on=0
			prob(i)=prob(i)+1
	endif ;xc
	if (badpix_centroid eq 'yes') then begin
	    go_on=0
	    if not silent then printf,1,i,xin(i),yin(i),' bad pix in centroid at iteration ',(count-1)
	    prob(i)=prob(i)+4
	    badpix_prob=badpix_prob+1
	endif; badpix_centroid		
    endwhile
    
    
    if (count gt 100) then begin 
    	centerprob=centerprob+1
    	if not silent then printf,1,i,xin(i),yin(i),' Too many centering iterations ',count
	prob(i)=prob(i)+2
    endif ;count

    if ( sqrt(  (xc-xin(i))^2+ (yc-yin(i))^2) gt mult*radius(i))	then begin
    	go_on=0
	centerprob=centerprob+1
    	if not silent then printf,1,i,xin(i),yin(i),' Centroid shift too large '
	prob(i)=prob(i)+2
    endif; sqrt
    endif ;center
badpix_mom='no'    
; check for saturated pixels and bad pixels
    if (go_on eq 1) then begin
	for county=round((yc-cut_rad(i)-1)),round((yc+cut_rad(i)+1)) do begin
    	for countx=round((xc-cut_rad(i)-1)),round((xc+cut_rad(i)+1)) do begin
	;for county=(yc-cut_rad(i)-1),(yc+cut_rad(i)+1) do begin
    	;for countx=(xc-cut_rad(i)-1),(xc+cut_rad(i)+1) do begin
	    dist=(sqrt((countx-xc)*(countx-xc)+(county-yc)*(county-yc)));
    	    if (dist le cut_rad(i)) then begin		
    	    	if  ((img(countx,county) gt saturation) or (img(countx,county) eq bad_val) or $ 
		(wt_image(countx,county) eq 0)) then begin
		go_on=0
		;centerprob=centerprob+1.
		prob(i)=prob(i)+4
		badpix_mom='yes'
		;print,i,countx,county,img(countx,county),wt_image(countx,county),saturation,bad_val
    	    endif;  img
	    endif;  dist			
	endfor
	endfor
    endif ;go_on	
if (badpix_mom eq 'yes') then begin
    badpix_prob=badpix_prob+1 
     ;print,i,badpix_prob
     printf,1,i,xin(i),yin(i),' Bad pixel(s) in centroiding'
endif

;find moments and center error
    ;xc=(xc+0.5)
    ;yc=(yc+0.5)
    npix=0.
    npp=0
    nmm=0
    npm=0
    nmp=0
    if (go_on eq 1) then begin
	sum_ixxxx=0 & sum_ixxxy=0 & sum_ixxyy=0 &sum_ixyyy=0 & sum_iyyyy=0
	sum_xxdI=0 & sum_xydI=0 & sum_yydI=0
        ;for county=round(yc-cut_rad(i)-1),round(yc+cut_rad(i)+1) do begin
    	;for countx=round(xc-cut_rad(i)-1),round(xc+cut_rad(i)+1) do begin
	for county=(yc-cut_rad(i)-1),(yc+cut_rad(i)+1) do begin
    	for countx=(xc-cut_rad(i)-1),(xc+cut_rad(i)+1) do begin	 
	    distx=fix(countx)+0.5-xc
	    disty=fix(county)+0.5-yc
	    ;above two lines changed jrhodes 10-2003
	    ;also changed several lines down to measure distx, disty
	    dist=sqrt( (fix(countx)+0.5-xc)^2 + (fix(county)+0.5-yc)^2)
    	    ;dist=(sqrt((countx-xc)^2+(county-yc)^2))
	    if (dist le cut_rad(i)) then begin		
			npix=npix+1.
			g_f=exp(-(dist*dist)/(2*radius(i)*radius(i)));
			
			;xx_int=xx_int+g_f*(countx-xc)*(countx-xc)*(img(countx,county)-back(i));
			;yy_int=yy_int+ g_f*(county-yc)*(county-yc)*(img(countx,county)-back(i));
			;xy_int=xy_int+g_f*(countx-xc)*(county-yc)*(img(countx,county)-back(i));	
		    	;I_err=sqrt( dI*dI+((img(countx,county)-back(i)))/(gain*exp_time))
    	    	    	;sum_xxdI=sum_xxdI+g_f*g_f*I_err*I_err*(countx-xc)*(countx-xc)
    	    	    	;sum_yydI=sum_yydI+g_f*g_f*I_err*I_err*(county-yc)*(county-yc)
    	    	    	;sum_xydI=sum_xydI+g_f*g_f*I_err*I_err*(county-yc)*(countx-xc)
			
			xx_int=xx_int+g_f*(distx)*(distx)*(img(countx,county)-back(i));
			yy_int=yy_int+ g_f*(disty)*(disty)*(img(countx,county)-back(i));
			xy_int=xy_int+g_f*(distx)*(disty)*(img(countx,county)-back(i));	
		    	I_err=sqrt( dI*dI+((img(countx,county)-back(i)))/(gain*exp_time))
    	    	    	sum_xxdI=sum_xxdI+g_f*g_f*I_err*I_err*(distx)*(distx)
    	    	    	sum_yydI=sum_yydI+g_f*g_f*I_err*I_err*(disty)*(disty)
    	    	    	sum_xydI=sum_xydI+g_f*g_f*I_err*I_err*(disty)*(distx)
			sum_int=sum_int+ g_f*(img(countx,county)-back(i))
sum_ixxxx=sum_ixxxx+g_f*(img(countx,county)-back(i))* (distx)^4
sum_ixxxy=sum_ixxxy+g_f*(img(countx,county)-back(i))*((distx)^3)*(disty)
sum_ixxyy=sum_ixxyy+g_f*(img(countx,county)-back(i))*((distx)^2)* (disty)^2
sum_ixyyy=sum_ixyyy+g_f*(img(countx,county)-back(i))*((disty)^3)*(distx);
sum_iyyyy=sum_iyyyy+g_f*(img(countx,county)-back(i))*(disty)^4
    	    endif;  dist			
	endfor
	endfor
	ixx=xx_int/sum_int;
	iyy=yy_int/sum_int;
	ixy=xy_int/sum_int;
	ixxxx=sum_ixxxx/sum_int;
	ixxxy=sum_ixxxy/sum_int;
	ixxyy=sum_ixxyy/sum_int;
	ixyyy=sum_ixyyy/sum_int;
	iyyyy=sum_iyyyy/sum_int;
;print,ixx,iyy,sum_int
	sum_int=0;
    	xc_err=(sqrt(sum_xxdI))/sum_int
	yc_err=(sqrt(sum_yydI))/sum_int

	e1=(ixx-iyy)/(ixx+iyy);
	e2=(2*ixy)/(ixx+iyy);
	alpha=0.5*atan(e2,e1);

;find moment errors and covariances
	sum_xx_ixx=0;
	sum_yy_iyy=0;
	sum_xy_ixy=0;
 	sum_cov_xx_xy=0;
 	sum_cov_xx_yy=0;
 	sum_cov_xy_yy=0;
	;saturation=0
	
	for county=(yc-cut_rad(i)-1),(yc+cut_rad(i)+1) do begin
    	for countx=(xc-cut_rad(i)-1),(xc+cut_rad(i)+1) do begin	 
	    distx=countx-xc
	    disty=county-yc
	    dist=sqrt( (fix(countx)+0.5-xc)^2 + (fix(county)+0.5-yc)^2)
	    if (dist le cut_rad(i)) then begin		
	    	g_f=exp(-(dist*dist)/(2*radius(i)*radius(i)))
	    	sum_int=sum_int+ g_f*(img(countx,county)-back(i))
			;I_err=sqrt( dI*dI+((img(countx,county)-back(i)))/gain)
			I_err=sqrt( dI*dI+((img(countx,county)-back(i)))/(gain*exp_time))
    	    	    	;sum_xx_ixx+=g_f*g_f*I_err*I_err*pow( pow(countx-xc,2)-ixx,2);
    	    	    	;sum_yy_iyy+=g_f*g_f*I_err*I_err*pow( pow(county-yc,2)-iyy,2);
    	    	    	;sum_xy_ixy+=g_f*g_f*I_err*I_err*pow( (countx-xc)*(county-yc) -ixy,2);
    	    	    	;sum_cov_xx_yy+=g_f*g_f*((countx-xc)*(countx-xc)-ixx)*((county-yc)*(county-yc)-iyy);
    	    	    	;sum_cov_xx_xy+=g_f*g_f*((countx-xc)*(countx-xc)-ixx)*((county-yc)*(countx-xc)-ixy);
    	    	    	;sum_cov_xy_yy+=g_f*g_f*((countx-xc)*(county-yc)-ixy)*((county-yc)*(county-yc)-iyy);
			sum_xx_ixx=    sum_xx_ixx+    g_f*g_f*I_err*I_err*((countx-xc)^2-ixx)^2
    	    	    	sum_yy_iyy=    sum_yy_iyy+    g_f*g_f*I_err*I_err*((county-yc)^2-iyy)^2
    	    	    	sum_xy_ixy=    sum_xy_ixy+    g_f*g_f*I_err*I_err*((countx-xc)*(county-yc) -ixy)^2
    	    	    	sum_cov_xx_yy= sum_cov_xx_yy+ g_f*g_f*I_err*I_err*((countx-xc)*(countx-xc)-ixx)*((county-yc)*(county-yc)-iyy)
    	    	    	sum_cov_xx_xy= sum_cov_xx_xy+ g_f*g_f*I_err*I_err*((countx-xc)*(countx-xc)-ixx)*((county-yc)*(countx-xc)-ixy)
    	    	    	sum_cov_xy_yy= sum_cov_xy_yy+ g_f*g_f*I_err*I_err*((countx-xc)*(county-yc)-ixy)*((county-yc)*(county-yc)-iyy)

		endif ;dist			
	endfor ;county
	endfor;  countx 
	ixx_err=sqrt(sum_xx_ixx)/sum_int;
	iyy_err=sqrt(sum_yy_iyy)/sum_int;
	ixy_err=sqrt(sum_xy_ixy)/sum_int;
	covxx_yy=sum_cov_xx_yy/(sum_int*sum_int);
	covxx_xy=sum_cov_xx_xy/(sum_int*sum_int);
	covxy_yy=sum_cov_xy_yy/(sum_int*sum_int);
    	sum_int=0
endif else begin ;go on
    ixx=-99 & iyy=-99 & ixy=0 & ixxxx=0 & ixxxy=0
    ixxyy=0 & ixyyy=0 & iyyyy=0
    ixx_err=0 & iyy_err=0 & ixy_err=0
    covxx_yy=0 & covxx_xy=0 & covxy_yy=0
    xc_err=0 & yc_err=0
endelse
    ;test(i)=xy_test/npix
    x(i)=xc & y(i)=yc & xx(i)=ixx & yy(i)=iyy & xy(i)=ixy 
    xxxx(i)=ixxxx & xxxy(i)=ixxxy & xxyy(i)=ixxyy & xyyy(i)=ixyyy & yyyy(i)=iyyyy
    xxerr(i)=ixx_err & yyerr(i)=iyy_err & xyerr(i)=ixy_err
    c_xx_yy(i)=covxx_yy & c_xx_xy(i)=covxx_xy & c_xy_yy(i)=covxy_yy
    xcerr(i)=xc_err & ycerr(i)=yc_err
endfor
;area_ab=area
;area=cat.isoarea_image
a=cat.a_image
b=cat.b_image
theta=cat.theta_image
mag=cat.mag_auto
e1=(xx-yy)/(xx+yy)
e2=2*xy/(xx+yy)

;e1_err=sqrt((pow(ixx+iyy,-2))*(ixx_err*ixx_err*pow(1-e1,2)+iyy_err*iyy_err*pow(1+e1,2)-2*(1-e1*e1)*covxx_yy));
;e2_err=sqrt((pow(ixx+iyy,-2))*((ixx_err*ixx_err+iyy_err*iyy_err+2*covxx_yy)*e2*e2+4*(ixy_err*ixy_err-e2*(covxx_xy+covxy_yy))))
e1_err=sqrt( ((xx+yy)^(-2))*(xxerr*xxerr*(1-e1)^2 +yyerr*yyerr*(1+e1)^2-2*(1-e1*e1)*c_xx_yy))
e2_err=sqrt( ((xx+yy)^(-2))*((xxerr*xxerr+yyerr*yyerr+2*c_xx_yy)*e2*e2+4*(xyerr*xyerr-e2*(c_xx_xy*c_xy_yy))))


moms={sexnum:sexnum,x:x,y:y,xx:xx,yy:yy,xy:xy,xxxx:xxxx,xxxy:xxxy,$
xxyy:xxyy,xyyy:xyyy,yyyy:yyyy,radius:radius,sn:sn,back:back,$
class:class,area:area,area_ab:area_ab,det_area:det_area,flags:flags,$
a:a,b:b,theta:theta,mag:mag,prob:prob,test:test,$
xxerr:xxerr,yyerr:yyerr,xyerr:xyerr,$
c_xx_yy:c_xx_yy,c_xx_xy:c_xx_xy,c_xy_yy:c_xy_yy,$
e1:e1,e2:e2,e1_err:e1_err,e2_err:e2_err,$
xcerr:xcerr,ycerr:ycerr,kron:kron,fwhm:fwhm,mu_max:mu_max,$
mag_iso:mag_iso,magerr_iso:magerr_iso,mag_auto:mag_auto,magerr_auto:magerr_auto,$
x_sex:x_sex,y_sex:y_sex,theta_world:theta_world,theta_image:theta_image}
save,filename=outfile,moms
if offedge gt 0 or centerprob gt 0 or badpix_prob gt 0 then $
print,offedge,' are off edge and ',centerprob,' have center problems',badpix_prob,' have bad pix'
printf,1,offedge,' are off edge and ',centerprob,' have center problems',badpix_prob,' have bad pix'
close,1
END
