pro tinytim_create, n_stars,                 $
                    X=x,                     $
                    Y=y,                     $
                    FOCUS_RANGE=focus_range, $
                    PIXEL_SCALE=pixel_scale, $
                    OUTPUT_DIR=output_dir,   $
                    SAME_STAR=same_star,     $
                    WAVELENGTH=wavelength,   $
                    FILTER=filter,           $
                    OVER=over,               $
                    PLOTIT=plotit,           $
                    SIMAGE=simage,           $
                    MOMS=moms,               $
                    DECIMAL=decimal,         $
                    RRG=rrg,                 $
                    SHAPELETS=shapelets,     $
                    RAW=raw,                 $
                    VERTEX=vertex,           $
                    FLAT_CCD=flat_ccd,       $
                    DOUBLE_CONVOLUTION=double_convolution,$
                    EXACT_POSITION=exact_position, $
                    FITSNAME=FITSNAME,       $
                    COORDFILE=COORDFILE
;+      
; NAME:
;      TINYTIM_CREATE
;
; CATEGORY:
;      Mimicry of ACS (COSMOS) images.
;
; PURPOSE:
;      Makes big images containing lots of stars by repeatedly running
;      Tiny Tim, at various positions on the CCD and at various focus values.
;
; INPUTS:
;      N_STARS     - Number of stars in the grid. e.g. tinytim_create,30
;                    will create an image containing a 30x30 grid of stars.
;
; OPTIONAL INPUTS:
;      X and Y     - Positions of stars (two vectors of same length).
;      FOCUS_RANGE - Create images for several focus values, in integer
;                    values between focus_range[0] microns and 
;                    focus_range[1] microns.
;      OUTPUT_DIR  - Directory where files will be written out. Warning:
;                    file names do not change. If this directory already
;                    contains the output of a previous run, it will be
;                    overwritten!
;      PIXEL_SCALE - Size of final output pixel, in arcseconds. Default 0.05.
;      WAVELENGTH  - Positive values are used to indicate the wavelength for
;                    monochromatic observations, at which the PSF is to be
;                    calculated. Negative numbers indicate the TinyTim stellar
;                    spectral type e.g. -13 would give a G8V star.
;      FILTER      - Used in conjunction with stellar spectral types, this
;                    string should be recognised by TinyTim as an ACS filter.
;      SAME_STAR   - Create a grid of identical stars, by passing the
;                    same focus and CCD positions to Tiny Tim. Focus is
;                    specified using FOCUS_RANGE, and the others should be
;                    entered as SAME_STAR in the format [chip#,x,y].
;      DECIMAL     - Add this part-integer values to focus position(s).
;      OVER        - Oversampling factor. Determines the precision with which
;                    it is possible to locate a star within a pixel.
;                    Tint Tim performs oversampling only up to 10 times. Any
;                    more is obtained from that with cubic interpolation.
;                    To obtain positions accurate to within 1/100th of a
;                    pixel, set over=20 and the /exact_position flag.
;
;      FITSNAME     -NAME OF THE OUTPUT FITS
; KEYWORDS:
;      RAW         - Creates simulated images of ACS data before
;                    the application of DRIZZLE and its associated 
;                    correction for astrometric distortions in the camera.
;      EXACT_POSIT - Place stars at the exact position calculated or input.
;                    Otherwise, they are shifted to the centre of the nearest
;                    pixel.
;      VERTEX      - Place stars at the vertex of four adjacent pixels instead.
;      FLAT_CCD    - Remove TinyTim model of CCD height variation.
;      RRG         - Also run RRG and output a moms catalogue of the stars.
;      SHAPELETS   - Also perform a shapelet decomposition and output a
;                    shapecat catalogue of the stars.
;
; OUTPUTS:
;      Simulated starfield image written to disc.
;
; OPTIONAL OUTPUTS:
;      SIMAGE      - Simulated image returned in this variable.
;      MOMS        - RRG moment catalogue returned in this variable.
;
; UNWANTED OUTPUTS:
;      Temporary files are written to disc, called .temp*. These are removed
;      automatically at the end. But if the routine is interrupted, they will
;      clutter up disc space. Their names feature a random number, so that
;      many instances of this routine can be run simultaneously without
;      overwriting these temporary files.
;
; EXAMPLE USE:
;      setenv TINYTIM /usr/local/optical/tinytim
;      idl
;      IDL> tinytim_create, 120, focus=[-10,5], /rrg
;      IDL> tinytim_create,pixel_scale=0.03,over=10,output_dir="small_pixels",/rrg,/shapelets,focus=-9

; MODIFICATION HISTORY:
;      Dec 15 - Changed so that focus_range is now a list not range DH
;      Oct 05 - Chip height variations removed by RM.
;      Sep 05 - Size of weight function for RRG/shapelets generalised by RM.
;      Aug 05 - PIXEL_SCALE option finalised by RM.
;      Apr 05 - PIXEL_SCALE option debugged by RM. It is still highly untested.
;      Mar 05 - X, Y, OVER and SAME_STAR optional inputs added by RM.
;      Mar 05 - Order of chips swapped over in (raw) _flt files by RM.
;      Mar 05 - Creation of raw (still distorted images) reintroduced by RM.
;      Feb 05 - DRIZZLE-like convolution with final pixel added by RM.
;      Jan 05 - Stars moved to post-distortion positions by RM.
;      Jan 05 - RRG output added by Jason Rhodes.
;      Jan 05 - Focus difference between CCDs incorporated by RM.
;      Dec 04 - Bug fixed in positions of stars passed to Tiny Tim by RM.
;      Dec 04 - Written by Richard Massey.
;-
codeDir = '/Users/DavidHarvey/Library/Code/IDL/tinytim/IDL'
if keyword_set(coordfile) then begin
     print, "READING IN COORDINATE FILE: ",coordfile
     readcol, coordfile, x, y, $
              format='D,D', /silent
     print,"READ ",n_elements(x)," COORDINATES"
     if n_elements(x) eq 0 then return
endif


; Parse input options
compile_opt idl2
acspix=[4096,4096]                                                   ; Number of pixels in an ACS CCD.
if not keyword_set(n_stars) then n_stars=fix((acspix[0]-68.)/68.)    ; Number of stars along one side of the grid (the TinyTim images are 34x34 pixels).
n_stars=fix(n_stars) & n_stars=n_stars+(n_stars mod 2)               ; Make sure this is an even integer, so it can be split easily along the chip boundary.
if not keyword_set(over) then over=5                                 ; Oversampling factor for Tiny Tim output - used internally only.
over=fix(over) & if over lt 2 then message,"Not oversampling enough!" ; Must oversample by an integer!

if not keyword_set(filter) then filter="f814w"                       ; Filter used to take observation (only used with stellar spectral types)
if not keyword_set(focus_range) then focus_range=[0,0]               ; Range of focues values to be computed during this process.
if n_elements(focus_range) lt 2 then focus_range=[1,1]*focus_range   ; Cope with single values of focus.
if not keyword_set(output_dir) then output_dir="." ; Where to write the output file(s).
if not file_test(output_dir,/directory) then file_mkdir,output_dir   ; Make sure there's somewhere for the data to go
if not file_test(output_dir,/directory) then message,"Output directory does not exist and cannot be created!"
if not keyword_set(pixel_scale) then pixel_scale=0.05                ; Size of final pixels in arcseconds (only used if drizzling)
if keyword_set(raw) then pixel_scale=0.05                            ; Doesn't make sense otherwise 
if keyword_set(raw) then simagepix=acspix else simagepix=fix([4600,4600]*0.05/pixel_scale) ; Number of pixels in the simulated stargrid.
;if keyword_set(same_star) then message,"SAME_STAR doesn't work!"    ; Warning that I'd never quite got this to do what I wanted
if pixel_scale eq 0.05 then rrg_min_rad=3.5 $
  else if pixel_scale eq 0.03 then rrg_min_rad=6.0 $
  else rrg_min_rad=3.5*0.05/pixel_sacle
if pixel_scale eq 0.05 then beta=2.0 else beta=2.0*0.05/pixel_scale
n_max=20
spawn,"hostname",hostname                                            ; Obtain name of processor.
spawn,"echo $$",unique & unique=".temp_"+unique                      ; Obtain a unique number from process ID, with which to prepend all temporary files.
;spawn,"echo $$",unique & unique="~/data/.temp/"+unique               ; Obtain a unique number from process ID, with which to prepend all temporary files.
;spawn,"setenv TINYTIM /usr/local/optical/tinytim" ; (This needs to be run from outside IDL - kept here for reference only)
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



; Tidy up directory to prepare for future temporary files
file_delete,unique+"_runme",/allow_nonexistent    
file_delete,unique+"_parameter_file",/allow_nonexistent    
file_delete,unique+"_parameter_file_focus",/allow_nonexistent    
file_delete,unique+"_psf_image00.fits",/allow_nonexistent 
file_delete,unique+"_psf_image00_psf.fits",/allow_nonexistent 
file_delete,unique+"_psf_image.tt3",/allow_nonexistent 



; Decide star positions
if keyword_set(x) and keyword_set(y) then begin
  ; Check that input positions lie within the bounds of the CCDs
  x=0>x<(acspix[0]-1) & y=0>y<(acspix[1]-1)
  ; Count the stars!
  n_stars=n_elements(x)<n_elements(y)
  n_j=1
  
endif else begin
  ; Make an artificial grid of stars
  shapelets_make_xarr, [n_stars,n_stars], x, y
  x=x*(acspix[0])/((n_stars)>1)+acspix[0]/2
  y=y*(acspix[1])/((n_stars)>1)+acspix[1]/2
  n_j=n_stars
endelse
if keyword_set(raw) then begin
  ; Just imaging the raw CCD image
  ydist=y
  xdist=x
endif else begin
  ; Determine where the stars move, after correction for optical distortions
  xdist=acs_map_xy(x,y,pixel_scale=pixel_scale,offset=drizzle_offset)
  ydist=xdist.y
  xdist=xdist.x
endelse


; Make stars
for ifocus=0,n_elements(focus_range)-1 do begin
   focus = focus_range[ifocus]
  ; Create an empty array to contain the simulated image
  simage=dblarr(simagepix[0],simagepix[1])

  ; Create some empty variables to extend into the RRG catalogue
  if keyword_set(rrg) then begin
    n=n_stars*n_stars
    xpos=0. & ypos=0. &   xx=0. &   yy=0. &   xy=0. 
    xxxx=0. & xxxy=0. & xxyy=0. & xyyy=0. & yyyy=0. 
  endif

  for i=0,n_stars-1 do begin
    for j=0,n_j-1 do begin
      ; Locate star within raw image
      if y[i,j] gt 2048. then chip=1 else chip=2

      
      ; Progress report
      print,"Focus position "+strtrim(string(focus),2)+"."+decimal_string+"um, "+$
           ;"CCD chip #"+strtrim(string(chip),2)+", "+$
            "star #"+strtrim(string(i*n_j+j+1),2)+"/"+strtrim(string(n_stars*n_j),2)+$
            ", running on "+hostname+"."
            
      ; Run Tiny Tim (tiny1)
      openw,lun,unique+"_runme",/get_lun
      printf,lun,"/Users/DavidHarvey/Library/Code/IDL/tinytim/tinytim/tiny1 "+unique+"_parameter_file << EOF" ; Rin Tiny1
      printf,lun,strtrim(string(15),2)                    ; Select ACS
      printf,lun,strtrim(string(chip),2)                  ; Select which CCD chip
      printf,lun,strtrim(string(fix(x[i,j])),2)           ; x-position of star
      printf,lun,strtrim(string(fix(y[i,j] mod 2048.)),2) ; y-position of star
      if wavelength lt 0 then begin
        printf,lun,filter                                 ; Select filter
        printf,lun,strtrim(string(1),2)                   ; Select stellar spectrum 
        printf,lun,strtrim(string(-wavelength),2)         ; 13 would give a G8V star
      endif else begin
        printf,lun,"MONO"                                 ; Faster to use only one wavelength 
        printf,lun,strtrim(string(wavelength),2)          ; Wavelength [nanometres]
      endelse
      printf,lun,strtrim(string(3.0),2)                   ; Diameter of psf image
      printf,lun,unique+"_psf_image"                      ; Select output filename
      printf,lun,"EOF"
      close,lun
      free_lun,lun
      
      spawn,"source "+unique+"_runme",junk
      file_delete,unique+"_runme"
      
      
      ; Test we're using a known version of Tiny Tim
      if not stregex("*Tiny Tim v6.3*",junk[1]) then begin
        message,"Using unknown version of Tiny Tim. Wanted v6.3, but got:",/info
        print,junk[0:8]
        stop
      endif
      
      
      ; Adjust focus position
      ;spawn,"nedit "+unique+"_parameter_file"
      wcl=file_lines(unique+"_parameter_file") & wcl=wcl[0]
      openr,lun_in,unique+"_parameter_file",/get_lun
      openw,lun_out,unique+"_parameter_file_focus",/get_lun
      temp_string=""
      focus_zerothorder=0
      for line=1,wcl do begin
        readf,lun_in,temp_string
        if keyword_set(same_star) and strmatch(temp_string,"# Camera ID number*") then begin
          printf,lun_out,strtrim(string(24),2)+"  # Camera ID number"
        endif else if keyword_set(same_star) and strmatch(temp_string,"# Position 1*") then begin
          printf,lun_out,strtrim(string(1024),2)+" "+strtrim(string(512),2)+"  # Position 1"
        endif else if strmatch(temp_string,"# Focus*") and keyword_set(flat_ccd) then begin
          printf,lun_out,temp_string
          printf,lun_out,"#"
          printf,lun_out,0.00,0.00,0.00,0.00,0.00,0.00
          printf,lun_out,0.00,0.00,0.00,0.00,0.00,0.00
          printf,lun_out,0.00,0.00,0.00,0.00,0.00,0.00
          printf,lun_out,0.00,0.00,0.00,0.00,0.00,0.00
          printf,lun_out,0.00,0.00,0.00,0.00,0.00,0.00
          printf,lun_out,0.00,0.00,0.00,0.00,0.00,0.00
          line=line+7
          for l=1,2 do readf,lun_in,temp_string
          ;focus_zerothorder+=float((strsplit(temp_string,/EXTRACT))[0])
          for l=1,5 do readf,lun_in,temp_string
        endif else if strmatch(temp_string,"*Z4 = Focus for center of ACS/WFC field*") then begin
          focus_zeropoint=float((strsplit(temp_string,/EXTRACT))[0])
          printf,lun_out,focus_zeropoint+focus_zerothorder+(focus+decimal)*0.011," # Z4 = Focus for centre of ACS/WFC field"
          ;if focus eq 0 and decimal eq 0 then begin
          ;  printf,lun_out,temp_string
          ;endif else begin
          ;  case chip of
          ;    1: printf,lun_out,0.028+focus_zerothorder+(focus+decimal)*0.011," # Z4 = Focus for centre of ACS/WFC field"
          ;    2: printf,lun_out,0.048+focus_zerothorder+(focus+decimal)*0.011," # Z4 = Focus for centre of ACS/WFC field"
          ;  endcase
          ;endelse
        endif else if keyword_set(raw) then begin
          printf,lun_out,temp_string
        endif else if strmatch(temp_string,"# X,Y -> V2 transform coefficients*") then begin
          printf,lun_out,temp_string
          printf,lun_out,"#"
          printf,lun_out,0.00,0.05
          printf,lun_out,0.00,0.00,0.00
          printf,lun_out,0.00,0.00,0.00,0.00
          printf,lun_out,0.00,0.00,0.00,0.00,0.00
          line=line+5 & for l=1,5 do readf,lun_in,temp_string
        endif else if strmatch(temp_string,"# X,Y -> V3 transform coefficients*") then begin
          printf,lun_out,temp_string
          printf,lun_out,"#"
          printf,lun_out,0.05,0.00
          printf,lun_out,0.00,0.00,0.00
          printf,lun_out,0.00,0.00,0.00,0.00
          printf,lun_out,0.00,0.00,0.00,0.00,0.00
          line=line+5 & for l=1,5 do readf,lun_in,temp_string
        endif else if strmatch(temp_string,"# V2,V3 -> X*") then begin
          printf,lun_out,temp_string
          printf,lun_out,"#"
          printf,lun_out,0.00,20.0
          printf,lun_out,0.00,0.00,0.00
          printf,lun_out,0.00,0.00,0.00,0.00
          printf,lun_out,0.00,0.00,0.00,0.00,0.00
          line=line+5 & for l=1,5 do readf,lun_in,temp_string
        endif else if strmatch(temp_string,"# V2,V3 -> Y*") then begin
          printf,lun_out,temp_string
          printf,lun_out,"#"
          printf,lun_out,20.0,0.00
          printf,lun_out,0.00,0.00,0.00
          printf,lun_out,0.00,0.00,0.00,0.00
          printf,lun_out,0.00,0.00,0.00,0.00,0.00
          line=line+5 & for l=1,5 do readf,lun_in,temp_string
        endif else begin
          printf,lun_out,temp_string
        endelse
      endfor
      close,lun_in
      free_lun,lun_in
      close,lun_out
      free_lun,lun_out
      file_delete,unique+"_parameter_file"    
      ;spawn,"nedit "+unique+"_parameter_file_focus"


      ; Run the rest of Tiny Tim (tiny2 and tiny3)
      spawn,"/Users/DavidHarvey/Library/code/IDL/tinytim/tinytim/tiny2 "+unique+"_parameter_file_focus",junk   ; This creates a 134x134 pixel image unique+"_psf_image00_psf.fits", with 0.023054 arcsec/pixel.
      spawn,"/Users/DavidHarvey/Library/code/IDL/tinytim/tinytim/tiny3 "+unique+"_parameter_file_focus sub="+$
            strtrim(string(over<10),2),junk                ; With over=5, this creates a 371x371 pixel image, with 0.01 arcsec/pixel.


      ; Read in the postage stamp image and delete temporary files
      fits_read,unique+"_psf_image00.fits",pstamp,pstamp_header
      if keyword_set(plotit) then plt_image,alog10(float(pstamp)>1e-8),/fr,/col
      file_delete,unique+"_parameter_file_focus"
      file_delete,unique+"_psf_image.tt3"
      file_delete,unique+"_psf_image00_psf.fits"
      file_delete,unique+"_psf_image00.fits"


      ; Oversample more than Tiny Tim's 10x maximum, if necessary by cubic
      ; interpolation
      if over gt 10 then begin
        pstamp_size=round(size(pstamp,/DIM)*over/10.)
        if (over mod 10) eq 0 then begin
          pstamp=rebin(pstamp,pstamp_size[0],pstamp_size[1])
        endif else begin
          pstamp=congrid(pstamp,pstamp_size[0],pstamp_size[1],/center,/cubic)
        endelse
        if keyword_set(plotit) then plt_image,alog10(float(pstamp)>1e-8),/fr,/col
      end


      ; Obtain charge diffusion kernel
      n_header=size(pstamp_header,/DIMENSIONS)
      kernel=[[float((strsplit(pstamp_header[n_header-4],/extract))[1:3])],$
              [float((strsplit(pstamp_header[n_header-3],/extract))[1:3])],$
              [float((strsplit(pstamp_header[n_header-2],/extract))[1:3])]]

 
      ; Distort charge diffusion kernel to mimic correction for astrometric
      ; distortions
      if keyword_set(raw) then begin
        ; If we haven't distorted the PSF, catch up with calculation of the
        ; star's sub-pixel location
        pstamp_size=size(pstamp,/DIM)
        max_pstamp=max(pstamp,grid_centre)
        grid_centre=[grid_centre mod pstamp_size[0],grid_centre/pstamp_size[0]]
      endif else begin
        ; Determine linear approximation to astrometric distortion, locally
        case chip of
          1: begin
               ; Coefficients obtained from *_flt_coeffs2.dat, output by DRIZZLE
               xcoeffs=[2048-56.69626992,0.99648395,0.039235942,8.443148e-06,-4.9350914e-06,1.8714837e-06,-4.7000765e-10,-6.3905282e-11,-5.2871972e-10,-5.8953743e-12,2.4971845e-14,4.9862677e-15,3.637933e-14,-1.8702128e-14,1.405061e-14]
               ycoeffs=[2048-1030.25854442,0.029775126,1.0053584,-1.4995461e-06,6.0549338e-06,-7.2055625e-06,7.1958641e-11,-5.1743686e-10,-8.3696034e-11,-4.2292829e-10,-1.7454337e-14,-1.9344815e-15,-2.8682819e-14,1.3391103e-14,1.0777943e-15]
             end
          2: begin
               ; Coefficients obtained from *_flt_coeffs2.dat, output by DRIZZLE
               xcoeffs=[2048+28.69390174,0.98437314,0.045339576,8.2503964e-06,-7.1088858e-06,1.7728154e-06,-4.5976284e-10,-1.2408655e-10,-5.29678e-10,4.2119219e-11,1.7393012e-14,3.2388574e-14,9.2331066e-14,-2.0957229e-14,1.3739822e-13]
               ycoeffs=[2048+1047.39000180,0.042595268,0.97152006,-2.51594e-06,5.9160766e-06,-9.3900933e-06,6.7942709e-11,-4.3467827e-10,-5.3234262e-11,-3.9818745e-10,-1.6672274e-14,-3.5575914e-14,-1.0733564e-13,2.7255259e-14,-1.680405e-13] 
             end
        endcase
        xpower=[0,1,0,2,1,0,3,2,1,0,4,3,2,1,0]
        ypower=[0,0,1,0,1,2,0,1,2,3,0,1,2,3,4]
        ; Find distortion matrix
        Atli=where(xpower eq 1)
        Atl=total((xcoeffs*y[i,j]^ypower)[where(xpower eq 1)])
        Atr=total((xcoeffs*x[i,j]^xpower)[where(ypower eq 1)])
        Abl=total((ycoeffs*y[i,j]^ypower)[where(xpower eq 1)])
        Abr=total((ycoeffs*x[i,j]^xpower)[where(ypower eq 1)])
        A=[[Atl,Atr],[Abl,Abr]]
        ; Find inverse distortion matrix
        DetA=(Atl*Abr)-(Abl*Atr)
        AItl=Abr/DetA
        AItr=-Abl/DetA
        AIbl=-Atr/DetA
        AIbr=Atl/DetA
        AI=[[AItl,AItr],[AIbl,AIbr]]


        ; Crop zeros around edge of postage stamp, for speed later on
        if keyword_set(plotit) then plt_image,alog10(float(pstamp)>1e-8),/fr,/col
        print,"Cropping"
        cropx=where(total(pstamp,2) ne 0) & cropx=[min(cropx),max(cropx)]
        cropy=where(total(pstamp,1) ne 0) & cropy=[min(cropy),max(cropy)]
        pstamp=pstamp[cropx[0]:cropx[1],cropy[0]:cropy[1]]
        if keyword_set(plotit) then plt_image,alog10(float(pstamp)>1e-8),/fr,/col
        
        
        ; Calculate original (distorted) pixels
        print,"Pixellating in distorted frame"
        pstamp_size=size(pstamp,/DIM)
        grid=fltarr(pstamp_size[0],pstamp_size[1],2)
        max_pstamp=max(pstamp,grid_centre)
        grid_centre=[grid_centre mod pstamp_size[0],grid_centre/pstamp_size[0]]
        for pi=0,pstamp_size[0]-1 do grid[pi,*,0]=(pi-grid_centre[0]+1)/float(over)
        for pi=0,pstamp_size[1]-1 do grid[*,pi,1]=(pi-grid_centre[1]+1)/float(over)
        grid_dist=fltarr(pstamp_size[0],pstamp_size[1],2)
        for oi=0,pstamp_size[0]-1 do begin
          for oj=0,pstamp_size[1]-1 do begin
            grid_position=0>round(([AItl*oi+AItr*oj,AIbl*oi+AIbr*oj]))<[pstamp_size[0]-1,pstamp_size[1]-1]
            grid_dist[oi,oj,*]=grid[grid_position[0],grid_position[1],*]
          endfor
        endfor
        grid_dist[*,*,0]=grid_dist[*,*,0]+0.5 ; Not sure why the 0.5 is necessary, but it does get the star to the centre of [0,0] pixel
        if keyword_set(vertex) then begin
          grid_dist=round(grid_dist+0.5)
        endif else begin
          grid_dist=round(grid_dist)
        endelse
        
        
        ; Pixellate in original (distorted) pixels
        grid_dist_xrange=[min(grid_dist[*,*,0]),max(grid_dist[*,*,0])]
        grid_dist_yrange=[min(grid_dist[*,*,1]),max(grid_dist[*,*,1])]
        small_pstamp_size=[grid_dist_xrange[1]-grid_dist_xrange[0]+1,$
                           grid_dist_yrange[1]-grid_dist_yrange[0]+1]
        small_pstamp=fltarr(small_pstamp_size)
        for pi=grid_dist_xrange[0],grid_dist_xrange[1] do begin
          for pj=grid_dist_yrange[0],grid_dist_yrange[1] do begin
            this_pixel=where(grid_dist[*,*,0] eq pi and $
                             grid_dist[*,*,1] eq pj,n_this)
            if n_this gt 0 then $
              small_pstamp[pi-grid_dist_xrange[0],pj-grid_dist_yrange[0]]=mean(pstamp[this_pixel]) ;else print,pi,pj
          endfor
        endfor
        if keyword_set(plotit) then plt_image,alog10(float(small_pstamp)>1e-8),/fr,/col


        ; Convolve image with charge diffusion kernel
        print,"Convolving with charge diffusion kernel"
        small_pstamp=convol(small_pstamp,kernel/total(kernel),/center,/edge_truncate)
        if keyword_set(plotit) then plt_image,alog10(float(small_pstamp)>1e-8),/fr,/col


        ;Put back into oversampled pstamp
        print,"Oversampling again"
        for pi=grid_dist_xrange[0],grid_dist_xrange[1] do begin
          for pj=grid_dist_yrange[0],grid_dist_yrange[1] do begin
             this_pixel=where(grid_dist[*,*,0] eq pi and $
                              grid_dist[*,*,1] eq pj,n_this)
             if n_this gt 0 then pstamp[this_pixel]=$
                small_pstamp[pi-grid_dist_xrange[0],pj-grid_dist_yrange[0]]
          endfor
        endfor
        if keyword_set(plotit) then plt_image,alog10(float(pstamp)>1e-8),/fr,/col

      
        ; Convolve with RAW ACS pixel size
        ; REALLY NOT SURE THAT THIS IS A GOOD IDEA...
        if keyword_set(double_convolution) then begin
        ; print,"Convolving with initial pixel to mimic data acquisition(?)"
        ; pixel=fltarr(over*0.05/pixel_scale,over*0.05/pixel_scale)+(1./(over*0.05/pixel_scale)^2)
          print,"Convolving with final pixel to mimic DRIZZLE(?)"
          pixel=fltarr(over,over)+(1./over^2)
          pstamp=convol(pstamp,pixel,/center,/edge_truncate)
        endif
        

        ; Convolve with final pixel size (to mimic the effect of DRIZZLE)
        ; NOT SURE IF I SHOULD DO THIS! WELL, THE STARS ARE TOO SMALL WITHOUT IT...
        print,"Convolving with final pixel to mimic DRIZZLE(?)"
        pixel=fltarr(over,over)+(1./over^2)
        pstamp=convol(pstamp,pixel,/center,/edge_truncate)

 
     endelse
      
      
      ; Shift around a bit to get the centre of the star to the desired place 
      ; within a final pixel
      print,"Shifting to desired place within a pixel"
      grid_centre=(grid_centre+0.5)/float(over*0.05/pixel_scale)
      ;print,"grid_centre",grid_centre
      delta_pixel=0.05/pixel_scale/over
      ;print,"delta_pixel",delta_pixel
      if keyword_set(exact_position) then begin
        desired_position=[xdist[i,j],ydist[i,j]] mod 1
      endif else begin
        if keyword_set(vertex) then begin
          desired_position = [0. ,0. ] ; Aim for corner of pixel
        endif else begin
          desired_position = [0.5,0.5] ; Aim for centre of pixel
        endelse
      endelse
      ;print,"desired_position",desired_position
      shift=(over-round(((grid_centre-desired_position) mod 1)/delta_pixel))>0


      ; Add rows/columns of pixels to bottom/left to move centre of star
      pstamp_size=size(pstamp,/DIM)
      if shift[0] gt 0 then pstamp=[fltarr(shift[0],pstamp_size[1]),pstamp]
      if shift[1] gt 0 then pstamp=[[fltarr(pstamp_size[0]+shift[0],shift[1])],$
                                    [pstamp]]
      grid_centre=grid_centre/delta_pixel+shift
      ;print,"grid_centre",grid_centre
      
      ; Add rows/columns of pixels to top/right to leave an integer number x(over)
      if abs(abs(2*((over*pixel_scale/0.05) mod 1)-1)-1) gt 1e-4 then message,"Interpolation won't work correctly for this amount of oversampling!"
      pstamp_size=size(pstamp,/DIM)
      pstamp_temp=fltarr(pstamp_size+round(over*pixel_scale/0.05)-$
                  (pstamp_size mod round(over*pixel_scale/0.05)))
      pstamp_temp[0:pstamp_size[0]-1,0:pstamp_size[1]-1] = $
      pstamp_temp[0:pstamp_size[0]-1,0:pstamp_size[1]-1] + pstamp
      pstamp=pstamp_temp
      if keyword_set(plotit) then plt_image,alog10(float(pstamp)>1e-8),/fr,/col
      
      
      ; Rebin to final pixel scale
      print,"Pixellating at final pixel scale"
      pstamp=boxave(pstamp,round(over*pixel_scale/0.05))*(round(over*pixel_scale/0.05)^2)
      if keyword_set(plotit) then plt_image,alog10(float(pstamp)>1e-8),/fr,/col
      grid_centre=grid_centre/float(over*pixel_scale/0.05) 
      ;print,"grid_centre",grid_centre
      
      
      ; Convolve raw image with charge diffusion kernel
      if keyword_set(raw) then begin
        print,"Convolving with charge diffusion kernel"
        pstamp=convol(pstamp,kernel/total(kernel),/center,/edge_truncate)
      endif


      ; Renormalise flux
      if total(pstamp) gt 0 then pstamp=pstamp/total(pstamp)


      ; Measure moments using RRG
      if keyword_set(rrg) then begin
        print,"Measuring shape moments using RRG"
        rrg_pstamp=(strtrim(string(unique),2)+"_"+strtrim(string(i),2)+"_"+strtrim(string(j),2)+".fits")[0]
        rrg_cat   =(strtrim(string(unique),2)+"_"+strtrim(string(i),2)+"_"+strtrim(string(j),2)+".cat")[0]
        rrg_mom   =(strtrim(string(unique),2)+"_"+strtrim(string(i),2)+"_"+strtrim(string(j),2)+".mom")[0]
        rrg_log   =(strtrim(string(unique),2)+"_"+strtrim(string(i),2)+"_"+strtrim(string(j),2)+".mom.log")[0]
        writefits,rrg_pstamp,pstamp
        spawn,"/usr/local/optical/sextractor/bin/sex "+rrg_pstamp+$
              " -c rrg.sex -CATALOG_NAME "+rrg_cat+$
              " -PARAMETERS_NAME rrg.param",junk,junk_error
        tinytim_rrg,rrg_pstamp,rrg_cat[0],rrg_mom,min_rad=rrg_min_rad,/stellar
        restore,rrg_mom;,/verbose
        file_delete,rrg_pstamp
        file_delete,rrg_cat
        file_delete,rrg_mom
        file_delete,rrg_log
        xpos=[xpos,xdist[i,j]]
        ypos=[ypos,ydist[i,j]]
        xx=[xx,moms.xx[0]];(over*pixel_scale/0.05)^2]
        xy=[xy,moms.xy[0]];(over*pixel_scale/0.05)^2]
        yy=[yy,moms.yy[0]];(over*pixel_scale/0.05)^2]
        xxxx=[xxxx,moms.xxxx[0]];(over*pixel_scale/0.05)^4]
        xxxy=[xxxy,moms.xxxy[0]];(over*pixel_scale/0.05)^4]
        xxyy=[xxyy,moms.xxyy[0]];(over*pixel_scale/0.05)^4]
        xyyy=[xyyy,moms.xyyy[0]];(over*pixel_scale/0.05)^4]
        yyyy=[yyyy,moms.yyyy[0]];(over*pixel_scale/0.05)^4]   
        ;print,xx[1:*]
        ;print,xy[1:*]
        ;print,yy[1:*]
      endif
      
      
      ; Measure shapelet coefficients
      if keyword_set(shapelets) then begin
        decomp=shapelets_decomp(pstamp,beta,n_max,centre=grid_centre)
	      decomp.x=[x[i,j],y[i,j]]
        if i+j eq 0 then begin
	      shapecat={name:"TinyTim model stars at HST focus="+strtrim(string(focus),2),$
	                type:"shapecat",                       	 $
	                n:1L,                                  	 $
             	    maxn_max:n_max,			   	 $
             	    n_coeffs:long(decomp.n_coeffs),              $
             	    polar:0B,				   	 $
             	    x:[[x[i,j]],[y[i,j]]], 		   	 $
             	    beta:beta,	    			   	 $
             	    n_max:n_max,	    		   	 $
             	    coeffs:transpose(decomp.coeffs),       	 $
             	    coeffs_error:transpose(decomp.coeffs_error), $
             	    flag:intarr(1,2),		    	  	 $
             	    flag_interpret:["OK","OK"],	          	 $
             	    chisq:decomp.chisq[1],                	 $
		    sextractor:0B,                        	 $
		    morphology:0B,                        	 $
		    shear_estimates:0B}
	endif else shapelets_add,shapecat,decomp
      endif


      ; Find the centre of star in final pixel coords (may not be quite at
      ; old grid_centre. Could also use this modulo desired_position).
      ;grid=fltarr(pstamp_size[0],pstamp_size[1],2)
      ;max_pstamp=max(pstamp,grid_centre)
      ;grid_centre=[grid_centre mod pstamp_size[0],grid_centre/pstamp_size[0]]


      ; Locate postge stamp within larger image
      print,"Placing star into the final image"
      pstamp_size=size(pstamp,/dimensions)
      left=round(xdist[i,j]-grid_centre[0])
      if left lt 0 then begin
          pstamp=pstamp[-left:*,*]
          pstamp_size[0]=pstamp_size[0]+left
          left=0
      endif
      right=left+pstamp_size[0]-1
      if right gt simagepix[0]-1 then begin
        pstamp=pstamp[0:pstamp_size[0]+simagepix[0]-2-right,*]
        pstamp_size[0]=pstamp_size[0]+simagepix[0]-1-right
        right=simagepix[0]-1
      endif
      bottom=round(ydist[i,j]-grid_centre[1])
      if bottom lt 0 then begin
        pstamp=pstamp[*,-bottom:*]
        pstamp_size[1]=pstamp_size[1]+bottom
        bottom=0
      endif
      top=bottom+pstamp_size[1]-1
      if top gt simagepix[1]-1 then begin
        pstamp=pstamp[*,0:pstamp_size[1]+simagepix[1]-2-top]
        pstamp_size[1]=pstamp_size[1]+simagepix[1]-1-top
        top=simagepix[1]-1
      endif
      ;print,"Actual position",round(xdist[i,j]-grid_centre[0])+grid_centre[0],$
      ;                        round(ydist[i,j]-grid_centre[1])+grid_centre[1]
      
      ; Add postage stamp to image array
      if right ge 0 and left lt simagepix[0] and top ge 0 and bottom lt simagepix[1] then $
        simage[left:right,bottom:top]+= pstamp     

    endfor
  endfor


  ; Decide output file name
  if not keyword_set( fitsname ) then $
     output_file=output_dir+path_sep()+"TinyTim_f"+strtrim(string(focus),2) $
  else $
     output_file=output_dir+path_sep()+fitsname
  
  if keyword_set(decimal) then output_file=output_file+"."+decimal_string
  if keyword_set(raw) then output_file=output_file+"_flt"
  
  
  ; Write simulated RRG catalogue
  if keyword_set(rrg) then begin
    print,"Writing RRG catalogue to file "+output_file+".moms"
    print,x
    print,xxxx
    
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
    save,filename=output_file+".moms",moms
  endif

  
  ; Write simulated shapelet catalogue
  if keyword_set(shapelets) then begin
    print,"Writing shapelet catalogue to file "+output_file+".shapecat"
    shapelets_write_shapecat,shapecat,output_file,/FULL_PATH
  endif
  

  ; Write simulated image
  print,"Writing image to file "+output_file+".fits"
  file_delete,output_file+".fits",/allow_nonexistent
  if keyword_set(raw) then begin
    ;rdfits_struct,"sample_flt.fits",imagestruct
    restore,/verb,codeDir+"/tinytim_rawheaders.sav"
    fits_open,output_file+".fits",fcb,/write
    fits_write,fcb,0,header0
    fits_write,fcb,simage[0:simagepix[0]-1,0:simagepix[1]/2-1],header1,extname="SCI",extver=1;"sci,1"
    fits_write,fcb,fltarr(simagepix[0],simagepix[1]/2),header2,extname="ERR",extver=1;"err,1"
    fits_write,fcb,fltarr(simagepix[0],simagepix[1]/2),header3,extname="DQ",extver=1;"dq,1"
    fits_write,fcb,simage[0:simagepix[0]-1,simagepix[1]/2:*],header4,extname="SCI",extver=2;"sci,1"
    fits_write,fcb,fltarr(simagepix[0],simagepix[1]/2),header5,extname="ERR",extver=2;"err,1"
    fits_write,fcb,fltarr(simagepix[0],simagepix[1]/2),header6,extname="DQ",extver=2;"dq,1"
    fits_close,fcb
  endif else writefits,output_file+".fits",simage

endfor

end
