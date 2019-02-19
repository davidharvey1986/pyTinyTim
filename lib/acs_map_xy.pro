function acs_map_xy, x, y, INVERSE=inverse,        $
                           ORDER=order,            $
                           NOSTORE=nostore,        $
                           OFFSET=offset,          $
                           CHIP=chip,              $
                           PIXEL_SCALE=pixel_scale

;+
; NAME:
;       ACS_MAP_XY
;
; PURPOSE:
;       Provides the coordinate tranformation to map between 
;       positions on raw ACS images, and those after correction
;       for astrometic distortions within the ACS camera.
;       Such a correction is applied by multidrizzle.
;
;       Note that only the DRIZZLE polynomial is calculated; the
;       distortion map adjusts positions by a further ~0.1 pixels.
;       This needs to be applied separately, but if the inverse
;       transformation is needed, at least an addition is easy to
;       invert!
;  
; CATEGORY:
;       ACS data reduction.
;
; CALLING SEQUENCE:
;       result=acs_map_xy(x_coords,y_coords)
;
; INPUTS:
;       X           - Object x coordinate, or an array of them.
;       Y           - Object y coordinate, or an array of them.
;
; OPTIONAL INPUTS:
;       PIXEL_SCALE - The output pixel scale as used by DRIZZLE.
;       OFFSET      - The [x,y] position of the bottom-left pixel after
;                     correction for geometric distortions. If an image
;                     consists of multiple, dithered exposures, this
;                     should vary for each exposure, reflecting the
;                     overall dither pattern.
;       CHIP        - Which CCD the objects are on. If omitted, it is
;                     assumed that the objects are all on a combined 
;                     4096x4096 image, with no gap between the CCDs.
;                     Chip #1 is taken to be the bottom one.
;       ORDER       - Maximum order of polynomial used for the 
;                     transformation. Default is set by the number of read
;                     ceofficients (which gives a quartic).
;
; KEYWORD PARAMETERS:
;       INVERSE     - Apply the inverse transformation, from the frame
;                     after correction to the raw CCD frame.
;
; OUTPUTS:
;       Returns a structure containing transformed x' and y' coords.
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
;       A full demonstration is provided in acs_map_xy_demo.pro.
;
; MODIFICATION HISTORY:
;       Nov 05 - Memory use with large arrays streamlined by RM.
;       Nov 05 - Chips distinguished with different pixels scales by RM.
;       Aug 05 - PIXEL_SCALE keyword added by RM.
;       Jan 05 - Inverse distortion coefficients supplied by Herve Aussel.
;       Jan 05 - Written by Richard Massey.
;-

COMPILE_OPT idl2

; Desired position of bottom-left pixel after mapping, to enable dithering
; (STScI standard notation seems to use the centre of the chip instead).
if not keyword_set(pixel_scale) then pixel_scale=0.05
if not keyword_set(offset) then begin
  if pixel_scale eq 0.05 then offset=[199,210] $
    else if pixel_scale eq 0.03 then offset=[247,224] $
    else begin
      offset=[0.,0.]
      message,"Default offset not known for this pixel scale! Assuming zero.",/info
    endelse
endif
n_positions=n_elements(x)
if n_elements(y) ne n_positions then message,"x and y arrays are of different sizes!"


; Set polynomial coefficients for transformation (or read in using read_coeffs.pro).
if not keyword_set(inverse) then begin
  ; Centre of the two chips.
  origin1=[2048.,1024.]
  origin2=[2048.,3072.]
  ; Coefficients obtained from *_flt_coeffs2.dat, output by DRIZZLE
  ;chip1xcoeffs=[2048-56.69626992,0.99648395,0.039235942,8.443148e-06,-4.9350914e-06,1.8714837e-06,-4.7000765e-10,-6.3905282e-11,-5.2871972e-10,-5.8953743e-12,2.4971845e-14,4.9862677e-15,3.637933e-14,-1.8702128e-14,1.405061e-14]
  ;chip1ycoeffs=[2048-1030.25854442,0.029775126,1.0053584,-1.4995461e-06,6.0549338e-06,-7.2055625e-06,7.1958641e-11,-5.1743686e-10,-8.3696034e-11,-4.2292829e-10,-1.7454337e-14,-1.9344815e-15,-2.8682819e-14,1.3391103e-14,1.0777943e-15]
  ;chip2xcoeffs=[2048+28.69390174,0.98437314,0.045339576,8.2503964e-06,-7.1088858e-06,1.7728154e-06,-4.5976284e-10,-1.2408655e-10,-5.29678e-10,4.2119219e-11,1.7393012e-14,3.2388574e-14,9.2331066e-14,-2.0957229e-14,1.3739822e-13]
  ;chip2ycoeffs=[2048+1047.39000180,0.042595268,0.97152006,-2.51594e-06,5.9160766e-06,-9.3900933e-06,6.7942709e-11,-4.3467827e-10,-5.3234262e-11,-3.9818745e-10,-1.6672274e-14,-3.5575914e-14,-1.0733564e-13,2.7255259e-14,-1.680405e-13] 
  ; Coefficients obtained from oas1611bj_idc.fits, input to DRIZZLE (almost identical, but with an overall scale factor needed to convert from pixels to arcseconds).
  chip1xcoeffs=[99.565188,4.9825169E-02,1.9619130E-03, 4.2225790E-07,-2.4670911E-07, 9.3607483E-08,-2.3507861E-11, -3.2000680E-12, -2.6439229E-11, -2.9675790E-13, 1.2487490E-15, 2.4934469E-16, 1.8191949E-15,-9.3522379E-16,  7.0261870E-16]/5D-2
  chip1ycoeffs=[50.887074,1.4886410E-03,5.0269939E-02,-7.4966152E-08, 3.0282621E-07,-3.6023110E-07, 3.6018710E-12, -2.5871029E-11, -4.1843920E-12, -2.1149299E-11,-8.7282642E-16,-9.6736218E-17,-1.4343210E-15, 6.6963921E-16,  5.3896481E-17]/5D-2
  chip2xcoeffs=[103.83470,4.9219720E-02,2.2672340E-03, 4.1262069E-07,-3.5540299E-07, 8.8666752E-08,-2.2995350E-11, -6.2186931E-12, -2.6492589E-11,  2.0806501E-12, 8.6975979E-16, 1.6196320E-15, 4.6171330E-15,-1.0479930E-15,  6.8707738E-15]/5D-2
  chip2ycoeffs=[154.76950,2.1297960E-03,4.8578139E-02,-1.2579410E-07, 2.9587099E-07,-4.6947380E-07, 3.4025090E-12, -2.1720251E-11, -2.6551810E-12, -1.9880131E-11,-8.3371838E-16,-1.7790191E-15,-5.3674559E-15, 1.3629340E-15, -8.4030801E-15]/5D-2
endif else begin
  ; Centre of the two chips after forward transformation.
  origin1=[2047.9090,1089.8271]+offset
  origin2=[2133.2993,3167.4756]+offset
  ; Reverse transformation coefficients, calculated using Herve Aussel's fit_polynomial_distortion.pro.
  chip2xcoeffs=[2048+17.690567,1.0098704,-0.041546981,-9.0331744e-06,7.5455644e-06,-3.5491556e-06,6.5557962e-10,-1.8975432e-11,5.0821254e-10,6.3986712e-10,-3.5403050e-14,1.2697593e-14,-1.4901954e-13,1.0700322e-13,-1.9272843e-13]
  chip2ycoeffs=[4096-1067.8955,-0.036814084,1.0101043,2.8356072e-06,-7.4504696e-06,1.0394685e-05,-1.9416512e-10,3.4142281e-10,-7.1805489e-11,-4.5094846e-10,2.4085682e-14,-2.2323542e-15,1.8504679e-13,-1.1194485e-13,2.6685584e-13]
  chip1xcoeffs=[2048+14.317421,1.0106205,-0.043983312,-8.9777910e-06,7.4831681e-06,-2.6797870e-06,6.3869379e-10,-3.7342433e-10,7.9598247e-10,-1.5421049e-10,-3.9667949e-14,2.1662802e-14,-6.4200504e-14,4.3270883e-14,-1.8717430e-14]
  chip1ycoeffs=[+1032.3554,-0.036663688,1.0124915,2.6735703e-06,-7.1984241e-06,9.1768566e-06,-1.6076053e-10,7.8555927e-10,-3.0590836e-10,6.1370412e-10,2.1460207e-14,-2.1748604e-14,5.9100563e-14,-3.2919786e-14,1.9975680e-14]
endelse
if not keyword_set(order) then order=!values.f_infinity
order=order<round((sqrt(9-8*(1-n_elements(chip1xcoeffs)))-3)/2.) ; i.e. quartic!


; Decide which chip the objects are on.
if not keyword_set(chip) then begin
  ; Objects are on some combined image (assumed to be 4096x4096, with no gap).
  if keyword_set(inverse) then begin
    ; Work out which is on which CCD, by first getting the distorted CCD boundary
    boundary=acs_map_xy([0,4095,0,4095],[2047,2047,2048,2048],order=order,pixel_scale=pixel_scale,offset=offset)
    boundary_x1=mean([boundary.x[0],boundary.x[2]])
    boundary_x2=mean([boundary.x[1],boundary.x[3]])
    boundary_y1=mean([boundary.y[0],boundary.y[2]])
    boundary_y2=mean([boundary.y[1],boundary.y[3]])
    boundary_m=(boundary_y2-boundary_y1)/(boundary_x2-boundary_x1)
    boundary_c=boundary_y1-boundary_m*boundary_x1
    chip1=where(y lt boundary_m*x+boundary_c,n_chip1)
    chip2=where(y ge boundary_m*x+boundary_c,n_chip2)
    ;chip1=where(y lt (2098+offset[1])*0.05/pixel_scale,n_chip1)
    ;chip2=where(y ge (2098+offset[1])*0.05/pixel_scale,n_chip2)
  endif else begin
    ; Boundary is a horizontal line between pixels 2047 and 2048
    chip1=where(y lt 2048.,n_chip1)
    chip2=where(y ge 2048.,n_chip2)
  endelse
  chipid=byte(x-x)+1B
  if n_chip2 gt 0 then chipid[chip2]=2B
endif else begin
  ; All objects are on one of the chips, specified by the user.
  if keyword_set(inverse) then message,"Not yet implemented this combination of options!"
  if n_elements(chip) gt 1 then message,"Not yet implemented this option!"
  case chip of
    1: begin
         n_chip1=n_elements(x)
         chip1=lindgen(n_chip1)
         n_chip2=0
       end
    2: begin
         n_chip2=n_elements(x)
         chip2=lindgen(n_chip2)
         n_chip1=0
         origin2[1]=origin2[1]-2048
       end
    else: message,"ACS only has two CCDs!"
  endcase
  chipid=byte(x-x+chip)
endelse


; Create some empty arrays of the right size, in which to put the final answer.
xprime=float(x-x)
yprime=float(y-y)


; Perform transformation via the polynomial model.
; THIS IS BLOODY SLOW IF THERE ARE A LOT OF POSITIONS!
if n_chip1 gt 0 then begin
  k=0
  for i=0,order do begin
    for j=0,i do begin
      xprime[chip1]+=chip1xcoeffs[k]*(x[chip1]-origin1[0])^(i-j)*(y[chip1]-origin1[1])^(j)
      yprime[chip1]+=chip1ycoeffs[k]*(x[chip1]-origin1[0])^(i-j)*(y[chip1]-origin1[1])^(j)
      k+=1
    endfor
  endfor
endif
delvarx,chip1
if n_chip2 gt 0 then begin
  k=0
  for i=0,order do begin
    for j=0,i do begin
      xprime[chip2]+=chip2xcoeffs[k]*(x[chip2]-origin2[0])^(i-j)*(y[chip2]-origin2[1])^(j)
      yprime[chip2]+=chip2ycoeffs[k]*(x[chip2]-origin2[0])^(i-j)*(y[chip2]-origin2[1])^(j)
      k+=1
    endfor
  endfor
endif
delvarx,chip2


; Free up some memory
if keyword_set(nostore) then begin
  delvarx,x
  delvarx,y
endif


; Shift to dithered position (polynomial fit has a strange convention to put bottom-
; left pixel at (-56.6,72.1) by default. Perhaps some other place stays invariant).
if not keyword_set(inverse) then begin
  xprime=(temporary(xprime)+56.605274)*0.05/pixel_scale+offset[0]
  yprime=(temporary(yprime)+72.085616)*0.05/pixel_scale+offset[1]
endif else begin
  xprime=temporary(xprime)/0.05*pixel_scale;-56.605274
  yprime=temporary(yprime)/0.05*pixel_scale;-72.085616
endelse


; Tell the world what we have found!
new_coeffs={chip:temporary(chipid)}
new_coeffs=create_struct(temporary(new_coeffs),"x",temporary(xprime))
new_coeffs=create_struct(temporary(new_coeffs),"y",temporary(yprime))
return,new_coeffs

end
