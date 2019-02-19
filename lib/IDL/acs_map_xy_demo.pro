pro acs_map_xy_demo

;+
; NAME:
;       ACS_MAP_XY_DEMO
;
; PURPOSE:
;       Demonstrates the coordinate transformation between 
;       positions on raw ACS images, and those after correction
;       for astrometic distortions within the ACS camera.
;       Such a correction is applied by multidrizzle.
;  
; CATEGORY:
;       ACS data reduction.
;
; CALLING SEQUENCE:
;       acs_map_xy_demo
;
; INPUTS:
;       None.
;
; OPTIONAL INPUTS:
;       None.
;
; KEYWORD PARAMETERS:
;       None.
;
; OUTPUTS:
;       Plots are drawn to current output device.
;
; SUBROUTINES USED:
;       acs_map_xy.pro.
;
; EXAMPLE:
;       IDL> acs_map_xy_demo
;
; MODIFICATION HISTORY:
;       Aug 05 - Updated by RM.
;       Jan 05 - Written by Richard Massey.
;-

COMPILE_OPT idl2

; Set plotting options
ngrid=50
xran=[-200,4600]
yran=[-200,4600]
psym=7
symsize=0.5
loadct,3,/silent
colour1=150
colour2=250
offset=[199,210] ; Dither offset
pixel_scale=0.03
dxran=xran*(0.05/pixel_scale)
dyran=yran*(0.05/pixel_scale)

; Manufacture a regular-spaced grid of points across the ACS CCDs
x=(findgen(ngrid)-0.5*(float(ngrid)-1)) # replicate(1.,ngrid/2)
y=replicate(1.,ngrid) # (findgen(ngrid/2)-0.5*(float(ngrid/2)-1))
chip1x=((x/float(ngrid))+0.5)*4096   & chip2x=((x/float(ngrid))+0.5)*4096
chip1y=((y/float(ngrid/2))+0.5)*2048 & chip2y=((y/float(ngrid/2))+1.5)*2048


; Plot the initial grid  
plot,[0,0],xran=xran,yran=yran,/nodata,/xstyle,/ystyle,/isotropic,$
     title="Regularly-spaced grid of points"
oplot,[0,0],[-1,1]*1e5,psym=-3 & oplot,[-1,1]*1e5,[0,0],psym=-3
oplot,chip2x,chip2y,psym=psym,symsize=symsize,color=colour2
oplot,chip1x,chip1y,psym=psym,symsize=symsize,color=colour1
read,"Type a number to continue: ",junk


; Apply geometric distortion
distorted_chip1=acs_map_xy(chip1x,chip1y,offset=offset,pixel_scale=pixel_scale)
distorted_chip2=acs_map_xy(chip2x,chip2y,offset=offset,pixel_scale=pixel_scale)
; Another option that would work just as well
;distorted_chip1=acs_map_xy(chip1x,chip1y,offset=offset,chip=1)
;distorted_chip2=acs_map_xy(chip2x,chip2y-2048,offset=offset,chip=2)
; A final option that does the same
;distorted_chips=acs_map_xy([chip1x,chip2x],[chip1y,chip2y],offset=offset)


; Plot the distorted grid
plot,[0,0],xran=dxran,yran=dyran,/nodata,/xstyle,/ystyle,/isotropic,$
     title="Distorted grid"
oplot,[0,0],[-1,1]*1e5,psym=-3 & oplot,[-1,1]*1e5,[0,0],psym=-3
oplot,distorted_chip2.x,distorted_chip2.y,psym=psym,symsize=symsize,color=colour2
oplot,distorted_chip1.x,distorted_chip1.y,psym=psym,symsize=symsize,color=colour1
read,"Type a number to continue: ",junk


; UN-apply geometric distortion
undistorted_distorted_chip2=acs_map_xy(distorted_chip2.x,distorted_chip2.y,offset=offset,pixel_scale=pixel_scale,/inverse)
undistorted_distorted_chip1=acs_map_xy(distorted_chip1.x,distorted_chip1.y,offset=offset,pixel_scale=pixel_scale,/inverse)


; Plot the distorted grid
plot,[0,0],xran=xran,yran=yran,/nodata,/xstyle,/ystyle,/isotropic,$
     title="Undistorted distorted grid"
oplot,[0,0],[-1,1]*1e5,psym=-3 & oplot,[-1,1]*1e5,[0,0],psym=-3
oplot,undistorted_distorted_chip2.x,undistorted_distorted_chip2.y,psym=psym,symsize=symsize,color=colour2
oplot,undistorted_distorted_chip1.x,undistorted_distorted_chip1.y,psym=psym,symsize=symsize,color=colour1
read,"Type a number to continue: ",junk


; Plot the difference between this and the original positions
plot,[0,0],xran=xran,yran=yran,/nodata,/xstyle,/ystyle,/isotropic,$
     title="Difference"
oplot,[0,0],[-1,1]*1e5,psym=-3 & oplot,[-1,1]*1e5,[0,0],psym=-3
for i=0,ngrid^2/2-1 do begin
 oplot,[chip2x[i],undistorted_distorted_chip2.x[i]],[chip2y[i],undistorted_distorted_chip2.y[i]],color=colour2,psym=-3
 oplot,[chip1x[i],undistorted_distorted_chip1.x[i]],[chip1y[i],undistorted_distorted_chip1.y[i]],color=colour1,psym=-3
endfor

end
