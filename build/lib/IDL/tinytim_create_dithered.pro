pro tinytim_create_dithered, OVER=over, FOCUS=focus

;+      
; NAME:
;      TINYTIM_CREATE
;
; CATEGORY:
;      Mimicry of ACS (COSMOS) images.
;
; PURPOSE:
;      Makes a TinyTim model at several dither positions
;
; INPUTS:
;      Requires dither_*.txt files each containing a list of star positions.
;
; OPTIONAL INPUTS:
;      OVER  - How much to oversample (DEFAULT=50x)
;              This also determines the accuracy within a pixel to which
;              a star can be placed.
;      FOCUS - What focus values to calculate (DEFAULT=[-10,5]), in integer
;              steps.
;
; OUTPUTS:
;      Simulated starfield images written to disc.
;
; SUBROUTINES REQUIRED:
;      tinytim_create
;
; MODIFICATION HISTORY:
;      Apr 05 - Written by Richard Massey.
;-

; Count files containg star positions in each dithered pointing
n_dithers=n_elements(file_search("dither_*.txt"))

; Set default options
if n_elements(over)  eq 0 then over=50
if n_elements(focus) eq 0 then focus=[-10,5]

; Loop over each dither position
for i=1,n_dithers do begin

  ; Read in coordinates for one dither position dither position
  filename="dither_"+strtrim(string(i),2)+".txt"
  readcol,filename,x,y,chip
  chip1=where(chip eq 1)
  y[chip1]=y[chip1]+2048 

  ; Make the image at that dither position
  tinytim_create, X=x, Y=y, /RAW, /EXACT_POSITION, $
                  OVER=over, FOCUS=focus, OUTPUT_DIR=filename+"_output"

endfor


end
