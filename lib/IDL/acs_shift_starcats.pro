pro acs_shift_starcats

; Richard Massey - Jan 2005
;
; Change the positions in RRG catalogues of star moments, from _undist.fits
; images, to the positions they would have had after astrometric distortions.
; Writes out a new file to current directory.

; Specify where to look
indir="Jason2/"
outdir="rrg_catalogues/"

; Find all RRG catalogues
files=strmid(file_search(indir,"*_undist.fits.mom",count=n_files),strlen(indir))

; Loop over all files
for j=0,n_files-1 do begin

  ; Load old catalogue
  restore,indir+files[j];,/verbose

  ; Calculate object positions that would have been found after distortion corrn.
  dist=acs_map_xy(moms.x, moms.y)
  
  ; Initialise a new moments catalogue
  distmoms={x:dist.x,y:dist.y,x_undist:moms.x,y_undist:moms.y} 

  ; Append old variables to the catalogue
  names=tag_names(moms)
  for i=0,n_tags(moms)-1 do begin
    case strupcase(names[i]) of
      "X":
      "Y":
      else: distmoms=create_struct(distmoms,names[i],moms.(i))
    endcase
  endfor

  ; Write new catalogue
  moms=distmoms
  save,moms,filename=outdir+files[j];,/verbose

endfor

end
