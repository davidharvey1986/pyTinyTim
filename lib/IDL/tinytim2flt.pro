pro tinytim2flt, image, tinytim_name, $
                 data_dir=data_dir, $
                 output_dir=output_dir

; A program to take the focys position in dataDir and copy it such that
; it resembles the flt of the guiven imageName

;image : the filename of the flt to copy to
;tinytim : the tinytim model to copy over
;output_dir : the output directory
  
  if not keyword_set(data_dir) then output_dir='.'
  if not keyword_set(output_dir) then output_dir='.'

  fits_info, data_dir+'/'+image, $
             n_ext=n_ext, extname=ext_name, $
             textout=2
  
  fits_open, output_dir+'/'+image, fcb, /write

  for iExt = 0, 1 do begin
     if iExt ne 3  then begin
        fits_read,  data_dir+'/'+image, $
                    data, header, $
                    exten_no=iExt, $
                    EXTNAME=ext_name[iExt]
               

        fits_write, fcb, data, header,$
                    EXTNAME=ext_name[iExt]
     endif
     
        
     
  endfor

  fits_close,fcb


end


  
