pro tinytim_flt, cluster

;PURPOSE : TO CREATE N NUMBER OF FLT FILES WITH
;(FOR NOW) A GRID OF STARS REPRESENTING THE PSG
;FOR EACH STAR
  pixel_scale=0.03/0.049
  nstars=30


  rootDir = '/Users/DavidHarvey/Documents/Work/CLASH_PSF/'
  dataDir = rootDir+'/clusters/'+cluster
  spawn,'ls '+dataDir,filters
  nstars=10
  print, filters

  for iFilter=0,n_elements(filters)-1 do begin
     
     readcol,dataDir+'/'+filters[iFilter]+'/'+cluster+'_clash.cat',$
        ID, ra_list, dec_list, format='D,D,D',/silent

     ra_list = ra_list[0:200]
     dec_list = dec_list[0:200]
     
     outputDir = dataDir+'/'+filters[iFilter]+'/TinyTim'
 
     ;; openw,luna,outputDir+'/catfile',/get_lun

     ;First read FocusArray to get the image names
     ;and the focus positions

     readcol,dataDir+'/'+filters[iFilter]+'/FocusArray.txt',$
             imageName, focus, format='A,I'
     
     for iImage=0, n_elements( imageName)-1 do begin

    
        sky2xy, dataDir+'/'+filters[iFilter]+'/'+$
               imageName[iImage]+'_flt.fits', $
               ra_list, dec_list, $
                x_list, y_list, /hst_flt
        
        ;; printf, luna, imageName[iImage]+'_flt.fits  '+outputDir+'/'+$
        ;;         imageName[iImage]+'_sci1.cat '+$
        ;;         outputDir+'/'+imageName[iImage]+'_sci2.cat'
        
        tinytim_create, nstars, $
                        focus_range=focus[iImage], $
                        x = x_list, y=y_list, $
                        pixel_scale=0.03, $
                        output_dir=outputDir, $
                        filter=filters[iFilter],$
                        /raw, fitsname=imageName[iImage]+'_TT',$
                        /exact_position

        ;; openw,lun,outputDir+'/'+imageName[iImage]+'_sci1.cat',/get_lun
        ;; openw,lun2,outputDir+'/'+imageName[iImage]+'_sci2.cat',/get_lun
        
        ;; for iStar=0, n_elements(x_list)-1 do begin
        ;;    if  y_list[iStar] gt 0 and y_list[iStar] lt 2048 $
        ;;       and x_list[iStar] gt 0 and x_list[iStar] lt 4096 then $
        ;;       printf, lun, x_list[iStar], y_list[iStar], 20, strtrim(iStar,2), $
        ;;               format='(F10.5,1X,F10.5,1X,F10.5,1X,I5)' $
        ;;    else if y_list[iStar] gt 2048 and y_list[iStar] lt 4096 $
        ;;       and x_list[iStar] gt 0 and x_list[iStar] lt 4096 then $
        ;;       printf, lun2, x_list[iStar], y_list[iStar]-2048, 20, strtrim(iStar,2), $
        ;;               format='(F10.5,1X,F10.5,1X,F10.5,1X,I5)' 
        ;; endfor
        ;; close, lun
        ;; close, lun2
        ;; free_lun, lun
        ;; free_lun, lun2
     endfor
     ;; openw, lun, outputDir+'/reference.cat' 
     ;; for i=0, n_elements(ra_list)-1 do $
     ;;    printf, lun, ra_list[i], dec_list[i], 20, strtrim(i,2), $
     ;;            format='(F10.5,1X,F10.5,1X,F10.5,1X,I3)'
     ;; close, lun
     ;; close,luna
     ;; free_lun,lun
     ;; free_lun,luna
  endfor
end
  
  
  
