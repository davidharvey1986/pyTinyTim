import numpy as np
import pyfits as py
import RRGtools as at
import pyRRG as rrg
import ipdb as pdb
import glob as glob
import os as os
'''
Generate a postage_stamp for each
object, 501x501 pixels across,
for the wht_image, sci_image and psf_image

'''

def postage_stamp( image_name, x, y, width=500, ext=0 ):
    '''
    return the postage stamp for the given
    image from x and y

    '''
        
    if width/2. != width/2:
        width_half = (width - 1)/2
        width_upper = width_half+1
    else:
        width_half = width/2
        width_upper = width_half
        
    image = py.open( image_name )[ext].data
    
    return image[ y-width_half:y+width_upper,
                   x-width_half:x+width_upper]
    
def tweak_pos(  image_name, ra, dec, width=500, ext=0, tweakMax=100 ): 

    '''
    return the postage stamp for the given
    image from x and y

    '''
    if width/2. != width/2:
        width_half = (width - 1)/2
        width_upper = width_half+1
    else:
        width_half = width/2
        width_upper = width_half
    x, y = at.deg2pix( image_name, [ra], [dec] )
    print('CONVERTED %0.5f %0.5f -> %0.5f %0.5f' %\
          (ra, dec, x, y))

    image = py.open( image_name )[ext].data

    
    y_int = np.int(np.floor(y))
    x_int = np.int(np.floor(x))
   
    tmp_image =  image[ y_int-width_upper:y_int+width_half,
                   x_int-width_upper:x_int+width_half]

    xgrid, ygrid = np.meshgrid( np.arange(width), np.arange(width) )

    #I need a prior on how far the tweaking is
    rGrid = np.sqrt((xgrid-width_upper)**2 + (ygrid-width_upper)**2)
    
    maxInPrior = np.max(tmp_image[rGrid<tweakMax])

    xnew = np.int(xgrid[maxInPrior== tmp_image] +\
                   ( x_int - width_upper ))
    ynew = np.int(ygrid[maxInPrior == tmp_image] + \
                  ( y_int - width_upper ))
   
    return xnew, ynew

def all_stamps( psf=True, hstFilter='F814W', width=301, outDir='output' ):
    '''
    Read the data file and get the name and the
    ra and dec of the object

    '''
    
    Lenses = glob.glob('/Users/DavidHarvey/Documents/Work/pyHST/SLACS/GAL*/GAL*'+hstFilter+'_drz_sci.fits')

    for iLens in Lenses:

        LensDir =  '/'.join(iLens.split('/')[0:-1])
        
        LensName = iLens.split('/')[-1].split('_')[0]
        postName =  '../'+outDir+'/'+LensName+'_'+hstFilter+'_post.fits'
        print postName
        if os.path.isfile(postName):
            continue                      
        ra_deg = py.open(iLens)[0].header['RA_APER']
        dec_deg = py.open(iLens)[0].header['DEC_APER']
        
        wht_name = LensDir+'/'+LensName+'_'+hstFilter+'_drz_wht.fits'
        psf_name = LensDir+'/TinyTim/redrizzle/'+LensName+'_'+hstFilter+'_TT_drz_sci.fits'

        if not os.path.isfile(wht_name):
            raise ValueError('Cant find %s weight file' % wht_name)
        
        xnew, ynew = tweak_pos( iLens, ra_deg, dec_deg, width=width )
        
        img_post =  postage_stamp( iLens, xnew, ynew, width=width )
        wht_post =  postage_stamp( wht_name, xnew, ynew, width=width )
        
        
        xpsf, ypsf = tweak_pos( psf_name, ra_deg, dec_deg, width=width )
        psf_post =  postage_stamp( psf_name, xpsf, ypsf, width=width )
        
        hdus = py.HDUList([py.PrimaryHDU(),
                            py.hdu.ImageHDU(img_post),
                            py.hdu.ImageHDU(wht_post),
                            py.hdu.ImageHDU(psf_post)])
                            
        
    
        skymed, skysig, skyskew = rrg.mmm( py.open(iLens)[0].data )

        print("SKY MED : %0.7f \nSKY SIG : %0.7f\n" % \
            (skymed, skysig ))

        hdus[0].header['SKYMED'] = (skymed, 'MEDIAN NOISE BKG NOISE')
        hdus[0].header['SKYSIG'] = (skysig, 'STANDARD DEVIATION IN NOISE' )
        hdus[0].header['SKYSKEW'] = (skyskew, 'SKEWNESS IN NOISE DISTRIBUTION')
    
        hdus.writeto( postName, clobber=True )
    
    

    
    
    
    
    
