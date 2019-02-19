import pyfits as py
import glob as glob
import numpy as np
import drizzlepac as drizzlepac
import os

def tinytim_change_header( cluster, dir_ext='', \
                               dataDir='./', tweaked=False ):

    '''
    Take a tinytim image and change the name to
    the corresponding flt name and then add the
    appropriate header

    KEYWORDS:
       TWEAKED : USE IMAGES THAT HAVE BEEN TWEAKED TO CREATE
                  PSFS THIS WAY THEY CAN BE DRIZZLED IMEDIATELY
                  ALTHOUGH THIS REQUIRES THE ENTIRE FIELD OF PSGS
                  TO BE ALREADY DRIZZLED AND HAVE EXISTING FILES
                  IN THE PRE-REQUISITE PATHS

    '''

  

    images, focus = np.loadtxt(dataDir+'/FocusArray.txt', \
                                dtype=('str'), unpack=True )


    for iImage in xrange(len(images)):
        if not tweaked:
            FLT = glob.glob(dataDir+'/'+images[iImage]+'*flt*')[0]
        else:
            FLT = dataDir+'/TinyTim/redrizzle/'+filter_string+\
                                '/'+images[iImage]+'_flt.fits'
                
        change_header(  dataDir+'/TinyTim/'+images[iImage]+'_TT_flt.fits', \
                            FLT, \
                            dataDir+'/TinyTim/'+images[iImage]+'_flt.fits')
                              

def change_header( fits_data, fits_header, out_file ):
    '''
    A function to change the header of the data file
    It take the fits header and chages the data in
    the fits_header filte

    INPTUS : FITS_HEADER : THE NAME OF THE FITS FILE
                        OF WHICH WE WANT THE HEADER FILES
             FITS_DATA : THE DATA TO BE COPIED OVER TO THE FITS_HEADER
             OTU_FILE : THE NAME IOF THE FILE TO BE WRITTEN

    OUTPUTS : WRITES THE NEW FLT TO THE OUT_FILE

    '''
    if not os.path.isfile(fits_header):
        print fits_header+' DOES NOT EXIST'
        return
    if not os.path.isfile(fits_data):
        print fits_data+' DOES NOT EXIST'
        return
    FLT_HDUs = py.open(fits_header)
    TinyTim_HDUs= py.open(fits_data)
    
    for iHDU in xrange(len(TinyTim_HDUs)):
    
        if (iHDU != 3) and (iHDU != 6) and (iHDU > 0):
            FLT_HDUs[ iHDU ].data = TinyTim_HDUs[iHDU].data.astype(np.float32)

    FLT_HDUs.writeto( out_file, clobber=True)
