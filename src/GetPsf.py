'''
Run RRG on the image and get the PSF image of the
exposure in question at the RA and dec of the target
'''

import pyRRG as RRG
import os as os
import pyfits as fits
from parseGalaxyList import *
import sys
from pyTinyTim import *
def main( galaxyList='Galaxies.lis'):
    GalaxyList = parseGalaxyList.parseGalaxyList( filename=galaxyList)

        
    for iGalaxy in GalaxyList:

        getPSF(iGalaxy)


def getPSF( galaxy ):
    TargetDir=galaxy['Target']+'/'
    TargetImage = galaxy['Target']+'_'+galaxy['Filter']+'_drz_sci.fits'
    print TargetImage
    if not os.path.isfile(TargetDir+'/FocusArray.txt'):
        measurePSF(TargetDir,TargetImage)

    psfFile = TargetDir+'/TinyTim/redrizzle/'+\
      galaxy['Target']+'_'+galaxy['Filter']+'_TT_drz_sci.fits'
    
    if not os.path.isfile(psfFile):
        print galaxy['Target'], galaxy['RA'], galaxy['DEC']
        tinytim_main( galaxy['Target'], \
                            galaxy['Filter'], \
                        ra=[galaxy['RA']], \
                        dec=[galaxy['DEC']],\
                        dataDir=TargetDir)
                      
def measurePSF( TargetDir,TargetImage ):
    codeDir = os.getcwd()
    os.chdir(TargetDir)
    try:
        RRG.main(TargetImage)
        os.chdir(codeDir)
    except:
        os.chdir(codeDir)
        raise
 
                    
