import glob as glob
import pyfits as py
import astro_tools as at
import numpy as np
import ipdb as pdb
def tinytim_radec_grid( cluster, n_stars=25, datadir=None):
    '''
    The purpose of this script is to create a grid of ra and dec positions in a grid
    on the final drizzled image.

    These will be passed into the tinytim_flt.py script.

         cluster : string of the clsuter in quertion

    ketwords :
         n_stars : how many stars per side

    '''

    if datadir is None:
        datadir='/Users/DavidHarvey/Documents/Work/CLASH_PSF/clusters/'+cluster

    #Get the hst filter lists
    hst_filter = glob.glob( datadir+'/F*' )[-1].split('/')[-1]

    hst_image = cluster+'_'+hst_filter+'_drz_sci.fits'
    
    
    drz_image = py.open( datadir+'/'+hst_filter+'/'+hst_image)


    x_vec = np.linspace( 0, drz_image[0].header['NAXIS1'], n_stars )
    y_vec = np.linspace( 0, drz_image[0].header['NAXIS2'], n_stars )

    x_grid, y_grid = np.meshgrid( x_vec, y_vec )

    
    ra, dec = at.pix2deg( datadir+'/'+hst_filter+'/'+hst_image, x_grid.reshape( n_stars**2 ), y_grid.reshape( n_stars**2) )



    
    return ra, dec
    

    
        
