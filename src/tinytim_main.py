from pyTinyTim import *
import os as os
import glob as glob
from tinytim_flt import tinytim_flt
from tinytim_change_header import tinytim_change_header
from tinytim_drizzle import tinytim_drizzle
def tinytim_main( cluster, filterName, ra=None, dec=None,
                      pixel_scale=0.03,
                      drizzle_kernel='square',\
                      dataDir='.'):
    '''

    This is the final script that determines a grid of star positions
    for the final output of PSFs.

    This will
    1. Run the grid_radec
    2. Generate flts with the stars on
    3. Redrizzle for a final psf



    '''

    #1. get coordinates
    if ra is None:
        ra, dec = tinytim_radec_grid( cluster, filterName,\
                                             data_dir=dataDir )

    #2. get stars
    tinytim_flt( cluster, filterName, dataDir=dataDir, \
                        ra=ra, dec=dec, pixel_scale=pixel_scale  )

    #Change the header
    tinytim_change_header( cluster, dataDir=dataDir )

    #. redrizzle
    tinytim_drizzle( cluster, filterName, dataDir=dataDir,
                            pixel_scale=pixel_scale,
                            drizzle_kernel=drizzle_kernel )

