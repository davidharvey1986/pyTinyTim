import os as os
import glob as glob
import ipdb as pdb
def tinytim_main( cluster, filterName, ra=None, dec=None,
                      pixel_scale=0.03,
                      drizzle_kernel='square'):
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
        ra, dec = tinytim_radec_grid( cluster )

    #2. get stars
    tinytim_flt( cluster, filterName, \
                        ra=ra, dec=dec, pixel_scale=pixel_scale  )

    #Change the header
    tinytim_change_header( cluster )

    #. redrizzle
    tinytim_drizzle( cluster,
                            pixel_scale=pixel_scale,
                            drizzle_kernel=drizzle_kernel )




