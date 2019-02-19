import tinytim as TT
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

    #0 Set up environment
    #setup_environment( cluster)
    #1. get coordinates
    if ra is None:
        ra, dec = TT.tinytim_radec_grid( cluster )

    #2. get stars
    TT.tinytim_flt( cluster, filterName, \
                        ra=ra, dec=dec, pixel_scale=pixel_scale  )

    #Change the header
    TT.tinytim_change_header( cluster )

    #. redrizzle
    TT.tinytim_drizzle( cluster,
                            pixel_scale=pixel_scale,
                            drizzle_kernel=drizzle_kernel )



def setup_environment( cluster, datadir=None, ttdir=None, 
                       cluster_type='relaxed', hst_filter='F814W' ):

    if datadir is None:
        datadir='/Users/DavidHarvey/Documents/Work/Mergers/data/'+cluster_type+'/'+cluster+'/shape/'
    if ttdir is None:
        ttdir='/Users/DavidHarvey/Documents/Work/CLASH_PSF/clusters/'+cluster

    os.system('mkdir -p '+ttdir+'/'+hst_filter)
    os.system('cp '+datadir+'/*flt*.fits '+ttdir+'/'+hst_filter)
    os.system('cp '+datadir+'/'+cluster+'*drz_sci*.fits '+ttdir+'/'+hst_filter)
    os.system('cp '+datadir+'/FocusArray.txt '+ttdir+'/'+hst_filter)
        
    
