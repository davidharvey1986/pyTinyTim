from pyHST import drizzle as drizzle
import tinytim_change_header as tt_flt
import ipdb as pdb
import glob as glob
import numpy as np
from subprocess import call
import os as os
import numpy as np
import pyfits as fits

def tinytim_drizzle( cluster, combine_type='iminmed', output=None,
                         drizzle_kernel='square',
                         pixel_scale=0.03, dataDir='./'):

        
        
    #Make sure it is a filter!
    
    TT_dir = 'TinyTim/'

    os.system( 'mkdir -p '+TT_dir+'/redrizzle' )
        
 
    #TWeak the fake flts back using those in the
    #actual drizzle process
    
    #I dont need to tweak since i use the header information
    #from the flt files.
    
    os.environ['jref'] = TT_dir
        
    input_str = TT_dir+'/*q_flt.fits'

    if output is None:
        output = TT_dir+'/redrizzle/'+cluster+'_TT'
    thresh=1.0
    search_rad=1.0
    
    if fits.__version__ != '3.1.6':
        raise ImportError('Not the correct version of pyfits, needs 3.1.6')
    if np.__version__ != '1.11.0':
        raise ImportError('Not the correct version of numpy, needs 1.11.0')
        
    drizzle.astrodrizzle.AstroDrizzle( input_str, \
                                            output=output, \
                                            final_scale=pixel_scale, \
                                            final_pixfrac=0.8, \
                                            final_wcs=True, \
                                            combine_type='iminmed', \
                                            skysub=True,
                                            driz_separate=False,
                                            static=False, \
                                            median=False,
                                            blot=False,
                                            driz_cr=False,
                                            final_kernel=drizzle_kernel)
        
    os.system("rm -fr "+TT_dir+"/*mask*")
