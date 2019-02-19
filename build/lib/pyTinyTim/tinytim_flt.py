import drizzlepac 
import numpy as np
import ipdb as pdb
import glob as glob
import tinytim_run as tinytim_run
import tinytim_change_header as tt_header
import nearest_neighbour as NN
import os as os
from RRGtools import *
import re
def tinytim_flt( cluster, filter_string, ra=None, dec=None, tweaked=False,
                     pixel_scale=0.03, dataDir='./' ):
    '''

    A function that creates fake flt files that are popualted
    with stars at the positions of the objects given in a catalogue

    Each flt will have the same header as the true flt from the
    data

    Will output N FLT images of PSFS with focus positions as
    determined by RRG_psf_cor

    So therefore requries some pre working

    1. Drizzled image of the cluster and reduced FLT files
       in the folder

           /Users/DavidHarvey/Documents/Work/CLASH_PSF/CLUSTER/FILTER

    2. A catalogue of object positions at which to calcualte the PSF
       named

          CLUSTER_clash.cat

    3. A file called FocusArray.txt that is produced during the shape
       measurement process of rrg. This file contains N FLT lines with
       image name (minus any extension such as _flt.fits) and then the
       estiamted focus positions of HST ACS at the time of obseravtion


    4. Once I have these in the correct folders, it should be good to go.

    REQUIRES : TINYTIM, TINYTIM_CREATE.PRO, DRIZZLEPAC

    KEYWORDS :
        RA, DEC : COODAITNES GIVEN BY USER TO USE, IF NONE
                  IT WILL LOOK FOR A CATALOGUE
                  
        TWEAKED : USE IMAGES THAT HAVE BEEN TWEAKED TO CREATE
                  PSFS THIS WAY THEY CAN BE DRIZZLED IMEDIATELY
                  ALTHOUGH THIS REQUIRES THE ENTIRE FIELD OF PSGS
                  TO BE ALREADY DRIZZLED AND HAVE EXISTING FILES
                  IN THE PRE-REQUISITE PATHS

    UPDATES: Changed the datadir paths to current working directory.
             And removed loop through filters

    '''

    images, focus = np.loadtxt(dataDir+'FocusArray.txt', \
                                dtype=('str'), unpack=True )

    if ra is None or dec is None:
        ra, dec = np.loadtxt( dataDir+'/'+cluster+'_clash.cat', \
                                  usecols=(1,2), unpack=True )

    outputDir = dataDir+'/TinyTim/'
        
    nImages = len(focus)
        
    for iImage in xrange(nImages):
        if nImages == 1:
            image_use = images
        else:
            image_use = images[iImage]
                
        fits =  glob.glob(dataDir+'/'+image_use+'*flt*')
        if len(fits) != 1:
            raise ValueError('Unexpceted number of filenames: %s', \
                                 len(fits))

        else:
            fits = fits[0]
                
        
        x, y = deg2pix_flt( fits, ra, dec, postage_stamp=50.)
        
            
            
        if len(x) == 0 or len(y) == 0:
            print 'No galaxies in field'
            continue

        wavelength=re.findall(r'\d+',filter_string)[0]
        
        tinytim_run.run( x, y, focus_range=focus[iImage], \
                             pixel_scale=pixel_scale, \
                             output_dir="'"+outputDir+"'", \
                            filter="'"+filter_string+"'",\
                            raw=1, fitsname="'"+image_use+"_TT'",\
                            exact_position=1, wavelength=wavelength)



def postage_stamp( cluster, postage_size=80, dataDir='./' ):
    '''
    Create postage stamps of FLTs

    To do this the plan is to take the positions of all
    the ra and dec, turn them into x and y
    
    Then find where the nearest neighbours are
    use only those stars that are far enough away from others that
    they can be used.

    With those objects that have close neighbours, I will choose
    one of the neighbours and append that to the fits image

    Then i will loop back around and do the other pair in that image

    And if it is a triplet then move back around and creat that one.

    So I will out put many fits images, each one where the closest neighbour
    is greater than the required distance
    '''
    
    ra_hms = []
    dec_hms = []

    images, focus = np.loadtxt(dataDir+'/FocusArray.txt', \
                                dtype=('str'), unpack=True )

    ra, dec = np.loadtxt( dataDir+'/'+cluster+'_clash.cat', \
                                usecols=(1,2), unpack=True )

    nImages = len(focus)
        
    for iImage in xrange(len(images)):
            
        outputDir = dataDir+'/TinyTim/'+images[iImage]
        if not os.path.isfile( outputDir ):
            os.system("mkdir -p "+outputDir)

        if nImages == 1:
            image_use = images
        else:
            image_use = images[iImage]

        fits =  dataDir+'/'+image_use+'_flt.fits'
            
        x_todo, y_todo = deg2pix_flt( fits, ra, dec)
        n_iterations = 0
            
        while len(x_todo) > 0:
                
            n_iterations +=1

                
            x_drizzle, y_drizzle, x_run, y_run = \
                    pick_isolated(  x_todo, y_todo, \
                                    postage_size=postage_size)
            x_todo_run = []
            y_todo_run = []

            for iTodo in xrange(len(x_run)):
                #Now loop through each one and make
                #the point is not close to another
                iDist = np.sqrt((x_run[iTodo]-x_drizzle)**2 + \
                                    (y_run[iTodo]-y_drizzle)**2)

                if np.all(iDist > postage_size ):
                    x_drizzle = np.append(x_drizzle, x_run[iTodo]) 
                    y_drizzle = np.append(y_drizzle, y_run[iTodo])
                else:
                    x_todo_run = np.append( x_todo_run, x_run[iTodo])
                    y_todo_run = np.append( y_todo_run, y_run[iTodo])

                ra_run, dec_run = pix2deg_flt(fits, x_drizzle, y_drizzle)

                #Catalogue which galaxies have gone in.

                np.savetxt( outputDir+"/"+image_use+"_"+str(n_iterations)+".cat",\
                             np.transpose([ra_run, dec_run]), fmt='%10.5f',delimiter=',')
    
        x_todo = np.array(x_todo_run)
        y_todo = np.array(y_todo_run)

def pick_isolated( x, y, postage_size=None ):
    
    distance, closest_index = NN.nearest_neighbour(x, y)

    isolated_x = x[ distance > postage_size ]
    isolated_y = y[ distance > postage_size ]
    iso_index = np.arange(len(distance))[ distance > postage_size ]

    x_todo = np.delete(x, iso_index)
    y_todo = np.delete(y, iso_index)

    return isolated_x, isolated_y, x_todo, y_todo
    

def remove_edges( x, y, fits):
    '''
    Remove the postage stamps that are too close to the edges of
    the field

    '''
