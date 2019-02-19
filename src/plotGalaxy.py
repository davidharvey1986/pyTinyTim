'''
Plot the three layers of the data reduction
'''

import pyfits as fits
from matplotlib import pyplot as plt
import numpy as np

def main(galaxy='GAL-0364-52000-084'):

    psfFile='../10494_201/'+galaxy+'_F814W_post.fits'

    hduFile = fits.open(psfFile)
    print hduFile
    vmin = [ -2., 0., -4.3]
    vmax = [ 0.1, 1., -0.5]
    for i in xrange(len(hduFile[1:])):
        figure = plt.figure()
        dataLin = hduFile[i+1].data
        dataLin/=np.max(dataLin)
        if i!=1:
            print np.median(dataLin[dataLin != 0])
            dataLin[ dataLin ==0 ] = np.median(dataLin[dataLin != 0])
            data = np.log10(dataLin)
        else:
            data = dataLin

        print np.median(data), np.max(data)
        plt.imshow( data, vmin=vmin[i], vmax=vmax[i] )
        plt.axis('off')
        plt.savefig('PlotGalaxyExt_'+str(i+1)+'.pdf')
    
    
    plt.show()
