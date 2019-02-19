import numpy as np
import RRGtools as at

def parseGalaxyList( filename='Galaxies.lis', hstFilter='F814W'):
    '''
    Parse the list of galaxies copied from the 
    abstract search online
    '''
    
    
    data = np.loadtxt(filename, dtype=object,skiprows=1)
    filterList = data[:,16]

    dtype = [('Target', object), ('Filter',object), \
                 ('RA', float), ('DEC', float)]
    parsedData = np.array([], dtype=dtype)
    
    for iImage, iFilter in enumerate(filterList):

        if 'F814W' in iFilter:
            RaHMS = data[iImage,2]+':'+data[iImage,3]+':'+data[iImage,4]
            DecHMS = data[iImage,5]+':'+data[iImage,6]+':'+data[iImage,7]

            raDEG, decDEG = at.hmstodd(RaHMS, DecHMS)
            iEntry = np.array((data[iImage,1],hstFilter, raDEG,decDEG), dtype=dtype)
            parsedData = np.append(parsedData, iEntry)


    return parsedData
    
