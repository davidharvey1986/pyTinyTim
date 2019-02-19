import numpy as np
import RRGtools as at

def parseGalaxyList( filename='Galaxies.lis'):
    '''
    Parse the list of galaxies copied from the 
    abstract search online
    '''
    
    
    data = np.loadtxt(filename, dtype=object,skiprows=1)
    filterList = data[:,16]

    dtype = [('Target', object), ('Filter',object), \
                 ('RA', float), ('DEC', float)]
    parsedData = np.array([], dtype=dtype)
    
    for iImage in xrange(data.shape[0]):
        
        RaHMS = data[iImage,2]+':'+data[iImage,3]+':'+data[iImage,4]
        DecHMS = data[iImage,5]+':'+data[iImage,6]+':'+data[iImage,7]
    
        raDEG, decDEG = at.hmstodd(RaHMS, DecHMS)

        BothHstFilters = data[iImage,16].split(';')
        whichNotClear = [ not 'CLEAR' in i for i in BothHstFilters]
        hstFilter=bothHstFilters[whichNotClear]
        print hstFilter

        iEntry = np.array((data[iImage,1],hstFilter, raDEG,decDEG), dtype=dtype)
        parsedData = np.append(parsedData, iEntry)


    return parsedData
    
