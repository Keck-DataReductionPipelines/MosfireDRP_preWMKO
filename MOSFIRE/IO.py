'''
MOSFIRE Input/Output Utility Code
Written March 2, 2011 by npk

Provides tools to read fits files and parse their headers.
'''

import numpy as np
import pyfits as pf
import CSU
from pyraf import iraf

import unittest

def readfits(path):
    '''Read a fits file from path and return a tuple of (header, data, 
    Target List, Science Slit List (SSL), Mechanical Slit List (MSL),
    Alignment Slit List (ASL)).'''
    hdulist = pf.open(path)
    header = hdulist[0].header
    data = hdulist[0].data

    return (header, data)

def readheader(path):
    '''Reads a header (only) from a fits file'''

    return pf.getheader(path)

def readfits_all(path):
    '''Read a fits file from path and return a tuple of (header, data, 
    Target List, Science Slit List (SSL), Mechanical Slit List (MSL),
    Alignment Slit List (ASL)).'''
    hdulist = pf.open(path)
    header = hdulist[0].header
    data = hdulist[0].data
    targs = hdulist[1].data
    ssl = hdulist[2].data
    msl = hdulist[3].data
    asl = hdulist[4].data

    ssl = ssl[ssl.field("Slit_Number") != ' ']
    msl = msl[msl.field("Slit_Number") != ' ']

    return (header, data, targs, ssl, msl, asl)

def parse_header_for_bars(header):
    '''Parse {header} and convert to an array of CSU bar positions in mm. If the positon is negative it means the barstat is not OK'''

    poss = []
    posfmt = "B%2.2iPOS"
    statfmt = "B%2.2iSTAT"
    for i in range(1,CSU.numbars+1):
        p = posfmt % i
        s = statfmt % i
        pos = np.float32(header[p])
        if header[s] != 'OK':
            pos *= -1
        poss.append(pos)

    if len(poss) != CSU.numbars:
        raise CSU.MismatchError("Found %i bars instead of %i" % (lens(poss), CSU.numbars))
        

    return np.array(poss)

def imcombine(filelist, out, bpmask=None, reject="none"):
    '''Convenience wrapper around IRAF task imcombine
    
    reject: none, minmax, sigclip, avsigclip, pclip'''
    iraf.images()

    filelist = ["%s[0]" % f for f in filelist]
    pars = iraf.imcombine.getParList()
    iraf.imcombine.unlearn()

    s = ("%s," * len(filelist))[0:-1]
    s = s % tuple(filelist)

    t = iraf.imcombine(s, out, Stdin=filelist, Stdout=1,
            reject=reject)


    iraf.imcombine.setParList(pars)


    

class TestIOFunctions(unittest.TestCase):

    def setUp(self):
        pass

    def test_readfits(self):
        self.assertTrue(True)

if __name__ == '__main__':
    unittest.main()
