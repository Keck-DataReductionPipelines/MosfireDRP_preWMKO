'''
MOSFIRE Input/Output Utility Code
Written March 2, 2011 by npk

Provides tools to read fits files and parse their headers.
'''

import pyfits as pf
import numpy as np
import unittest


import os
import pdb

import CSU

def load_edges(maskname, band, options):
    ''' Load the slit edge functions '''
    path = os.path.join(options["outdir"], maskname)
    fn = os.path.join(path, "slit-edges_{0}.npy".format(band))

    return np.load(fn)

def load_lambdacenter(fnum, maskname, options):
    ''' Load the wavelength coefficient functions '''
    path = os.path.join(options["outdir"], maskname)
    fn = os.path.join(path, "lambda_center_coeffs_{0}.npy".format(fnum))

    ld = np.load(fn)

    for i in xrange(1,47):
        assert(ld[i-1]['slitno'] == i)

    return ld

def load_lambdadata(fnum, maskname, band, options):
    ''' Load the wavelength coefficient functions '''
    path = os.path.join(options["outdir"], maskname)
    fn = os.path.join(path, "lambda_coeffs_{0}.npy".format(fnum))

    print fn
    ld = np.load(fn)
    print ld[0]['slitno']

    for i in xrange(1,47):
        assert(ld[i-1]['slitno'] == i)

    return ld

def load_lambdaslit(fnum, maskname, band, options):
    ''' Load the wavelength coefficient functions '''
    path = os.path.join(options["outdir"], maskname)
    fn = os.path.join(path, "lambda_solution_{0}.fits".format(fnum))

    print fn

    return readfits(fn)



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


def readmosfits(path, extension=None):
    '''Read a fits file written by MOSFIRE from path and return a tuple of 
    (header, data, Target List, Science Slit List (SSL), Mechanical Slit 
    List (MSL), Alignment Slit List (ASL)).
    
    Note, the extension is typically not used, only used if the detector server
    does not append slit extension.
    '''

    print path

    hdulist = pf.open(path)
    header = hdulist[0].header
    data = hdulist[0].data

    if extension is not None:
        hdulist = pf.open(extension)

    try:
        targs = hdulist[1].data
        ssl = hdulist[2].data
        msl = hdulist[3].data
        asl = hdulist[4].data
    except:
        print "Improper MOSFIRE FITS File: %s" % path

    ssl = ssl[ssl.field("Slit_Number") != ' ']
    msl = msl[msl.field("Slit_Number") != ' ']
    asl = asl[asl.field("Slit_Number") != ' ']

    bs = CSU.Barset()
    bs.set_header(header, ssl=ssl, msl=msl, asl=asl, targs=targs)

    return (header, data, bs)

def readscitbl(path):

    print path

    hdulist = pf.open(path)
    header = hdulist[0].header
    try:
        targs = hdulist[1].data
        ssl = hdulist[2].data
        msl = hdulist[3].data
        asl = hdulist[4].data
    except:
        print "Improper MOSFIRE FITS File: %s" % path

    return header, targs, ssl, msl, asl


def parse_header_for_bars(header):
    '''Parse {header} and convert to an array of CSU bar positions in mm. If 
    the positon is negative it means the barstat is not OK'''

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
        raise CSU.MismatchError("Found %i bars instead of %i" % 
                (lens(poss), CSU.numbars))
        

    return np.array(poss)


def imcombine(filelist, out, bpmask=None, reject="none", nlow=None, nhigh=None):
    '''Convenience wrapper around IRAF task imcombine
    
    reject: none, minmax, sigclip, avsigclip, pclip'''

    #TODO: REMOVE Iraf and use python instead. STSCI Python has
    # A builtin routine.
    from pyraf import iraf
    iraf.images()

    filelist = ["%s[0]" % f for f in filelist]
    pars = iraf.imcombine.getParList()
    iraf.imcombine.unlearn()

    s = ("%s," * len(filelist))[0:-1]
    s = s % tuple(filelist)
    if reject == 'minmax':
        iraf.imcombine.nlow = 1
        iraf.imcombine.nhigh = 1
        t = iraf.imcombine(s, out, Stdin=filelist, Stdout=1,
            reject=reject, masktype='badvalue', maskvalue=1, nlow=nlow,
            nhigh=nhigh)
    else:
        t = iraf.imcombine(s, out, Stdin=filelist, Stdout=1,
            reject=reject, masktype='badvalue', maskvalue=1)

    iraf.imcombine.setParList(pars)


    

class TestIOFunctions(unittest.TestCase):

    def setUp(self):
        pass

    def test_readfits(self):
        self.assertTrue(True)

if __name__ == '__main__':
    unittest.main()
