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

import MOSFIRE
import CSU
import Options

theBPM = None # the Bad pixel mask

def badpixelmask():
    global theBPM

    path = Options.path_bpm

    if theBPM is None:
        hdulist = pf.open(path)
        header = hdulist[0].header
        theBPM = hdulist[0].data

    return theBPM

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

    return ld

def load_lambdadata(fnum, maskname, band, options):
    ''' Load the wavelength coefficient functions '''
    path = os.path.join(options["outdir"], maskname)
    fn = os.path.join(path, "lambda_coeffs_{0}.npy".format(fnum))

    ld = np.load(fn)

    return ld

def load_lambdaoutwards(fnum, maskname, band, options):
    ''' Load the wavelength coefficient functions '''
    path = os.path.join(options["outdir"], maskname)
    fn = os.path.join(path, "lambda_outwards_coeffs_{0}.npy".format(fnum))

    ld = np.load(fn)

    return ld

def load_lambdamodel(fnum, maskname, band, options):
    ''' Load the wavelength coefficient functions '''
    path = os.path.join(options["outdir"], maskname)
    fn = os.path.join(path, "lambda_mask_coeffs_{0}.npy".format(fnum))

    ld = np.load(fn)
    return ld

def load_lambdaslit(fnum, maskname, band, options):
    ''' Load the wavelength coefficient functions '''
    path = os.path.join(options["outdir"], maskname)
    fn = os.path.join(path, "lambda_solution_{0}.fits".format(fnum))

    print fn

    return readfits(fn, options)

def writefits(img, maskname, fname, options, header=None, bs=None,
        overwrite=False):
    '''Convenience wrapper to write MOSFIRE drp-friendly FITS files'''

    hdu = pf.PrimaryHDU(img)
    path = os.path.join(options["outdir"], maskname)
    if not os.path.exists(path):
        print("Output directory '%s' does not exist. The DRP will attempt" 
                "to create this directory." % path)
        os.mkdir(path)

    fn = os.path.join(path, fname)

    if header is None: header = {"DRPVER": MOSFIRE.__version__}
    else: header.update("DRPVER", MOSFIRE.__version__)

    if header is not None:
        for k in header.keys():
            if hdu.header.has_key(k): continue

            k = k.rstrip()
            if len(k) <= 8:
                hdu.header.update(k, header[k])
            else:
                hdu.header.update("hierarch " + k, header[k])

    if overwrite:
        try: 
            os.remove(fn)
            print "Removed old file '%s'%" % (fn)
        except: pass

    print "Wrote to '%s'" % (fn)
    hdu.writeto(fn)



def readfits(path, use_bpm=False):
    '''Read a fits file from path and return a tuple of (header, data, 
    Target List, Science Slit List (SSL), Mechanical Slit List (MSL),
    Alignment Slit List (ASL)).'''
    hdulist = pf.open(path)
    header = hdulist[0].header
    data = hdulist[0].data

    if use_bpm:
        theBPM = badpixelmask()
        data = np.ma.masked_array(data, theBPM, fill_value=0)

    return (header, data)



def readheader(path):
    '''Reads a header (only) from a fits file'''

    return pf.getheader(path)

def read_drpfits(maskname, fname, options):
    '''Read a fits file written by the DRP'''

    path = os.path.join(options["outdir"], maskname, fname)

    hdulist = pf.open(path)
    output = []

    for hdu in hdulist:
        output.append(hdu.header)

        if hdu.header.has_key("DRPVER"):
            hasdrpver = hdu.header.has_key("DRPVER") 

            if hasdrpver:
                itsver = hdu.header["DRPVER"]
                if itsver != MOSFIRE.__version__:
                    raise Exception("The file requested '%s' uses DRP version %f "
                        "but the current DRP version is %f. There might be an "
                        "incompatibility" % (path, itsver, MOSFIRE.__version__))

            else:
                raise Exception("The file requested '%s' does not seem to be "
                        "the result of this DRP. This should never be the "
                        " case.")

        output.append(hdu.data)


    return output

def fname_to_path(fname, options):
    '''Take a filename like m120507_0123, parse date, and return full path'''
    months = {"01": "jan", "02": "feb", "03": "mar", "04": "apr", "05": "may",
        "06": "jun", "07": "jul", "08": "aug", "09": "sep", "10": "oct",
        "11": "nov", "12": "dec"}

    try:
        fdate = fname.split("m")[1][0:6]
        yr, mn, dy = "20" + fdate[0:2], fdate[2:4], int(fdate[4:6])
        month = months[mn]
    except:
        print "Could not parse date out of file name: %s" % (fname)

    path = os.path.join(options["indir"], yr + month + "%2.2i" % dy)
    if not os.path.exists(os.path.join(path, fname)):
        path = os.path.join(options["indir"], yr + month + "%2.2i" % (dy-1))

        if not os.path.exists(path):
            raise Exception("Could not find file '%s' in '%s' out of parsed "
                "%s, %s, %s" % (fname,
                options["indir"], yr, month, dy))

    return path


def readmosfits(fname, options, extension=None):
    '''Read a fits file written by MOSFIRE from path and return a tuple of 
    (header, data, Target List, Science Slit List (SSL), Mechanical Slit 
    List (MSL), Alignment Slit List (ASL)).
    
    Note, the extension is typically not used, only used if the detector server
    does not append slit extension.
    '''

    path = os.path.join(fname_to_path(fname, options), fname)

    hdulist = pf.open(path)
    header = hdulist[0].header
    data = hdulist[0].data

    theBPM = badpixelmask()
    data = np.ma.masked_array(data, theBPM)

    if extension is not None:
        hdulist = pf.open(extension)

    try:
        header = hdulist[0].header
        targs = hdulist[1].data
        ssl = hdulist[2].data
        msl = hdulist[3].data
        asl = hdulist[4].data
    except:
        raise Exception("Improper MOSFIRE FITS File: %s" % path)

    if np.abs(header["REGTMP1"] - 77) > .05:
        raise Exception("The temperature of the detector is %f where it "
                "should be 77.000 deg. Please notify Keck support staff." %
                header["REGTMP1"])

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
        if (header[s] != 'OK') and (header[s] != 'SETUP'):
            pos *= -1
        poss.append(pos)

    if len(poss) != CSU.numbars:
        raise CSU.MismatchError("Found %i bars instead of %i" % 
                (lens(poss), CSU.numbars))
        

    return np.array(poss)


def imcombine(filelist, out, options, bpmask=None, reject="none", nlow=None,
        nhigh=None):

    '''Convenience wrapper around IRAF task imcombine
    
    reject: none, minmax, sigclip, avsigclip, pclip'''

    #TODO: REMOVE Iraf and use python instead. STSCI Python has
    # A builtin routine.
    from pyraf import iraf
    iraf.images()

    path = fname_to_path(filelist[0], options)
    filelist = [os.path.join(path, "%s[0]" % f) for f in filelist]
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
