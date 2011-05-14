'''

===================
MOSFIRE Flat Fields 
===================




npk April 14th 2011

'''

import os
import time
import unittest

import numpy as np
import pylab as pl
import scipy, scipy.ndimage
import pyfits

import MOSFIRE
from MOSFIRE import Fit, IO, Options

__version__ = "19Apr2011"

reload(Fit)
reload(IO)
reload(Options)

#from IPython.Shell import IPShellEmbed
#start_shell = IPShellEmbed()

def tester(band):
    options = Options.flat
    path = os.path.join(options["indir"], "m110326_%4.4i.fits")

    if band == 'Y':
            handle_flats([path % 3242, path % 3243, path % 3244], 
                            "npk_calib3_q1700_pa_0", band, options)

    if band == 'H':
            handle_flats([path % 3248, path % 3249, path % 3250], 
                            "npk_calib3_q1700_pa_0", band, options)

    if band == 'K':
            handle_flats([path % 3251, path % 3252, path % 3253], 
                            "npk_calib3_q1700_pa_0", band, options)

    if band == 'J':
            handle_flats([path % 3245, path % 3246, path % 3247], 
                            "npk_calib3_q1700_pa_0", band, options)

def tester2(band):
    options = Options.flat
    path = os.path.join(options["outdir"], "npk_calib3_q1700_pa_0")
    res = np.load(os.path.join(path, "slit-edges_%s.npy" % band))
    save_ds9_edges(res, options)

def tester3(band):
    options = Options.flat
    path = os.path.join(options["outdir"], "npk_calib3_q1700_pa_0", 
                    "combflat_2d_%s.fits" % band)
    res = np.load(os.path.join(path, "slit-edges_%s.npy" % band))


def handle_flats(flatlist, maskname, band, options):
    '''
    handle_flats is the primary entry point to the Flats module.

    handle_flats takes a list of individual exposure FITS files and creates:
    1. A CRR, dark subtracted, pixel-response flat file.
    2. A set of polynomials that mark the edges of a slit

    Inputs:
    flatlist: 
    maskname: The name of a mask
    band: A string indicating the bandceil

    Outputs:

    file {maskname}/flat_2d_{band}.fits -- pixel response flat
    file {maskname}/edges.np
    '''

    try:
            os.mkdir(os.path.join(options["outdir"], maskname))
    except OSError:
            pass

    # Imcombine
    combine(flatlist, maskname, band, options)
    (header, data, bs) = IO.readmosfits(flatlist[0])
    path = os.path.join(options["outdir"], "npk_calib3_q1700_pa_0", 
                    "combflat_2d_%s.fits" % band)
    (header, data) = IO.readfits(path)

    # Edge Trace
    results = find_and_fit_edges(data, bs.ssl, options)
    results[-1]["maskname"] = maskname
    results[-1]["band"] = band
    np.save(os.path.join(options["outdir"], maskname, 
            "slit-edges_%s" % band), results)
    save_ds9_edges(results, options)

    # Generate Flat
    out = os.path.join(options["outdir"], maskname, 
                    "pixelflat_2d_%s.fits" % (band))
    make_pixel_flat(data, results, options, out)

def make_pixel_flat(data, results, options, outfile):
    '''
    Convert a flat image into a flat field
    '''

    def pixel_min(y): return np.floor(np.min(y))
    def pixel_max(y): return np.ceil(np.max(y))

    def collapse_flat_box(dat):
        '''Collapse data to the spectral axis (0)'''
        v = np.median(dat, axis=0).ravel()

        return v

    tick = time.clock()
    xs = np.arange(Options.npix)
    flat = np.ones(shape=(Options.npix, Options.npix))

    for result in results[0:-1]:
        print result["Target_Name"]
        bf = result["bottom"]
        tf = result["top"]

        top = pixel_min(tf(xs))
        bottom = pixel_max(bf(xs))

        print "Bounding top/bottom: %i/%i" % (bottom, top)

        v = collapse_flat_box(data[bottom:top,:])
        ok = v > 180

        v = np.poly1d(np.polyfit(xs[ok],v[ok],
            options['flat-field-order']))(xs).ravel()

        for i in np.arange(bottom-1, top+1):
                flat[i,:] = v

    for r in range(len(results)-2):
        print r
        first = results[r]
        second = results[r+1]

        tf = first["bottom"]
        bf = second["top"]

        for i in range(Options.npix):
            top = np.floor(tf(i))
            bottom = np.ceil(bf(i))
            
            data[top:bottom, i] = flat[top:bottom,i]

    lowsn = data<225
    flat[lowsn] = data[lowsn]

    hdu = pyfits.PrimaryHDU((data/flat).astype(np.float32))
    if os.path.exists(outfile):
            os.remove(outfile)
    hdu.writeto(outfile)

    print "Pixel flat took %i s" % (time.clock()-tick)

def save_ds9_edges(results, options):
    '''
    Create a .ds9 file that saves the fit slit edge positions determined
    by find_and_fit_edges
    '''

    ds9 = ''

    W = Options.npix
    delt = Options.npix/30.
    
    S = 1
    for i in range(len(results) - 1):
        res = results[i]

        top = res["top"]
        bottom = res["bottom"]
        for i in np.arange(W/delt):
            x = delt * i
            sx = x + 1
            ex = x + delt + 1

            sy = top(sx) + 1
            ey = top(ex) + 1
            ds9 += "line(%f, %f, %f, %f) # fixed=1 edit=0 move=0 rotate=0 delete=0" % (sx, sy, ex, ey)

            if i == W/2:
                    ds9 += " # text={S%2.0i (%s)}" % (S, 
                                    res["Target_Name"])

            ds9 += "\n"

            sy = bottom(sx) + 1
            ey = bottom(ex) + 1
            ds9 += "line(%f, %f, %f, %f) # fixed=1 edit=0 move=0 rotate=0 delete=0 color=blue\n" % (sx, sy, ex, ey)

    band = results[-1]["band"]
    fn = os.path.join(options["outdir"], results[-1]["maskname"], 
                    "slit-edges_%s.ds9" % band)
    try:
            f = open(fn,'w')
            f.write(ds9)
            f.close()
    except IOError:
            print "IO Error"
            raise
    except:
            raise

def combine(flatlist, maskname, band, options):
    '''
    combine list of flats into a flat file'''

    out = os.path.join(options["outdir"], maskname, "combflat_2d_%s.fits" 
                    % (band))
    if os.path.exists(out):
            os.remove(out)

    IO.imcombine(flatlist, out, reject="sigclip")


def find_edge_pair(data, y, roi_width):
    '''
    find_edge_pair finds the edge of a slit pair in a flat

    data[2048x2048]: a well illuminated flat field [DN]
    y: guess of slit edge position [pix]

    Moves along the edge of a slit image
            - At each location along the slit edge, determines
            the position of the demarcations between two slits

    Outputs:
    xposs []: Array of x positions along the slit edge [pix]
    yposs []: The fitted y positions of the "top" edge of the slit [pix]
    widths []: The fitted delta from the top edge of the bottom [pix]
    scatters []: The amount of light between slits


    The procedure is as follows
    1: starting from a guess spatial position (parameter y), march
        along the spectral direction in some chunk of pixels
    2: At each spectral location, construct a cross cut across the
        spatial direction; select_roi is used for this.
    3: Fit a two-sided error function Fit.residual_disjoint_pair
        on the vertical cross cut derived in step 2.
    4: If the fit fails, store it in the missing list, otherwise
        store the position and fitted values in xposs, yposs, and widths.
    5: In the vertical cross-cut, there is a minimum value. This minimum
        value is stored as a measure of scattered light.

    Another procedure is used to fit polynomials to these fitted values.
    '''

    def select_roi(data, roi_width):
        v = data[y-roi_width:y+roi_width, xp-2:xp+2]
        v = np.median(v, axis=1) # Axis = 1 is spatial direction

        return v


    yposs = []
    widths = []
    xposs = []
    xposs_missing = []
    scatters = []

    #1 
    rng = np.linspace(10, 2040, 100).astype(np.int)
    for i in rng:
        xp = i
        #2
        v = select_roi(data, roi_width)
        xs = np.arange(len(v))

        # TODO: The number 200 is hardcoded and essentially made up.
        # A smarter piece of code belongs here.
        if np.median(v) < 200:
            xposs_missing.append(xp)
            continue

        #3
        ff = Fit.do_fit(v, residual_fun=Fit.residual_disjoint_pair)
        fit_ok = 0 < ff[4] < 4
        if fit_ok:
            (sigma, offset, mult1, mult2, add, width) = ff[0]

            xposs.append(xp)
            yposs.append(offset - roi_width)
            widths.append(width)

            between = offset + width/2
            if 0 < between < len(v)-1:
                scatters.append(np.min(v[between-2:between+2])) # 5
            else:
                scatters.append(np.nan)
        else:
            xposs_missing.append(xp)
            print "Skipping: %i" % (xp)

    
    return map(np.array, (xposs, xposs_missing, yposs, widths, scatters))

def fit_edge_poly(xposs, xposs_missing, yposs, widths, order):
    '''
    fit_edge_poly fits a polynomial to the measured slit edges.
    This polynomial is used to extract spectra.

    fit_edge_poly computes a parabola, and fills in missing data with a 
    parabola

    input-
    xposs, yposs [N]: The x and y positions of the slit edge [pix]
    widths [N]: the offset from end of one slit to beginning of another 
            [pix]
    order: the polynomial order
    '''

    # First fit low order polynomial to fill in missing data
    fun = np.poly1d(Fit.polyfit_clip(xposs, yposs, 2))
    wfun = np.poly1d(Fit.polyfit_clip(xposs, widths, 2))

    xposs = np.append(xposs, xposs_missing)
    yposs = np.append(yposs, fun(xposs_missing))
    widths = np.append(widths, wfun(xposs_missing))


    # Remove any fits that deviate wildly from the 2nd order polynomial
    ok = np.abs(yposs - fun(xposs)) < 5
    if not ok.any():
            raise Exception("Flat is not well illuminated? Cannot find edges")

    # Now refit to proper order
    fun = np.poly1d(Fit.polyfit_clip(xposs[ok], yposs[ok], order))
    wfun = np.poly1d(Fit.polyfit_clip(xposs[ok], widths[ok], order))
    res = fun(xposs[ok]) - yposs[ok]
    sd = np.std(res)
    ok = np.abs(res) < 2*sd

    return (fun, wfun, res, sd, ok)



def find_and_fit_edges(data, ssl, options):
    '''
    Given a flat field image, find_and_fit_edges determines the position
    of all slits.

    The function works by starting with a guess at the location for a slit
    edge in the spatial direction(options["first-slit-edge"]). 
    
    Starting from the guess, find_edge_pair works out in either direction, 
    measuring the position of the (e.g.) bottom of slit 1 and top of slit 2:


    ------ pixel y value = 2048

    Slit 1 data

    ------ (bottom)
    deadband
    ------ (top)

    Slit N pixel data ....

    ------- (bottom) pixel = 0

    --------------------------------> Spectral direction


    1. At the top of the flat, the slit edge is defined to be a pixel value
    2. The code guesses the position of the bottom of the slit, and runs
            find_edge_pair to measure slit edge locations.
    3. A low-order polynomial is fit to the edge locations with
            fit_edge_poly
    4. The top and bottom of the current slit, is stored into the
            result list.
    5. The top of the next slit is stored temporarily for the next
            iteration of the for loop.
    6. At the bottom of the flat, the slit edge is defined to be pixel 4.


    options:
    options["edge-order"] -- The order of the polynomial [pixels] edge.
    options["edge-fit-width"] -- The length [pixels] of the edge to 
            fit over

    '''
    y = 2028
    toc = 0
    slits = []
    top = [0., np.float(Options.npix)]
    numslits = np.round(np.array(ssl.field("Slit_length"), 
        dtype=np.float) / 7.02)


    results = []
    result = {}
    result["Target_Name"] = ssl[0].field("Target_Name")

    # 1
    result["top"] = np.poly1d([2027])

    for target in range(len(ssl) - 1):
        y -= 44.25 * numslits[target]
        print "-------------==========================-------------"
        print "Finding Slit Edges for %s starting at %4.0i" % (
                        ssl[target].field("Target_Name"), y)
        tock = time.clock()

        (xposs, xposs_missing, yposs, widths, scatters) = \
                        find_edge_pair(data, y, 
                                        options["edge-fit-width"])

        (fun, wfun, res, sd, ok) = fit_edge_poly(xposs, 
                        xposs_missing, yposs, widths, 
                        options["edge-order"])

        bottom = fun.c.copy() 
        top = wfun.c.copy() + fun.c.copy()
        bottom[-1] += y 
        top[-1] += y 

        #4
        result["xposs"] = xposs
        result["yposs"] = yposs
        result["bottom"] = np.poly1d(bottom)
        result["sd"] = sd
        result["ok"] = ok
        result["scatter"] = scatters

        results.append(result)

        #5
        result = {}
        result["Target_Name"] = ssl[target].field("Target_Name")
        result["top"] = np.poly1d(top)


        print fun
        print "Clipped Resid Sigma: %5.3f P2V: %5.3f" % (
                        np.std(res[ok]), res[ok].max()-res[ok].min())

        tic = time.clock()

        print " .... %4.2f s elapsed." % (tic - tock)
        print

    #6
    result["bottom"] = np.poly1d([3])
    results.append(result)

    results.append({"version": options["version"]})

    return results


class TestFlatsFunctions(unittest.TestCase):
    def setUp(self):

            pass

    def test_trace_edge(self):
            (header, data1, targs, ssl, msl, asl) = \
                            IO.readfits_all("/users/npk/desktop/c9/m110326_3242.fits")
            data = data1

            ssl = ssl[ssl.field("Slit_Number") != ' ']
            numslits = np.round(np.array(ssl["Slit_length"], 
                    dtype=np.float) / 7.02)

            for i in range(len(ssl)):
                    print ssl[i].field("Target_Name"), numslits[i]

if __name__ == '__main__':
    unittest.main()
