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

import pdb

import MOSFIRE
from MOSFIRE import Fit, IO, Options, CSU, Wavelength, Filters, Detector

__version__ = 0.1

#from IPython.Shell import IPShellEmbed
#start_shell = IPShellEmbed()

def handle_flats(flatlist, maskname, band, options, extension=None):
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

    tick = time.time()

    path = os.path.join(options["outdir"], maskname)
    try: os.mkdir(path)
    except OSError: pass

    # Check
    bpos = np.ones(92) * -1
    for fname in flatlist:

        hdr, dat, bs = IO.readmosfits(fname, options)
        try: bs0
        except: bs0 = bs

        if np.any(bs0.pos != bs.pos):
            raise Exception("Barsets do not seem to match")

        if hdr["filter"] != band:
            raise Exception("Filter name %s does not match header filter name "
                    "%s in file %s" % (band, hdr["filter"], fname))
        for i in xrange(len(bpos)):
            b = hdr["B{0:02d}POS".format(i+1)]
            if bpos[i] == -1:
                bpos[i] = b
            else:
                if bpos[i] != b:
                    raise Exception("Bar positions are not all the same in "
                            "this set of flat files")
    bs = bs0
    # Imcombine
    print "Attempting to combine: ", flatlist
    combine(flatlist, maskname, band, options)

    print "Combined '%s' to '%s'" % (flatlist, maskname)
    path = os.path.join(options["outdir"], maskname,
                    "combflat_2d_%s.fits" % band)
    (header, data) = IO.readfits(path, use_bpm=True)

    print "Flat written to %s" % path


    # Edge Trace
    results = find_and_fit_edges(data, header, bs, options)
    results[-1]["maskname"] = maskname
    results[-1]["band"] = band
    np.save(os.path.join(options["outdir"], maskname, 
            "slit-edges_{0}".format(band)), results)
    save_ds9_edges(results, options)

    # Generate Flat
    out = os.path.join(options["outdir"], maskname, 
                    "pixelflat_2d_%s.fits" % (band))
    make_pixel_flat(data, results, options, out, flatlist)
    print "Pixel flat took {0:6.4} s".format(time.time()-tick)

def make_pixel_flat(data, results, options, outfile, inputs):
    '''
    Convert a flat image into a flat field
    '''

    def pixel_min(y): return np.floor(np.min(y))
    def pixel_max(y): return np.ceil(np.max(y))

    def collapse_flat_box(dat):
        '''Collapse data to the spectral axis (0)'''
        v = np.median(dat, axis=0).ravel()

        return v

    flat = np.ones(shape=Detector.npix)

    hdu = pyfits.PrimaryHDU((data/flat).astype(np.float32))
    hdu.header.update("version", __version__, "DRP version")
    i = 0
    for flatname in inputs:
        nm = flatname.split("/")[-1]
        hdu.header.update("infile%2.2i" % i, nm)
        i += 1

    slitno = 0
    for result in results[0:-1]:
        slitno += 1

        hdu.header.update("targ%2.2i" % slitno, result["Target_Name"])

        bf = result["bottom"]
        tf = result["top"]
        try:
            hpps = result["hpps"]
        except:
            print "No half power points for this slit"
            hpps = [0, Detector.npix[0]]

        xs = np.arange(hpps[0], hpps[1])

        top = pixel_min(tf(xs))
        bottom = pixel_max(bf(xs))

        hdu.header.update("top%2.2i" % slitno, top)
        hdu.header.update("bottom%2.2i" % slitno, bottom)

        print "%s] Bounding top/bottom: %i/%i" % (result["Target_Name"],
                bottom, top)

        v = collapse_flat_box(data[bottom:top,hpps[0]:hpps[1]])

        x2048 = np.arange(Options.npix)
        v = np.poly1d(np.polyfit(xs,v,
            options['flat-field-order']))(x2048).ravel()

        for i in np.arange(bottom-1, top+1):
            flat[i,:] = v

    print "Producing Pixel Flat..."
    for r in range(len(results)-2):
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
    hdu.data = (data/flat).astype(np.float32)
    hdu.data = hdu.data.filled(1)
    if os.path.exists(outfile):
            os.remove(outfile)
    hdu.writeto(outfile)


def save_ds9_edges(results, options):
    '''
    Create a ds9 file that saves the fit slit edge positions determined
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
            if i == 10: txt=res["Target_Name"]
            else: txt=""

            ds9 += "line(%f, %f, %f, %f) # fixed=1 edit=0 move=0 rotate=0 delete=0 text={%s}\n" % (sx, sy, ex, ey, txt)
            ds9 += "line(%f, %f, %f, %f) # fixed=1 edit=0 move=0 rotate=0 delete=0\n" % (sx, sy, ex, ey)

            if i == W/2:
                    ds9 += " # text={S%2.0i (%s)}" % (S, 
                                    res["Target_Name"])

            ds9 += "\n"

            sy = bottom(sx) + 1
            ey = bottom(ex) + 1
            ds9 += "line(%f, %f, %f, %f) # fixed=1 edit=0 move=0 rotate=0 delete=0 color=blue\n" % (sx, sy, ex, ey)

        # Vertical line indicating half power points
        try:
            hpps = res["hpps"]
            sx = hpps[0] ; ex = hpps[0]
            sy = bottom(sx) ; ey = top(sx)
            ds9 += "line(%f, %f, %f, %f) # fixed=1 edit=0 move=0 rotate=0 delete=0\n" % (sx, sy, ex, ey)

            sx = hpps[1] ; ex = hpps[1]
            sy = bottom(sx) ; ey = top(sx)
            ds9 += "line(%f, %f, %f, %f) # fixed=1 edit=0 move=0 rotate=0 delete=0\n" % (sx, sy, ex, ey)
        except:
            continue
        
    band = results[-1]["band"]
    fn = os.path.join(options["outdir"], results[-1]["maskname"], 
                    "slit-edges_%s.reg" % band)
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

    IO.imcombine(flatlist, out, options, reject="minmax", nlow=1, nhigh=1)


def find_edge_pair(data, y, roi_width, hpps):
    '''
    find_edge_pair finds the edge of a slit pair in a flat

    data[2048x2048]: a well illuminated flat field [DN]
    y: guess of slit edge position [pix]
    hpps[2]: estimate of of spectral half power points [pix]

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
    rng = np.linspace(10, 2040, 50).astype(np.int)
    for i in rng:
        xp = i
        #2
        v = select_roi(data, roi_width)
        xs = np.arange(len(v))

        # TODO: The number 450 below is hardcoded and essentially made up.
        # A smarter piece of code belongs here.
        if (np.median(v) < 450) or (xp < hpps[0]) or (xp > hpps[1]):
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
                start = np.max([0, between-2])
                stop = np.min([len(v),between+2])
                scatters.append(np.min(v[start:stop])) # 5

                if False:
                    pl.figure(2)
                    pl.clf()
                    tmppix = np.arange(y-roi_width, y+roi_width)
                    tmpx = np.arange(len(v))
                    pl.scatter(tmppix, v)
                    pl.plot(tmppix, Fit.fit_disjoint_pair(ff[0], tmpx))
                    pl.draw()

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
    ok = np.abs(yposs - fun(xposs)) < 1
    if not ok.any():
            raise Exception("Flat is not well illuminated? Cannot find edges")

    # Now refit to user requested order
    fun = np.poly1d(Fit.polyfit_clip(xposs[ok], yposs[ok], order))
    wfun = np.poly1d(Fit.polyfit_clip(xposs[ok], widths[ok], order))
    res = fun(xposs[ok]) - yposs[ok]
    sd = np.std(res)
    ok = np.abs(res) < 2*sd


    # Check to see if the slit edge funciton is sane, 
    # if it's not, then we fix it.
    if np.abs(fun(0) - fun(2048)) > 15:
        print "Forcing a horizontal slit edge"
        fun = np.poly1d(np.median(yposs[ok]))
        wfun = np.poly1d(np.median(widths[ok]))


    return (fun, wfun, res, sd, ok)



def find_and_fit_edges(data, header, bs, options):
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

    # TODO: move hardcoded values into Options.py
    # y is the location to start
    y = 2034
    DY = 44.25

    toc = 0
    ssl = bs.ssl

    slits = []

    top = [0., np.float(Options.npix)]
    numslits = np.round(np.array(ssl["Slit_length"], 
        dtype=np.float) / 7.6)


    if np.sum(numslits) != CSU.numslits:
        raise Exception("The number of allocated CSU slits (%i) does not match "
                " the number of possible slits (%i)." % (np.sum(numslits),
                    CSU.numslits))
    results = []
    result = {}

    result["Target_Name"] = ssl[0]["Target_Name"]

    # 1
    result["top"] = np.poly1d([y])

    slitno = 1

    for target in xrange(len(ssl) - 1):

        y -= DY * numslits[target]

        slitno += numslits[target]

        print("%2.2i] Finding Slit Edges for %s starting at %4.0i. Slit "
                "composed of %i CSU slits" % ( target,
                    ssl[target]["Target_Name"], y, numslits[target]))

        tock = time.clock()

        hpps = Wavelength.estimate_half_power_points(slitno, header, bs)

        (xposs, xposs_missing, yposs, widths, scatters) = \
                        find_edge_pair(data, y, 
                                        options["edge-fit-width"], hpps)

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
        result["Target_Name"] = ssl[target+1]["Target_Name"]
        result["top"] = np.poly1d(top)

        hpps = Wavelength.estimate_half_power_points(slitno, header, bs)
        result["hpps"] = hpps

        tic = time.clock()

        print("    Clipped Resid Sigma: %5.3f P2V: %5.3f .... %4.2f s elapsed" % 
            (np.std(res[ok]), res[ok].max()-res[ok].min(),tic-tock))


    #6
    result["bottom"] = np.poly1d([3])
    result["hpps"] = hpps
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

            ssl = ssl[ssl["Slit_Number"] != ' ']
            numslits = np.round(np.array(ssl["Slit_length"], 
                    dtype=np.float) / 7.02)

            for i in range(len(ssl)):
                    print ssl[i]["Target_Name"], numslits[i]

if __name__ == '__main__':
    unittest.main()
