'''

-------------------------------------------------------------------------------
Created on April 23 2011 by npk, modified until April 24th

This file is used to test H-band wavelength calibration code. This experiment
is a continuation of that in drp5.py

You can find extensive notes in MOSFIRE DRP notebook #1 pg. 32 - 37

Lots of changes between drp5.py and drp6.py:
    - Fit.* functions now are better behaved:
        . fit sinbeta rather than use small angle approximation
        . More reasonable limits on parameters
        . Use pixel value 1024 as "zero point"; this provides more
            sane values for the incoming angle.
    - Created best guess functions for parameters alpha, beta, gamma, 
        and delta.
    - Learned that code does not deal well with blends or saturated lines. I
        will work on these two problems in drp7
'''

import numpy as np
import pylab as pl
import scipy as sp

from MOSFIRE import IO, Fit, Bspline

reload(IO)
reload(Fit)
reload(Bspline)

np.set_printoptions(precision=3)


# use the following file and id'd spectrum to guess


(header, data, targs, ssl, msl, asl) = IO.readfits_all(
    "/users/npk/desktop/c9/m110323_2718.fits")
band = 'H'


pl.ion()


def find_peaks(ll, spec, lines):
    xs = []
    sxs = []
    sigmas = []
    pix = np.arange(2048.)

    for lam in lines:
        f = 0.9985
        roi = (f*lam < ll) & (ll < lam/f)

        if not roi.any(): 
            xs.append(0.0)
            sxs.append(np.inf)
            sigmas.append(0.0)
            continue

        lsf = Fit.mpfitpeak(pix[roi], spec[roi], 
            error=np.sqrt(np.abs(spec[roi])))
    
        if lsf.perror is None:
            xs.append(0.0)
            sxs.append(np.inf)
            sigmas.append(0.0)
            continue

        if lsf.status < 0:
            xs.append(0.0)
            sxs.append(np.inf)
            sigmas.append(0.0)
            continue

        mnpix = np.min(pix[roi])
        mxpix = np.max(pix[roi])

        if (mnpix + 4) > lsf.params[1] < (mxpix-4):
            xs.append(0.0)
            sxs.append(np.inf)
            sigmas.append(0.0)
            continue

        if mnpix < 7:
            xs.append(0.0)
            sxs.append(np.inf)
            sigmas.append(0.0)
            continue

        if mxpix > 2040:
            xs.append(0.)
            sxs.append(np.inf)
            sigmas.append(0.0)
            continue

        xs.append(lsf.params[1])
        sxs.append(lsf.perror[1])
        sigmas.append(lsf.params[2])


    return [xs, sxs, sigmas]
    

def dofit(xs, slambdas, lines, alpha, sinbeta, gamma, delta, band):
    xs = np.array(xs)
    ok = np.isfinite(slambdas)
    lsf = Fit.do_fit_wavelengths(xs[ok], lines[ok], alpha, sinbeta, 
            gamma, delta, band, error=slambdas[ok])
    (order, alpha, sinbeta, gamma, delta) = lsf.params

    print "order alpha    sinbeta   gamma  delta"
    print "%1.0i %5.7f %3.5f %3.2e %5.2f" % (order, alpha, sinbeta, 
            gamma, delta)

    pix = np.arange(2048.)
    ll = Fit.wavelength_model(lsf.params, pix)

    return [np.abs((
        Fit.wavelength_model(lsf.params, xs[ok]) - lines[ok]))*1e4, 
            lsf.params]



def fit_2d_spec(data, pos, params, band):
    global DRAW2d
    ar_h_lines = np.array([1.465435, 1.474317, 1.505062, 1.517684, 1.530607, 
        1.533334, 1.540685, 1.590403, 1.599386, 1.618444, 1.644107,
        1.674465, 1.744967, 1.791961])

    Ar_K_lines = np.array([1.982291, 1.997118, 2.032256, 2.057443, 
        2.062186, 2.099184, 2.133871, 2.15409, 2.204558, 2.208321,
        2.313952, 2.385154])

    lines = ar_h_lines


    # do_fits is a closure and appens results
    def do_fits(params, pos, dy):
        (alpha, sinbeta, gamma, delta) = params[1:]
        thispos = pos - dy
        yposs.append(thispos)
        spec = np.median(data[thispos-3:thispos+3, :], axis=0)
        ll = Fit.wavelength_model(params, pix)

        [xs, sxs, sigmas] = find_peaks(ll, spec, lines)
        slambdas = d/order * np.array(sxs)

        [deltas, params] = dofit(xs, slambdas, lines, alpha, sinbeta,
                gamma, delta, band)
        results.append(np.array(params))
        resdelt.append(deltas)

    bmap = {"Y": 6, "J": 5, "H": 4, "K": 3}
    order = bmap[band]

    d = 1e3/110.5

    yposs = []
    results = []
    resdelt = []
    pix = np.arange(2048.)
    
    print "-------------------"
    initial = params
    # Search up
    for dy in range(0, 15):
        do_fits(params, pos, dy)

    # Search down
    for dy in range(-15,0):
        do_fits(params, pos, dy)

    return [yposs, np.array(results)]

DRAW = False
def fit_spec(data, pos, alpha, sinbeta, gamma, delta, band):
    global DRAW
    ar_h_lines = np.array([1.465435, 1.474317, 1.505062, 1.517684, 1.530607, 
        1.533334, 1.540685, 1.590403, 1.599386, 1.618444, 1.644107,
        1.674465, 1.744967, 1.791961])

    Ar_K_lines = np.array([1.982291, 1.997118, 2.032256, 2.057443, 
        2.062186, 2.099184, 2.133871, 2.15409, 2.204558, 2.208321,
        2.313952, 2.385154])

    lines = ar_h_lines

    bmap = {"Y": 6, "J": 5, "H": 4, "K": 3}
    order = bmap[band]
    d = 1e3/110.5

    pix = np.arange(2048)
    ll = Fit.wavelength_model([order, alpha, sinbeta, gamma, delta], pix)
    spec = np.median(data[pos-3:pos+3, :], axis=0)

    [xs, sxs, sigmas] = find_peaks(ll, spec, lines)

    slambdas = d/order * np.array(sxs)

    [deltas, params] = dofit(xs, slambdas, lines, alpha, sinbeta, gamma, 
            delta, band)

    return (params[1], params[2], params[3], params[4],
            np.median(deltas), np.std(deltas), pixel, sigmas)


alpha_pixel = np.poly1d([1.752e-14, 3.979e-9, -8.216e-6, 1.001])
if band == 'Y' or band == 'J':
    sinbeta_position = np.poly1d([0.0239, 36.2])
    sinbeta_pixel = np.poly1d([-2.578e-7, 0.00054, -0.2365])
    gamma_pixel = np.poly1d([-1.19e-19, 2.442e-16, 4.181e-13])
elif band == 'H' or band == 'K':
    sinbeta_position = np.poly1d([0.02365, 38.1])
    sinbeta_pixel = np.poly1d([-2.661e-7, 0.000553, -0.1956])
    gamma_pixel = np.poly1d([-1.19e-19, 2.442e-16, 6.181e-13])

delta_pixel = np.poly1d([-0.0158, 971.9])

def sinbeta_predict(position, pixel): return sinbeta_position(position) + \
    sinbeta_pixel(pixel)


pl.ion()

for i in range(4): 
    pl.figure(i)
    pl.clf()

bmap = {"Y": 6, "J": 5, "H": 4, "K": 3}
order = bmap[band]
np.seterr(all="ignore")
if True:
    positions = []
    alphas = []
    sinbetas = []
    gammas = []
    deltas = []
    MADs = []
    SDs = []
    pixels = []
    y0 = 2004.9
    for slit in range(35,36):
        os = "b%2.2ipos" % (slit*2-1)
        es = "b%2.2ipos" % (slit*2)
        csupos = (header[os] + header[es])/2.
        pixel = np.int(y0 - (slit-1) * 44.22)

        alpha = alpha_pixel(pixel)
        sinbeta = sinbeta_predict(csupos, pixel)
        gamma = gamma_pixel(pixel)
        delta = delta_pixel(pixel)


        print "------- Slit %i (%4.0i) @ %3.5f" % (slit, pixel, sinbeta)
        (newa, newb, newg, newd, MAD, SD, pixel, sigmas) = \
            fit_spec(data, pixel, alpha, sinbeta, gamma, delta, band)

        dalpha = dsinbeta = 0
        if MAD > .2:
            print "Searching for a better fit to improve MAD of %f" % MAD
            alphadirection = 1.0
            sinbetadirection = 1.0
            dsinbetas = np.sort(np.abs(np.linspace(0.15, -0.15, 35)))
            dalphas = np.array([0])
            for dalpha in dalphas:
                dalpha *= alphadirection
                alphadirection *= -1
                for dsinbeta in dsinbetas:
                    dsinbeta *= sinbetadirection
                    sinbetadirection *= -1

                    (newa, newb, newg, newd, MAD, SD, pixel, sigmas) = \
                            fit_spec(data, pixel, alpha+dalpha, 
                                    sinbeta+dsinbeta, gamma, delta, band)
                    if MAD < .2:
                        print "breaking at dalpha: %f, dsinbeta: %f" % (
                                dalpha, dsinbeta)
                        break
            
                if MAD < .2:
                    break

        [ypix, ress] = fit_2d_spec(data, pixel, [order, alpha+dalpha,
            sinbeta+dsinbeta, gamma, delta], band)

        1/0

        print
        print "RESULT: %3.3f %2.2e | %2.5f %5.2f" % ( csupos, MAD, newa, newb )
        print "====="
        positions.append(csupos)
        alphas.append(newa)
        sinbetas.append(newb)
        gammas.append(newg)
        deltas.append(newd)
        MADs.append(MAD)
        SDs.append(SD)
        pixels.append(pixel)

    alphas = np.array(alphas)
    sinbetas = np.array(sinbetas)
    gammas = np.array(gammas)
    deltas = np.array(deltas)

    MADs = np.array(MADs)
    SDs = np.array(SDs)
    positions = np.array(positions)
    pixels = np.array(pixels)
