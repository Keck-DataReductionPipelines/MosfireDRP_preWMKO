'''

Created on April 21 2011 by npk
Modified until April 23 2011 by npk

This file is used to test H-band wavelength calibration code

You can find extensive notes in MOSFIRE DRP notebook #1 pg. 25


NOTE THAT THIS FILE AND ITS NOTES ARE OBSOLETE. SEE NOTEBOOK
AND THE FILE drp6.py.

I determined that to first order the physical model for wavelengths
works well. The physical model looks like:

    lambda(pixel | alpha, beta) = alpha/(m sigma) 
        (sin(pixelsize/focallength pixel) - beta)

where alpha and beta are fit. m, sigma, pixelsize, and focal length are
known constants of MOSFIRE:
    m - dispersion order
    sigma - Grating is 110 lines per mm
    pixelsize - 18 micron
    focal length - 250 mm (camera)




Based on H band spectra taken in cooldown 9 I determined predicted values for
alpha and beta based on the CSU slit position (in mm), as well as pixel position
on the detector (in dimensionless pixels).

alpha_pixel maps pixel number (spatial direction) into an alpha (~ 1.0)
alpha_pixel = np.poly1d([3.2e-9, -6.4e-6, 1.0033])

beta_position maps CSU position mm to beta
beta_position = np.poly1d([0.02332, 33.8069])

beta_pixel maps pixel number (spatial direction) into a detla value added
    to the beta_position function.
beta_pixel = np.poly1d([-2.95e-7, 0.0006214, -0.2254])
def beta_predict(position, pixel): return beta_position(position) + \
    beta_pixel(pixel)


All of these assume y0 = 2004.9, and delta slit lengths of 44.25 pixels.

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
    "/users/npk/desktop/c9/m110323_2717.fits")

csupos = (header["b44pos"]+header["b43pos"])/2.


pl.ion()

DRAW = True
def fit_spec(data, pos, alpha, beta):
    global DRAW
    lines = np.array([1.49339, 1.49904, 1.51442, 1.51950, 1.52349, 1.53524, 
        1.54118,  1.6027, 1.62728, 1.64097, 1.71666])

    order = 4
    d = 1e3/110.5
    pix = np.arange(2048.)

    ll = alpha * d/order * (np.sin(18./(250e3) * pix) + np.radians(beta))

    spec = np.median(data[pos-3:pos+3, :], axis=0)
    if DRAW:
        pl.figure(2)
        pl.clf()
        pl.figure(3)
        pl.clf()
        pl.figure(1)
        pl.clf()
        pl.plot(ll, spec, '-+')
        for lam in lines:
            pl.axvline(lam, color='r')


        pl.title("Pixel pos: %i"%  pos)
        pl.figure(2)
        pl.clf()

    xs = []
    sxs = []
    pks = []
    rats = []
    for lam in lines:
        roi = ((lam-.002) < ll) & (ll < (lam+0.002))

        if not roi.any(): 
            xs.append(0.0)
            sxs.append(np.inf)
            continue

        lsf = Fit.mpfitpeak(pix[roi], spec[roi], 
            error=np.sqrt(np.abs(spec[roi])))
    
        if lsf.perror is None:
            xs.append(0.0)
            sxs.append(np.inf)
            continue

        if lsf.status < 0:
            xs.append(0.0)
            sxs.append(np.inf)
            continue

        #if DRAW:
        if False:
            pl.plot(pix[roi], spec[roi])
            pl.plot(pix[roi], spec[roi], 'o')
            pl.plot(pix[roi], Fit.gaussian(lsf.params, pix[roi]), 'r--')

        xs.append(lsf.params[1])
        sxs.append(lsf.perror[1])

        if False: 
            pl.axvline(lsf.params[1])

    if DRAW:
        pl.draw()

    def dofit(xs, sxs, alpha, beta, band):
        xs = np.array(xs)
        sxs = np.array(sxs)/1000.
        lsf = Fit.do_fit_wavelengths(xs, lines, alpha, beta, band, error=sxs)
        (order, alpha, beta, gamma, delta) = lsf.params

        print "alpha    beta   gamma  delta"
        print "%5.7f %3.5f %3.2e %5.2f" % (alpha, beta, gamma, delta)

        ll = Fit.wavelength_model(lsf.params, pix)
        if False:
            pl.figure(3)
            pl.clf()
            pl.plot(ll, spec)
            pl.title("Pixel %4.4f" % pos)
            for lam in lines:
                pl.axvline(lam, color='r')

            pl.draw()
    
        return [np.abs((Fit.wavelength_model(lsf.params, xs) - lines))/lines, 
                lsf.params]

    [delta, params] = dofit(xs, sxs, alpha, beta, 'H')
    return (params[1], params[2], np.median(delta), np.std(delta), pixel)


alpha_pixel = np.poly1d([3.2e-9, -6.4e-6, 1.0033])
beta_position = np.poly1d([0.02332, 33.8069])
beta_pixel = np.poly1d([-2.95e-7, 0.0006214, -0.2254])
def beta_predict(position, pixel): return beta_position(position) + \
    beta_pixel(pixel)

if False:
    (newa, newb, MAD) = fit_spec(data, 1024, alpha, beta)
    print "%4.4f %4.4f %4.4f" % (newa, newb, MAD)

pl.ion()

np.seterr(all="ignore")
if True:
    positions = []
    alphas = []
    betas = []
    MADs = []
    SDs = []
    pixels = []
    y0 = 2004.9
    for slit in range(1,47):
        os = "b%2.2ipos" % (slit*2-1)
        es = "b%2.2ipos" % (slit*2)
        csupos = (header[os] + header[es])/2.
        pixel = np.int(y0 - (slit-1) * 44.25)

        alpha = alpha_pixel(pixel)
        beta = beta_predict(csupos, pixel)

        print "------- Slit %i (%4.0i) @ %3.5f" % (slit, pixel, beta)
        (newa, newb, MAD, SD, pixel) = fit_spec(data, pixel, alpha, beta)
        dbeta = 0.0
        if MAD > 1e-4:
            print "Searching for Beta"
            direction = 1.0
            dbetas = np.sort(np.abs(np.linspace(0.3, -0.3, 20)))
            for dbeta in dbetas:
                dbeta *= direction
                direction *= -1
                (newa, newb, MAD, SD, pixel) = fit_spec(data, pixel, alpha, 
                        beta+dbeta)
                print "Beta: %4.2f yields MAD: %2.2e" % (beta+dbeta, MAD)
                if MAD < 1e-4:
                    print "breaking at dbeta: %f" % dbeta
                    break
            


        print
        print "RESULT: %3.3f %2.2e | %2.5f %5.2f" % ( csupos, MAD, newa, newb )
        print "====="
        positions.append(csupos)
        alphas.append(newa)
        betas.append(newb)
        MADs.append(MAD)
        SDs.append(SD)
        pixels.append(pixel)

    alphas = np.array(alphas)
    betas = np.array(betas)
    MADs = np.array(MADs)
    SDs = np.array(SDs)
    positions = np.array(positions)
    pixels = np.array(pixels)

if False:
    betas = np.arange(37.62, 37.68, .001)
    MADs = []
    foundbeta = []
    for beta in betas:
        (MAD, newbeta) = fit_spec(data, 1.000116, beta)

        print "----"
        print beta, MAD

        
        MADs.append(MAD)
        foundbeta.append(newbeta)

    MADs = np.array(MADs)
    foundbeta = np.array(foundbeta)
