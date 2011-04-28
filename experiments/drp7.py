'''
-------------------------------------------------------------------------------

Created on April 24 2011 by npk

File created to fix blended and saturated line problem found in drp6

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
    "/users/npk/desktop/c9/m110319_1949.fits")
band = 'H'

pl.ion()

DRAW = True
def fit_spec(data, pos, alpha, sinbeta, gamma, delta, band):
    global DRAW
    ar_h_lines = np.array([1.465435, 1.474317, 1.505062, 1.517684, 1.530607, 
        1.533334, 1.540685, 1.590403, 1.599386, 1.618444, 1.644107,
        1.674465, 1.744967, 1.791961])

    Ar_K_lines = np.array([1.982291, 1.997118, 2.032256, 2.057443, 
        2.062186, 2.099184, 2.133871, 2.15409, 2.204558, 2.208321,
        2.313952, 2.385154])

    Ne_J_lines = np.array([1.149125, 1.16719, 1.17227, 1.194655, 1.202994, 
        1.211564, 1.214306, 1.234677, 1.240622, 1.244273, 1.249108,
        1.27369, 1.280624, 1.296021, 1.301182, 1.321761, 1.327627,
        1.331685, 1.337077, 1.341026, 1.350788, 1.362638])

    Ne_J_lines = np.array([1.149125, 1.16719, 1.117227])

    AR_Y_lines = np.array([0.966043, 0.978718, 1.005481, 1.04809,
        1.067649, 1.07039, 1.088394, 1.110819])

    # NeAr H (1961, 1949)
    lines = np.array([1.465435, 1.474317, 1.505062, 1.517684, 1.530607, 
        1.533334, 1.540685, 1.590403, 1.599386, 1.618444, 1.644107,
        1.674465, 1.744967, 1.791961, \
        1.493386, 1.499041, 1.523487, 1.5352384, 1.5411803, 1.5608478,
        1.6027147, 1.6272797, 1.6409737, 1.6479254, 1.6793378,
        1.7166622])

    lines = np.array([ 1.540685, 1.5411803])


    lines.sort()

    bmap = {"Y": 6, "J": 5, "H": 4, "K": 3}
    order = bmap[band]
    d = 1e3/110.5
    pix = np.arange(2048.)

    ll = Fit.wavelength_model([order, pos, alpha, sinbeta, gamma, delta], pix)

    spec = np.median(data[pos-1:pos+1, :], axis=0)
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
        f = 0.9985
        roi = (f*lam < ll) & (ll < lam/f)

        if not roi.any(): 
            xs.append(0.0)
            sxs.append(np.inf)
            continue

        #lsf = Fit.mpfitpeak(pix[roi], spec[roi], 
         #   error=np.sqrt(np.abs(spec[roi])))


        lsf = Fit.mpfitpeak(pix[roi], spec[roi],
            error=np.sqrt(np.abs(spec[roi])))

        print lsf.params[2], lsf.params[1]
        print "LSF"

        pl.figure(6)
        pl.clf()
        pl.scatter(pix[roi], spec[roi])
        pl.plot(pix[roi], Fit.gaussian(lsf.params, pix[roi]))


        lsf = Fit.mpfitpeaks(pix[roi], spec[roi], 2,
            error=np.sqrt(np.abs(spec[roi])))

        print lsf.params[2], lsf.params[6], lsf.params[8]
        print "LSF"

        pl.figure(5)
        pl.clf()
        pl.scatter(pix[roi], spec[roi])
        pl.plot(pix[roi], Fit.multi_gaussian(lsf.params, pix[roi]))
    
        if lsf.perror is None:
            xs.append(0.0)
            sxs.append(np.inf)
            continue

        if lsf.status < 0:
            xs.append(0.0)
            sxs.append(np.inf)
            continue

        mnpix = np.min(pix[roi])
        mxpix = np.max(pix[roi])

        if (mnpix + 4) > lsf.params[1] < (mxpix-4):
            xs.append(0.0)
            sxs.append(np.inf)
            continue

        if mnpix < 7:
            xs.append(0.0)
            sxs.append(np.inf)
            continue

        if mxpix > 2040:
            xs.append(0.)
            sxs.append(np.inf)
            continue

        #if DRAW:
        if DRAW:
            pl.plot(pix[roi], spec[roi])
            pl.plot(pix[roi], spec[roi], 'o')
            pl.plot(pix[roi], Fit.gaussian(lsf.params, pix[roi]), 'r--')

        xval = lsf.params[1]
        xvals = lsf.perror[1]
        xs.append(xval)
        sxs.append(xvals)

        if DRAW: 
            if np.isfinite(xvals):
                pl.axvline(xval)

    if DRAW:
        pl.draw()

    def dofit(xs, slambdas, alpha, sinbeta, gamma, delta, band, y):
        xs = np.array(xs)
        ok = np.isfinite(slambdas)
        lsf = Fit.do_fit_wavelengths(xs[ok], lines[ok], alpha, sinbeta, 
                gamma, delta, band, y, error=slambdas[ok])
        (order, pixel_y, alpha, sinbeta, gamma, delta) = lsf.params

        print "order alpha    sinbeta   gamma  delta"
        print "%1.0i %5.7f %3.5f %3.2e %5.2f" % (order, alpha, sinbeta, 
                gamma, delta)

        ll = Fit.wavelength_model(lsf.params, pix)
        if DRAW:
            pl.figure(3)
            pl.clf()
            pl.plot(ll, spec)
            pl.title("Pixel %4.4f" % pos)
            for lam in lines:
                pl.axvline(lam, color='r')

            pl.draw()
    
        return [np.abs((
            Fit.wavelength_model(lsf.params, xs[ok]) - lines[ok]))*1e4, 
                lsf.params]

    slambdas = d/order * np.array(sxs)

    [deltas, params] = dofit(xs, slambdas, alpha, sinbeta, gamma, delta,
            band)

    if 0.4 < np.median(deltas) < 0.8:
        prev = np.median(deltas)
        [deltas, params] = dofit(xs, slambdas, params[1], params[2], 
                params[3], params[4], band)
        print "Previous MAD: %f is now %f" % (prev, np.median(deltas))

    return (params[1], params[2], params[3], params[4],
            np.median(deltas), np.std(deltas), pixel)


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
    for slit in range(1,47):
        os = "b%2.2ipos" % (slit*2-1)
        es = "b%2.2ipos" % (slit*2)
        csupos = (header[os] + header[es])/2.
        pixel = np.int(y0 - (slit-1) * 44.22)

        alpha = alpha_pixel(pixel)
        sinbeta = sinbeta_predict(csupos, pixel)
        gamma = gamma_pixel(pixel)
        delta = delta_pixel(pixel)


        print "------- Slit %i (%4.0i) @ %3.5f" % (slit, pixel, sinbeta)
        (newa, newb, newg, newd, MAD, SD, pixel) = \
            fit_spec(data, pixel, alpha, sinbeta, gamma, delta, band)


        1/0

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

                    (newa, newb, newg, newd, MAD, SD, pixel) = \
                            fit_spec(data, pixel, alpha+dalpha, 
                                    sinbeta+dsinbeta, gamma, delta, band)
                    if MAD < .2:
                        print "breaking at dalpha: %f, dsinbeta: %f" % (
                                dalpha, dsinbeta)
                        break
            
                if MAD < .2:
                    break


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
