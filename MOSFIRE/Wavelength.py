'''
==================================
MOSFIRE Wavelength Calibrations

This module is responsible for determining wavelength solutions.

The wavelength is represented as a physical model of MOSFIRE using the 
grating equation, the size of a pixel (18 micron), and the focal length
of the camera (250 mm), the lines per mm of the grating (110.5), and the
order of the grating (Y: 6, J: 5, H: 4, K: 3).

formally, a model is fit to a function:
 lambda(x, y | alpha, beta, gamma, delta) =

alpha d                                                                       3
------- cos(scale pixely)^-1 {sin(scale pixelx) - sin(beta)} + gamma (x-delta)
   m 

 where x and y are pixel values, alpha, beta, gamma, and delta are
INPUT:

OUTPUT:

npk April 26 2011

'''

import os

import numpy as np
import pylab as pl
import pyfits as pf

from MOSFIRE import CSU, Fit, IO, Options

__version__ = "27Apr2011"

reload(Options)
reload(CSU)
reload(IO)
reload(Fit)

def tester():
    handle_lambdas(['m110323_2718.fits'], 
            'npk_calib3_q1700_pa_0',
            Options.wavelength)
    
def handle_lambdas(imglist, maskname, options):
    '''
    handle_lambdas is the priamry entry point to the Wavelengths module.
    '''
    
    global bs
    path = os.path.join(options["outdir"], maskname)
    if not os.path.exists(path):
        raise Exception("Output directory '%s' does not exist. This " 
                "directory should exist." % path)

    for fname in imglist:
        fp = os.path.join(path, fname)

        mfits = IO.readfits_all(fp)
        fit_lambda(mfits, options)

def fit_lambda(mfits, options):
    global sol_2d
    np.seterr(all="ignore")
    
    (header, data, bs) = mfits
    linelist = pick_linelist(header)
    solutions = []
    lamout = np.zeros(shape=(2048, 2048))
    for slitno in range(1,47):
        print("-==== Fitting Slit %i" % slitno)
        parguess = guess_wavelength_solution(slitno, header, bs)
        sol_1d = fit_wavelength_solution(data, parguess, 
                linelist, options)

        sol_2d = fit_outwards(data, sol_1d, linelist, options)
        lamout = merge_solutions(lamout, slitno, parguess[4], bs, sol_2d, 
                options)

        sol = {"slitno": slitno, "center_sol": [sol_1d[1], sol_1d[2]], 
                "sigmas": sol_1d[3], "2d": sol_2d}
        solutions.append(sol)

        
        fn = "/Users/npk/desktop/c9_reduce/npk_calib3_q1700_pa_0/lam.fits"
        try: os.remove(fn)
        except: pass
        hdu = pf.PrimaryHDU(lamout)
        hdu.writeto(fn)
    return solutions


def guess_wavelength_solution(slitno, header, bs):
    '''Given a slit number guess the coefficient values
    return [order, y0, alpha, sinbeta, gamma, delta]
    '''

    band = header['filter'].rstrip()
    bmap = {"Y": 6, "J": 5, "H": 4, "K": 3}
    order = bmap[band]

    y0 = bs.csu_slit_to_pixel(slitno)
    csupos_mm = bs.csu_slit_center(slitno)

    print("Slit %i at %4.3f" % (slitno, csupos_mm))

    # The following values are determined through experimentation
    # with c9 data

    alpha_pixel = np.poly1d([1.257e-9, -2.440e-6, 9.965e-1])

    if band == 'Y' or band == 'J':
        sinbeta_position = np.poly1d([0.0239, 36.2])
        sinbeta_pixel = np.poly1d([-2.578e-7, 0.00054, -0.2365])
        gamma_pixel = np.poly1d([-1.19e-19, 2.442e-16, 4.181e-13])
    elif band == 'H' or band == 'K':
        sinbeta_position = np.poly1d([2.331e-2, 38.24])
        sinbeta_pixel = np.poly1d([-2.664e-7, 5.534e-4, -1.992e-1])
        gamma_pixel = np.poly1d([-1.630e-19, 3.461e-16, 5.807e-13])

    delta_pixel = np.poly1d([4.284e-5, -0.1145, 1219])


    return [alpha_pixel(y0),
            sinbeta_position(csupos_mm) + sinbeta_pixel(y0),
            gamma_pixel(y0),
            delta_pixel(y0),
            order,
            y0]

def find_known_lines(lines, ll, spec, options):
    inf = np.inf
    xs = []
    sxs = []
    sigmas = []

    pl.figure(4)
    pl.clf()
    pix = np.arange(2048.)
    for lam in lines:
        f = options["fractional-wavelength-search"]
        roi = (f*lam < ll) & (ll < lam/f)

        if not roi.any():
            xs.append(0.0)
            sxs.append(inf)
            continue

        
        lsf = Fit.mpfitpeak(pix[roi], spec[roi],
                error=1/np.sqrt(np.abs(spec[roi])))

        if (lsf.perror is None) or (lsf.status < 0):
            xs.append(0.0)
            sxs.append(inf)
            continue


        mnpix = np.min(pix[roi])
        mxpix = np.max(pix[roi])

        if (mnpix + 4) > lsf.params[1] < (mxpix-4):
            xs.append(0.)
            sxs.append(inf)
            continue


        if mnpix < 7:
            xs.append(0.0)
            sxs.append(inf)
            continue

        if mxpix > 2040:
            xs.append(0.0)
            sxs.append(inf)
            continue


        xs.append(lsf.params[1])
        sxs.append(lsf.perror[1])
        sigmas.append(lsf.params[2])

    return map(np.array, [xs, sxs, sigmas])

def fit_model_to_lines(xs, sxs, lines, parguess, options):

    ok = np.isfinite(sxs)

    if len(np.where(ok)[0]) < 3:
        return [[np.inf], parguess, None]

    slambda = sxs * dlambda_model(parguess)

    parinfo = [
        {'fixed': 0, 'value': parguess[0], 'parname': 'alpha', 'step': 1e-5,
            'limited': [0,0], 'limits': [0,0]},
        {'fixed': 0, 'value': parguess[1], 'parname': 'sinbeta', 
            'step': 1e-5, 'limited': [0,0], 'limits': [0, 0]},
        {'fixed': 0, 'value': parguess[2], 'parname': 'gamma','step': 1e-15,
            'limited': [1,1], 'limits': [0,20e-13]},
        {'fixed': 0, 'value': parguess[3], 'parname': 'delta', 'step': 1e-1,
            'limited': [1,1], 'limits': [0,2048]},
        {'fixed': 1, 'value': parguess[4], 'parname': 'order', 
            'limited': [0,0], 'limits': [3,7]},
        {'fixed': 1, 'value': parguess[5], 'parname': 'Y',
            'limited': [0,0], 'limits': [0,2048]}
    ]

    merit_function = Fit.mpfit_residuals(wavelength_model)
    lsf = Fit.mpfit_do(merit_function, xs[ok], lines[ok], 
            parinfo, error=slambda[ok])

    return [ np.abs((wavelength_model(lsf.params, xs[ok]) - lines[ok]))*1e4,
            lsf.params, lsf.perror]



def fit_wavelength_solution(data, parguess, lines, options, search_num=45):
    '''Tweaks the guessed parameter values and provides 1d lambda solution
    
    '''

    pix = np.arange(2048.)

    MAD = np.inf

    y0 = parguess[5]
    spec = np.median(data[y0-1:y0+1, :], 
        axis=0) # axis = 0 is spatial direction

    d = search_num*0.008
    dsinbetas = np.sort(np.abs(np.linspace(-d/2., d/2., search_num)))

    sinbetadirection = 1.0
    for dsinbeta in dsinbetas:
        dsinbeta *= sinbetadirection
        sinbetadirection *= -1

        pars = parguess
        pars[1] = parguess[1] + dsinbeta
        ll = wavelength_model(parguess, pix)
        [xs, sxs, sigmas] = find_known_lines(lines, ll, spec, options)
        [deltas, params, perror] = fit_model_to_lines(xs, sxs, lines, 
                pars, options)

        MAD = np.median(deltas)
        print("MAD: %3.3f A" % MAD)

        if MAD > 0.2: print "  search"
        else: break

    if MAD < 0.2:
        print("%3.5f %4.3f %3.3e %4.1f" % (params[0], params[1], params[2], 
            params[3]))
        return [deltas, params, perror, sigmas]
    else:
        print("Could not find parameters")
        return [[], parguess, None, []]

def fit_outwards(data, sol_1d, lines, options):

    def sweep(DY, N):
        deltas, params, perror, sigmas = sol_1d
        ret = []
        for deltay in range(N):
            params[5] += DY
            if (params[5] < 0) or (params[5] > 2047):
                print "Skipping out of bounds %i" % params[5]
                continue

            print("Fitting at %i" % params[5])
            [deltas, params, perror, sigmas] = fit_wavelength_solution(data, 
                    params, lines, options, search_num=10)
            success = True

            if (len(deltas) < 2) or (np.median(deltas) > .4):
                success = False

            ret.append([params, perror, np.median(deltas), success])

        return ret
    
    pix = np.arange(2048.)

    params_up = sweep(1, 20) # sweep up
    params_down = sweep(-1, 20) # sweep down

    params_down.reverse()
    params_down.extend(params_up)

    return params_down


def merge_solutions(lamout, slitno, order, bs, sol_2d, options):

    ff = filter(lambda x:
            (x[1] is not None) and
            (x[1][0] < 2e-5) and
            (x[2] < .1) and
            (x[3] == True), sol_2d)
    ar = np.array(map(lambda x: x[0], ff))

    pixel = ar[:,5]

    alpha = ar[:,0]
    beta = ar[:,1]
    gamma = ar[:,2]
    delta = ar[:,3]

    alphamodel = np.poly1d(np.polyfit(pixel, alpha, 1))
    betamodel  = np.poly1d(np.polyfit(pixel, beta, 1))
    gammamodel = np.poly1d(np.polyfit(pixel, gamma, 1))
    deltamodel = np.poly1d(np.polyfit(pixel, delta, 1))
    

    mn = np.min(pixel)-4
    if mn < 0: mn = 0
    mx = np.max(pixel)+4
    if mx > 2047: mx = 2047

    pixs = np.arange(2048.)
    for y in np.arange(np.int(mn), np.int(mx)+1):
        pars = [alphamodel(y),
                betamodel(y),
                gammamodel(y),
                deltamodel(y),
                order,
                y]
        ll = wavelength_model(pars, pixs)
        lamout[y,:] = ll


    return lamout

    
def dlambda_model(p):
    ''' Returns an approximate dlambda/dpixel '''
    x = 1024
    order = p[4]
    y = p[5]
    (alpha, sinbeta, gamma, delta) = p[0:4]
    sinbeta = np.radians(sinbeta)
    d = 1e3/110.5 # Groove spacing in micron
    pixelsize, focal_length = 18.0, 250e3 # micron
    scale = pixelsize/focal_length

    costerm = np.cos(scale * (y-1024))

    return scale/(order/d) * sinbeta / costerm


def wavelength_model(p, x):
    ''' Returns wavelength as function of pixel (x)
    
    The parameter list, p, contains
    p[0:3] -- alpha, beta, gamma, delta, model parameters
    p[4] -- the grating order.
    p[5] -- the pixel y position on detector [pix]

    x -- the x pixel position (dispersion direction)

    returns wavelength [micron]
    '''
    order = p[4]
    y = p[5]
    (alpha, sinbeta, gamma, delta) = p[0:4]
    sinbeta = np.radians(sinbeta)
    d = 1e3/110.5 # Groove spacing in micron
    pixelsize, focal_length = 18.0, 250e3 # micron
    scale = pixelsize/focal_length

    costerm = np.cos(scale * (y-1024))

    return alpha/(order/d) * 1/costerm * \
            (np.sin(scale * (x-1024)) + sinbeta) + \
            gamma * (x - delta)**3


def pick_linelist(header):
    band = header["filter"]
    Argon_on = header["pwstata8"] == 1
    Neon_on = header["pwstata7"] == 1

    assert(header["pwloca7"].rstrip() == "Neon")
    assert(header["pwloca8"].rstrip() == "Argon")

    lines = []
    if band == 'H':
        if Argon_on:
            lines.extend([1.465435, 1.474317, 1.505062, 1.517684, 1.530607, 
        1.533334, 1.540685, 1.590403, 1.599386, 1.618444, 1.644107,
        1.674465, 1.744967, 1.791961])
        if Neon_on:
            lines.extend([1.493386, 1.499041, 1.523487, 1.5352384, 
                    1.5411803, 1.5608478, 1.6027147, 1.6272797, 
                    1.6409737, 1.6479254, 1.6793378, 1.7166622])


    lines = np.array(lines)
    lines = np.sort(lines)

    return lines

    




if __name__ == "__main__":
    np.set_printoptions(precision=3)
    global bs, sol_2d
    pl.ion()
    tester()

    ff = filter(lambda x: (x[1] is not None) and (x[1][0] < 3e-5) and 
            (x[2] < .1) and (x[3] == True), sol_2d)

    ar = np.array(map(lambda x: x[0], ff))



