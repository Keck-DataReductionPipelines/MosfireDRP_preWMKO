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
npk   May  4 2011

'''

import os

import numpy as np
import pylab as pl
import pyfits as pf
from scipy.interpolate import interp1d

from MOSFIRE import CSU, Fit, IO, Options

__version__ = "27Apr2011"

reload(Options)
reload(CSU)
reload(IO)
reload(Fit)

#
# Glue code
#

    
def handle_lambdas(imglist, maskname, options):
    ''' handle_lambdas is the entry point to the Wavelengths module. '''
    
    path = os.path.join(options["outdir"], maskname)
    if not os.path.exists(path):
        raise Exception("Output directory '%s' does not exist. This " 
                "directory should exist." % path)

    for fname in imglist:
        fp = os.path.join(path, fname)

        mfits = IO.readmosfits(fp)
        fit_lambda(mfits, fname, maskname, options)
        apply_lambda(mfits, fname, maskname, options)

def fit_lambda(mfits, fname, maskname, options):
    '''Fit the two-dimensional wavelength solution to each science slit'''
    np.seterr(all="ignore")
    
    (header, data, bs) = mfits
    linelist = pick_linelist(header)
    solutions = []
    lamout = np.zeros(shape=(2048, 2048), dtype=np.float32)
    for slitno in range(1, 47):
        print("-==== Fitting Slit %i" % slitno)
        parguess = guess_wavelength_solution(slitno, header, bs)
        sol_1d = fit_wavelength_solution(data, parguess, 
                linelist, options)

        sol_2d = fit_outwards_xcor(data, sol_1d, linelist, options)
        lamout = merge_solutions(lamout, slitno, parguess[4], bs, sol_2d, 
                options)

        sol = {"slitno": slitno, "center_sol": [sol_1d[1], sol_1d[2]], 
                "sigmas": sol_1d[3], "2d": sol_2d, "lines": linelist,
                "csupos_mm": parguess[-1]}
        solutions.append(sol)

        path = os.path.join(options["outdir"], maskname)
        fn = os.path.join(path, "lambda_coeffs_{0}.npy".format(
            fname.replace(".fits","")))

        try: os.remove(fn)
        except: pass
        np.save(fn, solutions)

    return solutions

def apply_lambda(mfits, fname, maskname, options):
    '''Convert solutions into final output products'''

    (header, data, bs) = mfits

    band = header['filter'].rstrip()
    path = os.path.join(options["outdir"], maskname)
    fn = os.path.join(path, "slit-edges_{0}.npy".format(band))
    edgedata = np.load(fn)

    fn = os.path.join(path, "lambda_coeffs_{0}.npy".format(
        fname.replace(".fits","")))
    lambdadata = np.load(fn)
    
    slitedges = edgedata[0:-1]
    edgeinfo = edgedata[-1]

    # pixel_ and beta_set are lists containing the
    # measured values of beta for each scientific slit
    (pixel_set, alpha_set, beta_set, gamma_set, 
            delta_set) = mechanical_to_science_slit(bs, slitedges, lambdadata)

    print("{0}: Fitting across the mask...".format(maskname))
    print("{0}:\t...alpha ".format(maskname))
    alpha_lsf = fit_mask([0.99, 0, 0, 1024], pixel_set, 
            np.concatenate(alpha_set))
    print("{0}:\t...beta".format(maskname))
    beta_lsf = fit_mask([40., 0, 0, 1024], pixel_set, 
            np.concatenate(beta_set))
    print("{0}:\t...gamma".format(maskname))
    gamma_lsf = fit_mask([7e-13, 0, 0, 1024], pixel_set, 
            np.concatenate(gamma_set))
    print("{0}\t...delta".format(maskname))
    delta_lsf = fit_mask([1000, 0, 0, 1024.], pixel_set, 
            np.concatenate(delta_set))

    pixels = np.concatenate(pixel_set)

    fn = os.path.join(path, "mask_solution_{0}.npy".format(fname))
    try: os.remove(fn)
    except: pass
    np.save(fn, {"alpha_lsf": alpha_lsf, "beta_lsf": beta_lsf, 
        "gamma_lsf": gamma_lsf, "delta_lsf": delta_lsf,
        "pixels": pixel_set, "alphas": alpha_set, "betas": beta_set,
        "gammas": gamma_set, "deltas": delta_set})

    alphas, betas, gammas, deltas = map(np.concatenate, [alpha_set, beta_set,
        gamma_set, delta_set])
   
    def convfun(params, px, i):
        pp = np.copy(params[0:5])
        pp[4] = params[4+i]
        return mask_model(pp, [px])

    # write lambda
    lams = np.zeros((2048, 2048), dtype=np.float32)
    xx = np.arange(2048)
    for i in xrange(len(pixel_set)):
        edges = slitedges[i]
        top = edges["top"](xx)
        bottom = edges["bottom"](xx)
        px = np.arange(np.min(bottom), np.max(top))

        alphas = convfun(alpha_lsf.params, px, i)
        betas = convfun(beta_lsf.params, px, i)
        gammas = convfun(gamma_lsf.params, px, i)
        deltas = convfun(delta_lsf.params, px, i)

        cnt = 0
        for j in px:
            lams[j,:] = wavelength_model(
                    [alphas[cnt],
                        betas[cnt],
                        gammas[cnt],
                        deltas[cnt],
                        4,
                        j], xx)
            cnt += 1



    print("{0}: writing lambda".format(maskname))
    hdu = pf.PrimaryHDU(lams*1e4)
    fn = os.path.join(path, "lambda_solution_{0}.fits".format(fname))
    try: os.remove(fn)
    except: pass
    hdu.writeto(fn)

    print("{0}: rectifying".format(maskname))
    rectified = np.zeros((2048, 2048), dtype=np.float32)
    ll_fid = lams[1024, :]
    for i in xrange(2048):
        ll = lams[i,:]
        ss = data[i,:]

        f = interp1d(ll, ss, bounds_error=False)
        rectified[i,:] = f(ll_fid)

    hdu = pf.PrimaryHDU(rectified)
    fn = os.path.join(path, "rectified_{0}.fits".format(fname))
    try: os.remove(fn)
    except: pass
    hdu.writeto(fn)


    pl.figure(4)
    pl.clf()

    fvs = []
    poss= [2012, 1955, 1925, 1868, 1824, 1774, 1740, 1700, 1654, 1614, 1575,
            1519, 1476, 1433, 1396, 1300, 1073, 1036, 904, 858, 806, 725,
            682, 633, 591, 546, 497, 459, 408, 365, 316, 282, 235, 192]
    for i in poss:
        ll = lams[i,:]
        ss = data[i,:]

        roi = np.abs(ll-1.5997)<.001

        if not roi.any(): continue
        pl.plot(ll[roi]*1e4, ss[roi]/np.max(ss[roi]))

        pp = Fit.mpfitpeak(ll[roi]*1e4, ss[roi])
        fvs.append(pp.params[1])

    pl.ylim([-.1,1.1])
    pl.figure(5)
    pl.clf()
    pl.plot(poss, fvs,'*')

def mechanical_to_science_slit(bs, slitedges, lambdadata):
    '''Convert mechanical slit fits to fits accross the science slit'''
    all_pixs = []
    all_alphas = []
    all_betas = []
    all_gammas = []
    all_deltas = []

    for i in xrange(len(bs.ssl)):
        ss = bs.ssl[i]
        edges = slitedges[i]

        csuslits = bs.scislit_to_csuslit(i)

        scislit_pixs = []
        scislit_alphas = []
        scislit_betas = []
        scislit_gammas = []
        scislit_deltas = []
        for csuslit in csuslits:
            sol = lambdadata[csuslit-1]
            ff = filter(filter_2d_solutions, sol["2d"])
            assert(sol["slitno"] == csuslit)
            ar = np.array(map(lambda x:x[0], ff))

            pix = ar[:, 5]
            alpha = ar[:, 0]
            beta = ar[:, 1]
            gamma = ar[:, 2]
            delta = ar[:, 3]

            scislit_pixs.extend(pix)
            scislit_alphas.extend(alpha)
            scislit_betas.extend(beta)
            scislit_gammas.extend(gamma)
            scislit_deltas.extend(delta)

        (scislit_pixs, scislit_alphas, scislit_betas, scislit_gammas, 
                scislit_deltas) = map(np.array, [scislit_pixs, 
                    scislit_alphas, scislit_betas, scislit_gammas, 
                    scislit_deltas])

        srt = np.argsort(scislit_pixs)
        scislit_pixs = scislit_pixs[srt]
        scislit_alphas = scislit_alphas[srt]
        scislit_betas = scislit_betas[srt]
        scislit_gammas = scislit_gammas[srt]
        scislit_deltas = scislit_deltas[srt]

        all_pixs.append(scislit_pixs)
        all_alphas.append(scislit_alphas)
        all_betas.append(scislit_betas)
        all_gammas.append(scislit_gammas)
        all_deltas.append(scislit_deltas)

    return [all_pixs, all_alphas, all_betas, all_gammas, all_deltas]

def filter_2d_solutions(x):
    '''Select quality fits in two dimensional solution'''
    sds = x[1] # standard deviation of first fit parameter
    MAD = x[2] # Median absolute deviation of fit to data
    success = x[3] # Binary success criteria

    return (sds is not None) \
        and (sds[0] < 2e-5) \
        and (MAD < 0.1) \
        and (success)

def tester():
    handle_lambdas(['m110323_2718.fits'], 
            'npk_calib3_q1700_pa_0',
            Options.wavelength)

def apply_lambda_tester():
    fn = "/users/npk/desktop/c9_reduce/npk_calib3_Q1700_pa_0/m110323_2718.fits"
    mfits = IO.readmosfits(fn)

    apply_lambda(mfits, "m110323_2718", "npk_calib3_Q1700_pa_0", 
            Options.wavelength)

#
# Fitting Methods
#   

# Physical models for instrument
def param_guess_functions(band):
    '''Parameters determined from experimentation with cooldown 9 data'''
    alpha_pixel = np.poly1d([-8.412e-16, 3.507e-12, -3.593e-9, 
        6.303e-9, 0.9963])

    if band == 'Y' or band == 'J':
        sinbeta_position = np.poly1d([0.0239, 36.2])
        sinbeta_pixel = np.poly1d([-2.578e-7, 0.00054, -0.2365])
        gamma_pixel = np.poly1d([1.023e-25, -4.313e-22, 7.668e-17, 6.48e-13])
    elif band == 'H' or band == 'K':
        sinbeta_position = np.poly1d([2.331e-2, 38.24])
        sinbeta_pixel = np.poly1d([-2.664e-7, 5.534e-4, -1.992e-1])
        gamma_pixel = np.poly1d([1.033e-25, -4.36e-22, 4.902e-19, -8.021e-17,
            6.654e-13])

    delta_pixel = np.poly1d([-1.462e-11, 6.186e-8, -5.152e-5, -0.0396, 1193])

    return [alpha_pixel, sinbeta_position, sinbeta_pixel, 
            gamma_pixel, delta_pixel]
    
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

    [alpha_pixel, sinbeta_position, sinbeta_pixel, gamma_pixel, 
            delta_pixel] = param_guess_functions(band)

    return [alpha_pixel(y0),
            sinbeta_position(csupos_mm) + sinbeta_pixel(y0),
            gamma_pixel(y0),
            delta_pixel(y0),
            order,
            y0, 
            csupos_mm]

def find_known_lines(lines, ll, spec, options):
    inf = np.inf
    xs = []
    sxs = []
    sigmas = []

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

def fit_model_to_lines(xs, sxs, lines, parguess, options, fixed):

    ok = np.isfinite(sxs)

    if len(np.where(ok)[0]) < 3:
        return [[np.inf], parguess, None]

    slambda = sxs * dlambda_model(parguess)

    parinfo = [
        {'fixed': 0, 'value': parguess[0], 'parname': 'alpha', 'step': 1e-5,
            'limited': [0,0], 'limits': [0,0]},
        {'fixed': 0, 'value': parguess[1], 'parname': 'sinbeta', 
            'step': 1e-7, 'limited': [0,0], 'limits': [0, 0]},
        {'fixed': fixed, 'value': parguess[2], 'parname': 'gamma','step': 1e-13,
            'limited': [1,1], 'limits': [0,20e-13]},
        {'fixed': fixed, 'value': parguess[3], 'parname': 'delta', 'step': 1e-1,
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

def fit_wavelength_solution(data, parguess, lines, options, search_num=45,
        fixed=False):
    '''Tweaks the guessed parameter values and provides 1d lambda solution
    
    '''

    pix = np.arange(2048.)

    MAD = np.inf

    y0 = parguess[5]
    spec = np.median(data[y0-1:y0+1, :], 
        axis=0) # axis = 0 is spatial direction

    d = search_num*0.003
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
                pars, options, fixed)

        MAD = np.median(deltas)

        if MAD > 0.1: 
            pass # continue to search
        else: 
            break


    if MAD < 0.1:
        print("%3.5f %4.3f %3.3e %4.1f" % (params[0], params[1], params[2], 
            params[3]))
        return [deltas, params, perror, sigmas]
    else:
        print("Could not find parameters")
        return [[], parguess, None, []]

#
# Two dimensional wavelength fitting
#
def fit_outwards_xcor(data, sol_1d, lines, options):
    lags = np.arange(-5,5)
    pix = np.arange(2048.)

    print("outwards fitting...")
    def sweep(deltays):
        deltas, params0, perror, sigmas = sol_1d
        params0 = np.array(params0)
        y0 = params0[5]
        ret = []
        spec = np.median(data[y0-1:y0+1, :], axis=0)

        for deltay in deltays:
            params = params0.copy()
            params[5] = y0+deltay

            if (params[5] < 0) or (params[5] > 2047):
                print("%i: Skipping out of bounds %i" % (deltay, params[5]))
                continue

            spec2 = np.median(data[params[5]-1:params[5]+1, :], axis=0)
            xcs = []
            for lag in lags:
                xc = np.sum(spec * np.roll(spec2, lag))
                xcs.append(xc)

            fp = Fit.mpfitpeak(lags, np.array(xcs))
            spec = spec2

            params[1] -= np.degrees(fp.params[1] * 18/250e3)

            [deltas, params, perror, sigmas] = fit_wavelength_solution(
                    data, params, lines, options, search_num=5, fixed=False)

            success = True
            if (len(deltas) < 2) or (np.median(deltas) > .4):
                success = False

            print "%i: Success: %i" % (params[5], success)
            ret.append([params, perror, np.median(deltas), success])

        return ret

    pix = np.arange(2048.)

    params_up = sweep(range(0,20,2))
    params_down = sweep(range(-1,-20,-2))
    params = params_down

    params.reverse()
    params.extend(params_up)

    return params

def merge_solutions(lamout, slitno, order, bs, sol_2d, options):
    ff = filter(lambda x:
            (x[1] is not None) and
            (x[1][0] < 2e-5) and
            (x[2] < .1) and
            (x[3] == True), sol_2d)
    ar = np.array(map(lambda x: x[0], ff))

    pixels = ar[:,5]

    alphas = ar[:,0]
    betas = ar[:,1]
    gammas = ar[:,2]
    deltas = ar[:,3]

    alphamodel = np.poly1d(np.polyfit(pixels, alphas, 2))
    betamodel  = np.poly1d(np.polyfit(pixels, betas, 2))
    gammamodel = np.poly1d(np.polyfit(pixels, gammas, 2))
    deltamodel = np.poly1d(np.polyfit(pixels, deltas, 2))

    print "Alpha scatter: %3.6f" % np.std(alphas-alphamodel(pixels))
    print " Beta scatter: %3.4f" % np.std(betas-betamodel(pixels))
    print "Fitting pixels: ", pixels

    print betas
    print alphas

    mn = np.min(pixels)-4
    if mn < 0: mn = 0
    mx = np.max(pixels)+4
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


def fit_mask(pars, pixel_set, ys):
    '''Fit mask model function to the whole mask'''
    merit_fun = Fit.mpfit_residuals(mask_model)
    pars.extend(np.zeros(len(pixel_set)))
    pars = np.array(pars)

    parinfo = []
    for par in pars:
        parinfo.append({"fixed": 0, "value": par, "limited": [0, 0], 
            "limits": [0, 0]})

    return Fit.mpfit_do(merit_fun, pixel_set, ys, parinfo)

# Model Functions

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

def mask_model(p, xs):
    '''Fit a continuous smooth function to parameters in the mask.

    parameters:
        linear model is:
        x = xs - p[3]         2         3
        p[0] + p[1] x + p[2] x  + p[3] x  + discontinuity

        p[4:] -- [N] list of discontinuities
        '''

    cpix = p[0]
    cpar = p[1]
    radius_pix = p[2]
    radius_par = p[3]
    coeffs = p[4:]

    vals = []
    for i in xrange(len(coeffs)):
        x = np.array(xs[i]) - p[3]
        c = coeffs[i]
        y = p[0] + p[1] * x + p[2] * x*x + c
        vals.extend(y)

    return np.array(vals).ravel()

 
def plot_mask_fits(fname):
    def plotfun(px, ys, lsf, i):
        pl.plot(px, ys, '.')
        params = np.copy(lsf.params[0:5])
        params[4] = lsf.params[4+i]
        my = mask_model(params, [px])
        pl.plot(px, my)

    pl.figure(3)
    pl.clf()
    pl.subplot(2,2,1)
    pl.plot(pixels, alphas, '.')
    for i in xrange(len(pixel_set)):
        plotfun(pixel_set[i], alpha_set[i], alpha_lsf, i)

    pl.subplot(2,2,2)
    for i in xrange(len(pixel_set)):
        plotfun(pixel_set[i], beta_set[i], beta_lsf, i)

    pl.subplot(2,2,3)
    pl.plot(pixels, np.concatenate(gamma_set), '.')
    for i in xrange(len(pixel_set)):
        plotfun(pixel_set[i], gamma_set[i], gamma_lsf, i)

    pl.subplot(2,2,4)
    pl.plot(pixels, np.concatenate(delta_set), '.')
    for i in xrange(len(pixel_set)):
        plotfun(pixel_set[i], delta_set[i], delta_lsf, i)



def plot_data_quality(fname):
    global solution
    from matplotlib.backends.backend_pdf import PdfPages


    solutions = np.load(fname)

    pp = PdfPages("/users/npk/desktop/mp.pdf")
    filter_fun = (lambda x:
            (x[1] is not None) and
            (x[1][0] < 2e-5) and
            (x[2] < .1) and
            (x[3] == True))

    all_pix = []
    all_alphas = []
    all_betas = []
    all_gammas = []
    all_deltas = []
    for solution in solutions:
        sol_2d = solution["2d"]
        print "Slit: {0}".format(solution["slitno"])
        ff = filter(filter_fun, sol_2d)
        ar = np.array(map(lambda x: x[0], ff))

        pixels = ar[:,5]

        alphas = ar[:,0]
        betas = ar[:,1]
        gammas = ar[:,2]
        deltas = ar[:,3]
                

        all_pix.extend(pixels)
        all_alphas.extend(alphas)
        all_betas.extend(betas)
        all_gammas.extend(gammas)
        all_deltas.extend(deltas)

        alphamodel = np.poly1d(np.polyfit(pixels, alphas, 1))
        betamodel  = np.poly1d(np.polyfit(pixels, betas, 1))
        gammamodel = np.poly1d(np.polyfit(pixels, gammas, 1))
        deltamodel = np.poly1d(np.polyfit(pixels, deltas, 1))

        print "Scatters: {0:3.5} {1:3.5} {2:3.5} {3:3.5}".format(
                np.std(alphas-alphamodel(pixels)),
                np.std(betas-betamodel(pixels)),
                np.std(gammas-gammamodel(pixels)),
                np.std(deltas-deltamodel(pixels)),
                )

        pl.clf()
        pl.subplot(2,2,1)
        pl.title("Slit {0}".format(solution["slitno"]))
        pl.scatter(pixels, alphas)
        pl.plot(pixels, alphamodel(pixels))
        pl.ylabel(r'$\alpha$')

        pl.subplot(2,2,2)
        pl.scatter(pixels, betas)
        pl.plot(pixels, betamodel(pixels))
        pl.ylabel(r'$\beta$')

        pl.subplot(2,2,3)
        pl.scatter(pixels, gammas)
        pl.plot(pixels, gammamodel(pixels))
        pl.ylim([0,1e-12])
        pl.ylabel(r'$\gamma$')

        pl.subplot(2,2,4)
        pl.scatter(pixels, deltas)
        pl.plot(pixels, deltamodel(pixels))
        pl.ylabel(r'$\delta$')


        pp.savefig()

    [alpha_pixel, sinbeta_position, sinbeta_pixel, gamma_pixel, 
            delta_pixel] = param_guess_functions('H')
    pl.clf()
    pl.subplot(1,1,1)
    pl.scatter(all_pix, all_alphas, c=all_deltas)
    pl.plot(all_pix, alpha_pixel(all_pix), 'r')

    ff = np.poly1d(np.polyfit(all_pix, all_alphas, 4))
    pl.plot(all_pix, ff(all_pix))
    print "Alpha: ", ff
    pl.ylabel(r'$\alpha$')
    pp.savefig()

    pl.clf()
    delts = all_alphas - ff(all_pix)
    pl.scatter(all_pix, delts, c=all_gammas)
    pl.ylabel(r'$\Delta \alpha$')
    print "Scatter is {0} pixels".format(np.std(delts)*2048)
    pp.savefig()

    pl.clf()
    pl.scatter(all_pix, all_betas, s=.1)
    pl.ylabel(r'$\beta$')
    pp.savefig()

    pl.clf()
    pl.scatter(all_pix, all_gammas, c=all_gammas)
    pl.plot(all_pix, gamma_pixel(all_pix), 'r')
    ff = np.poly1d(np.polyfit(all_pix, all_gammas, 4))
    print "Gamma: ", ff

    pl.plot(all_pix, ff(all_pix), 'b')
    pl.ylabel(r'$\gamma$')
    pp.savefig()

    pl.clf()
    delta_pixel = np.poly1d([4.284e-5, -0.1145, 1219])
    pl.scatter(all_pix, all_deltas, c=all_gammas)
    pl.plot(all_pix, delta_pixel(all_pix), 'r')
    
    ff = np.poly1d(np.polyfit(all_pix, all_deltas, 4))
    print "Delta: ", ff

    pl.ylabel(r'$\delta$')
    pp.savefig()

    pp.close()

if __name__ == "__main__":
    np.set_printoptions(precision=3)
    global bs, sol_2d, solution
    if False:
        tester()

        ff = filter(lambda x: (x[1] is not None) and (x[1][0] < 3e-5) and 
                (x[2] < .1) and (x[3] == True), sol_2d)

        ar = np.array(map(lambda x: x[0], ff))
    else:
        plot_data_quality("/Users/npk/desktop/c9_reduce/npk_calib3_q1700_pa_0/lambda_coeffs_m110323_2718.npy")



