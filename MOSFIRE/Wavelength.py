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

 where x and y are pixel values, alpha, beta, gamma, and delta are measured
 parameters.

 scale is (pixel size) / (camera focal length)

-- Helper functions also exist for determining the on-order region of 
a spectrum --


INPUT:

OUTPUT:

npk Apr/May  2012 - Significant enhancements w/ first light data
npk April 26 2011
npk   May  4 2011

'''

import os
import time

import numpy as np
import pylab as pl
import pyfits as pf
from multiprocessing import Pool
from scipy.interpolate import interp1d
from matplotlib.widgets import Button

from MOSFIRE import CSU, Fit, IO, Options, Filters

#from IPython.Shell import IPShellEmbed
#ipshell = IPShellEmbed()

import pdb

__version__ = "1May2012"


MADLIMIT = 0.1
try:
    __IPYTHON__
    reload(Options)
    reload(CSU)
    reload(IO)
    reload(Fit)
    reload(Filters)
except:
    pass

#
# Glue code
#

    
def handle_lambdas(filelist, maskname, options):
    ''' handle_lambdas is the entry point to the Wavelengths module. '''
    
    path = os.path.join(options["outdir"], maskname)
    if not os.path.exists(path):
        raise Exception("Output directory '%s' does not exist. This " 
                "directory should exist." % path)

    for fname in filelist:
        fp = os.path.join(path, fname)

        mfits = IO.readmosfits(fp)
        fit_lambda(mfits, fname, maskname, options)
        apply_lambda(mfits, fname, maskname, options)

def fit_lambda(mfits, fname, maskname, options):
    '''Fit the two-dimensional wavelength solution to each science slit'''
    global bs, data, lamout, center_solutions
    np.seterr(all="ignore")
    
    fnum = fname.rstrip(".fits")
    
    print maskname
    path = os.path.join(options["outdir"], maskname)
    fn = os.path.join(path, "lambda_coeffs_{0}.npy".format(fnum))

    print "Writing to: ", fn

    (header, data, bs) = mfits
    band = header['filter'].rstrip()
    center_solutions = IO.load_lambdacenter(fnum, maskname, options)
    
    solutions = []
    lamout = np.zeros(shape=(2048, 2048), dtype=np.float32)

    tock = time.time()
    if True:
        p = Pool()
        solutions = p.map(fit_lambda_helper, range(1,47))
        p.close()


    tick = time.time()

    print "-----> Mask took %i" % (tick-tock)

    try: os.remove(fn)
    except: pass
    np.save(fn, solutions)

    return solutions


def fit_lambda_helper(slitno):
    '''This helper function exists for multiprocessing suport'''
    
    global bs, data, lamout, center_solutions

    slitidx = slitno-1

    tick = time.time()

    print("-==== Fitting Slit %i" % slitno)

    assert(center_solutions[slitidx]['slitno'] == slitno)
    sol_1d = center_solutions[slitidx]["sol_1d"]
    linelist = center_solutions[slitidx]["linelist"]


    sol_2d = fit_outwards_xcor(data, sol_1d, linelist, Options.wavelength,
            slitno)

    lamout = merge_solutions(lamout, slitno, sol_1d[1][4], bs, sol_2d, 
            Options.wavelength)


    sol = {"slitno": slitno, "center_sol": [sol_1d[1], sol_1d[2]], 
            "sigmas": sol_1d[3], "2d": sol_2d, "lines": linelist,
            "csupos_mm": sol_1d[-1]}
    print "%i] TOOK: %i" % (slitno, time.time()-tick)
    return sol


def fit_lambda_interactively(mfits, fname, maskname, options):
    '''Fit the two-dimensional wavelength solution to each science slit'''
    np.seterr(all="ignore")
    
    print maskname
    path = os.path.join(options["outdir"], maskname)
    fn = os.path.join(path, "lambda_center_coeffs_{0}.npy".format(
        fname.replace(".fits","")))

    print "Writing to: ", fn

    (header, data, bs) = mfits
    linelist = pick_linelist(header)
    
    try: solutions = np.load(fn)
    except: solutions = None

    lamout = np.zeros(shape=(2048, 2048), dtype=np.float32)

    tock = time.time()
    
    fig = pl.figure(1,figsize=(16,8))
    pl.ion()
    II = InteractiveSolution(fig, mfits, linelist, options, 1,
            solutions=solutions)
    pl.show()

    print len(II.solutions)

    print "save to: ", fn
    np.save(fn, np.array(II.solutions))

def apply_lambda(mfits, fname, maskname, options):
    '''Convert solutions into final output products'''

    (header, data, bs) = mfits

    band = header['filter'].rstrip()
    bmap = {"Y": 6, "J": 5, "H": 4, "K": 3}
    order = bmap[band]

    fnum = fname.rstrip(".fits")

    path = os.path.join(options["outdir"], maskname)
    edgedata = IO.load_edges(maskname, band, options)
    lambdadata = IO.load_lambdadata(fnum, maskname, band, options)
    print len(lambdadata)
    
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

    fn = os.path.join(path, "mask_solution_{0}.npy".format(
        fname.replace(".fits", "")))
    try: os.remove(fn)
    except: pass

    mask_fit_pars = [{"alpha_lsf": alpha_lsf, "beta_lsf": beta_lsf, 
        "gamma_lsf": gamma_lsf, "delta_lsf": delta_lsf,
        "pixels": pixel_set, "alphas": alpha_set, "betas": beta_set,
        "gammas": gamma_set, "deltas": delta_set}]
    np.save(fn, mask_fit_pars)

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
                        order,
                        j], xx)
            cnt += 1



    print("{0}: writing lambda".format(maskname))
    hdu = pf.PrimaryHDU(lams*1e4)
    fn = os.path.join(path, "lambda_solution_{0}".format(fname))
    try: os.remove(fn)
    except: pass
    hdu.writeto(fn)

    print("{0}: rectifying".format(maskname))
    dlam = np.median(np.diff(lams[1024,:]))
    hpp = Filters.hpp[band] 
    ll_fid = np.arange(hpp[0], hpp[1], dlam)
    nspec = len(ll_fid)

    rectified = np.zeros((2048, nspec), dtype=np.float32)

    for i in xrange(2048):
        ll = lams[i,:]
        ss = data[i,:]

        f = interp1d(ll, ss, bounds_error=False)
        rectified[i,:] = f(ll_fid)

    hdu = pf.PrimaryHDU(rectified)
    fn = os.path.join(path, "rectified_{0}".format(fname))
    try: os.remove(fn)
    except: pass
    hdu.writeto(fn)

    fvs = []
    poss= [2012, 1955, 1925, 1868, 1824, 1774, 1740, 1700, 1654, 1614, 1575,
            1519, 1476, 1433, 1396, 1300, 1073, 1036, 904, 858, 806, 725,
            682, 633, 591, 546, 497, 459, 408, 365, 316, 282, 235, 192]
    for i in poss:
        ll = lams[i,:]
        ss = data[i,:]

        roi = np.abs(ll-1.5997)<.001

        if not roi.any(): continue

        pp = Fit.mpfitpeak(ll[roi]*1e4, ss[roi])
        fvs.append(pp.params[1])

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
        
        print "SS#: %i, edges: %f" % (i, edges['top'](1024))

        csuslits = bs.scislit_to_csuslit(i)

        scislit_pixs = []
        scislit_alphas = []
        scislit_betas = []
        scislit_gammas = []
        scislit_deltas = []
        for csuslit in csuslits:
            sol = lambdadata[csuslit-1]
            ff = filter_2d_solutions(sol["2d"])

            assert(sol["slitno"] == csuslit)
            ar = np.array(map(lambda x:x[0], ff))

            try: 
                pix = ar[:, 5]
                alpha = ar[:, 0]
                beta = ar[:, 1]
                gamma = ar[:, 2]
                delta = ar[:, 3]
            except:
                #print "Skipping %i" % csuslit
                continue

            print alpha

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

def filter_2d_solutions(vec):
    '''Select quality fits in two dimensional solution'''


    def filter_fun(x):
        sds = x[1] # standard deviation of first fit parameter
        MAD = x[2] # Median absolute deviation of fit to data
        success = x[3] # Binary success criteria

        return (sds is not None) \
            and (sds[0] < 2e-5) \
            and (MAD < MADLIMIT) \
            and (success)

    return filter(filter_fun, vec)


#
# Fitting Methods
#   

# Physical models for instrument
def param_guess_functions(band):
    '''Parameters determined from experimentation with cooldown 9 data'''

    fudge_npk = 1.00100452
    fudge_npk = 1.0
    alpha_pixel = np.poly1d([-8.412e-16, 3.507e-12, -3.593e-9, 
        6.303e-9, 0.9963]) * fudge_npk

    # Note that these numbers were tweaked by hand by npk on 28 apr
    # they are not reliable. The fudge_* factors will need to change
    # or dissapear
    if band == 'Y' or band == 'J':
        fudge_npk = 0.008
        sinbeta_position = np.poly1d([0.0239, 36.2 + fudge_npk])
        sinbeta_pixel = np.poly1d([-2.578e-7, 0.00054, -0.2365])
        gamma_pixel = np.poly1d([1.023e-25, -4.313e-22, 7.668e-17, 6.48e-13])
    elif band == 'H' or band == 'K':
        fudge_npk = -0.05
        sinbeta_position = np.poly1d([2.331e-2, 38.24]) + fudge_npk
        sinbeta_pixel = np.poly1d([-2.664e-7, 5.534e-4, -1.992e-1])

        fudge_npk = 1
        gamma_pixel = np.poly1d([1.033e-25, -4.36e-22, 4.902e-19, -8.021e-17,
            6.654e-13]) * fudge_npk

    delta_pixel = np.poly1d([-1.462e-11, 6.186e-8, -5.152e-5, -0.0396,
        1193])  - 50

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
    
    lines = []

    if band == 'H':
        lines.extend([1.5056, 1.5833, 1.6692, 1.7653, 1.7880, 1.7994])
        return np.array(lines)
    if band == 'J':
        # From ccs oh_lines_j.txt
        lines = np.array([
            11538.7582 ,
            11591.7013 ,
            11627.8446 ,
            11650.7735 ,
            11696.3379 ,
            11716.2294 ,
            #11771.2773 ,
            11788.0779 ,
            11866.4924 ,
            11988.5382 ,
            12007.0419 ,
            12030.7863 ,
            12122.4957 ,
            12135.8356 ,
            12154.9582 ,
            12196.3557 ,
            12229.2777 ,
            12257.7632 ,
             12286.964 ,
            12325.9549 ,
            12351.5321 ,
            12400.8893 ,
             12423.349 ,
            12482.8503 ,
              12502.43 ,
            12589.2998 ,
            12782.9052 ,
            12834.5202 ,
            12905.5773 ,
            12921.1364 ,
            12943.1311 ,
            12985.5595 ,
            13021.6447 ,
             13052.818 ,
            13085.2604 ,
            13127.8037 ,
            13156.9911 ,
            13210.6977 ,
            13236.5414 ,
            13301.9624 ,
            13324.3509 ,
             13421.579])/1e4

    if band == 'K':
        # from ccs oh_lines_K.txt
        lines = np.array([
        19518.4784 ,
        19593.2626 ,
        19618.5719 ,
        19642.4493 ,
         19678.046 ,
        19701.6455 ,
        19736.4099 ,
        19751.3895 ,
        19771.9063 ,
        19839.7764 ,
        20008.0235 ,
        20193.1799 ,
        20275.9409 ,
         20339.697 ,
        20412.7192 ,
         20499.237 ,
        20563.6072 ,
         20729.032 ,
        20860.2122 ,
        20909.5976 ,
         #21033.062 ,
         #21115.889 ,
          #21156.24 ,
        21176.5323 ,
        #21232.4797 ,
        21249.5368 ,
        21279.1406 ,
        21507.1875 ,
        21537.4185 ,
        21580.5093 ,
        21711.1235 ,
        21802.2757 ,
         21873.507 ,
        21955.6857 ,
        22125.4484 ,
        22312.8204 ,
        22460.4183 ,
        22517.9267 ,
        22690.1765 ,
        22742.1907 ,
        22985.9156 ])/1e4

        #lines.extend([1.153879, 1.15917, #1.203083, 1.212249, 
            #1.222931, 1.23516,
            #1.242335, 1.268577, 1.284502, 1.29211, 1.268577, 1.278239, 1.29211,
            #1.302164, 1.308528, 1.323654, 
            #1.342156])

    lines = np.array(lines)
    return np.sort(lines)


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

def plot_mask_solution_ds9(fname, maskname, options):
    '''makes a ds9 region file guessing the wavelength solution'''


    inpath = options["indir"]
    path = os.path.join(options["outdir"], maskname)
    if not os.path.exists(path):
        raise Exception("Output directory '%s' does not exist. This " 
                "directory should exist." % path)

    fn = os.path.join(inpath, fname)
    (header, data, bs) = IO.readmosfits(fn)

    linelist = pick_linelist(header)

    ds9 = '''# Region file format: DS9 version 4.1
global color=red 
'''
    pix = np.arange(2048)
    colors = ["red", "blue"]
    cidx = 0


    for i in xrange(len(bs.ssl)):
        slits = bs.scislit_to_csuslit(i)

        print "Guessing: ", slits
        
        cidx = (cidx + 1) % 2
        color = colors[cidx]
        for slitno in slits:
            guess = guess_wavelength_solution(slitno, header, bs)
            ll = wavelength_model(guess, pix)

            if bs.is_alignment_slitno(slitno): color = 'green'
            
            for line in linelist:
                x = np.argmin(np.abs(ll - line))

                ds9 += "circle(%f, %f, 1) # color=%s text={}\n" % (x, guess[5],
                        color)

    path = Options.wavelength['outdir']


    fname = fname.rstrip(".fits")
    fn = os.path.join(path, maskname, ("guess_waves_%s.reg" % fname))
    try: os.remove(fn)
    except: pass

    try:
        f = open(fn, 'w')
        f.write(ds9)
        f.close()
    except:
        print "Could not write %s" % fn


def estimate_half_power_points(slitno, header, bs):
    ''' This helper function is used to determine the filter half-power points.
    This function is primarily used by the flat-field code to determine the 
    on order regions of an image.  '''

    band = header['filter'].rstrip()
    parguess = guess_wavelength_solution(slitno, header, bs)
    pix = np.arange(2048.)
    ll = wavelength_model(parguess, pix)


    hpp = Filters.hpp[band]
    return [ np.argmin(np.abs(ll-hpp[0])), np.argmin(np.abs(ll-hpp[1])) ]


def find_known_lines(lines, ll, spec, options):
    ''' 
    lines[N]: list of lines in wavelength units
    ll[2048]: lambda vector
    spec[2048]: spectrum vector (as function of lambda)
    options: wavelength options
    '''
    inf = np.inf
    xs = []
    sxs = []
    sigmas = []

    pix = np.arange(2048.)

    DRAW = False
    if DRAW:
        pl.ion()
        pl.clf()
        pl.figure(3, figsize=(23,4))
        pl.plot(spec)
        pl.xlim([0,500])
        pl.ylim([-50,2000])

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

        if DRAW:
            pl.axvline(lsf.params[1], color='red')
            pl.text(lsf.params[1], 1200, "%s,  %s" % (lam, lsf.params[1]),
                    fontsize=10,rotation=90)

    if DRAW:
        #pl.xlim([500, 550])
        pl.draw()

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
        {'fixed': fixed, 'value': parguess[2], 'parname': 'gamma','step': 1e-12,
            'limited': [1,1], 'limits': [-50e-13, 50e-13]},
        {'fixed': fixed, 'value': parguess[3], 'parname': 'delta', 'step': 1,
            'limited': [1,1], 'limits': [0, 2048]},
        {'fixed': 1, 'value': parguess[4], 'parname': 'order', 
            'limited': [0,0], 'limits': [3, 7]},
        {'fixed': 1, 'value': parguess[5], 'parname': 'Y',
            'limited': [0,0], 'limits': [0, 2048]}
    ]

    merit_function = Fit.mpfit_residuals(wavelength_model)
    
    lsf = Fit.mpfit_do(merit_function, xs[ok], lines[ok], 
            parinfo, error=slambda[ok])

    delt = np.abs(wavelength_model(lsf.params, xs[ok]) - lines[ok])*1e4

    xsOK = xs[ok]
    linesOK = lines[ok]


    return [ delt, lsf.params, lsf.perror]

def guesslims(spec):
    ''' Guess the spectral limits'''
    f = 1.1

    s = spec.copy()
    s.sort()
    return [-500, s[-10]*f]

class InteractiveSolution:

    header = None
    data = None
    bs = None
    parguess = None
    linelist0 = None
    foundlines = None
    options = None
    slitno = None
    spec = None
    good_solution = False

    ll = None
    pix = None
    xlim = None
    ylim = None
    MAD = None
    STD = None

    first_time = True


    def __init__(self, fig, mfits, linelist, options, slitno, solutions=None):
        self.header = mfits[0]
        self.data = mfits[1]
        self.bs = mfits[2]
        self.options = options
        self.linelist0 = linelist
        self.slitno = slitno
        self.fig = fig

        self.pix = np.arange(2048)

        if solutions is None:
            self.solutions = range(46)
        else:
            print solutions
            self.solutions = solutions
        self.cid = self.fig.canvas.mpl_connect('key_press_event', self)

        self.setup()

    
    def setup(self):
        self.parguess = guess_wavelength_solution(self.slitno, self.header,
                self.bs)
        self.linelist = self.linelist0
        y0 = self.parguess[5]

        S = self.solutions[self.slitno-1]
        if type(S) is not int: # previously setup
            self.MAD = S["MAD"]
            self.STD = S["STD"]
            self.parguess = S["params"]
            self.linelist = S["linelist"]
            self.foundlines = S["foundlines"]
            self.foundlinesig = S["foundlinesig"]
            self.y0 = self.parguess[5]
        else:
            self.MAD = self.STD = self.foundlines = self.linesig = None

        self.spec = np.median(self.data[y0-1:y0+1, :], axis=0) 
            # axis = 0 is spatial direction
        self.ll = wavelength_model(self.parguess, self.pix)
        self.xrng = [min(self.ll) * .99, max(self.ll)*1/.99]
        self.xlim = self.xrng
        self.redraw()

    def draw_found_lines(self):
        pl.subplot(2,1,1)

        pl.grid(True)
        xmin, xmax, ymin, ymax = pl.axis()
        if self.foundlines is not None:
            foundlams = wavelength_model(self.parguess, self.foundlines)
            ok = np.isfinite(self.foundlinesig) 

            for i in xrange(len(self.linelist)):
                if not ok[i]: continue
                D = (foundlams[i] - self.linelist[i])*1e4
                pl.axvline(foundlams[i], color='yellow', ymax=.75)
                pl.text(foundlams[i], 1500, "%1.2f" % D, rotation='vertical',
                        size=12)

            pl.subplot(2,1,2)
            pl.grid(True)

            pl.plot(self.linelist[ok], (foundlams[ok] - self.linelist[ok])*1e4, 
                    'o')


    def draw_residuals(self):
        pl.subplot(2,1,2)

        pl.xlabel(u"$\\lambda$ [$\\mu$ m]")
        pl.ylabel(u"$\\Delta$$\\lambda$ [$\\AA$]")
        gamma = self.parguess[2]
        delta = self.parguess[3]

        x = np.arange(2048)
        pl.plot(self.ll, gamma * (x - delta)**3)

        pl.xlim(self.xlim)

    def draw_vertical_line_marks(self):
        pl.subplot(2,1,1)
        xmin, xmax, ymin, ymax = pl.axis()
        i = 0
        for line in self.linelist:
            pl.axvline(line, color='red', linewidth=.5)

            pl.text(line, ymax*.75, "%5.1f" % (line*1e4), 
                    rotation='vertical', color='black')

            i = i+1
            fwl = self.options['fractional-wavelength-search']
            pl.plot([line*fwl,line/fwl], [0,0], linewidth=2)

    def redraw(self):
        pl.clf()

        pl.subplot(2,1,1)
        pl.subplots_adjust(left=.1,right=.95,bottom=.1,top=.90)
        pl.plot(self.ll, self.spec)

        if self.MAD is None:
            pl.title("[%i] Press 'f' to fit" % self.slitno)
        else:
            a,b,g,d = self.parguess[0:4]
            pl.title(u"[%i] Best fit STD: %0.2f $\AA$, MAD: %0.2f $\AA$: $\\alpha$%0.6f $\\beta$%3.1f $\\gamma$%1.4e $\\delta$%4.1f" % (self.slitno, self.STD, self.MAD, a, b, g, d))


        self.draw_vertical_line_marks()
        self.draw_found_lines()
        self.draw_residuals()

        pl.subplot(2,1,1)
        xmin, xmax, ymin, ymax = pl.axis()
        pl.xlim(self.xlim)
        pl.ylim([-1000, ymax*.8])



    def shift(self, x):
        theline = np.argmin(np.abs(x - self.linelist))

        delt = x - self.linelist[theline] 
        self.parguess[1] -= delt*10
        self.ll = wavelength_model(self.parguess, self.pix)
        self.redraw()

    def drop_point(self, x):
        theline = np.argmin(np.abs(x - self.linelist))
        self.linelist = np.delete(self.linelist, theline)
        if self.foundlines is not None:
            self.foundlines = np.delete(self.foundlines, theline)
            self.foundlinesig = np.delete(self.foundlinesig, theline)
        
    def __call__(self, event):
        kp = event.key
        x = event.xdata
        y = event.ydata

        if x is None: return
        if y is None: return

        if kp == 'x':
            self.xlim = [x*.985,x/.985]

            pl.ion()
            pl.subplot(2,1,1) ; pl.xlim(self.xlim)
            pl.subplot(2,1,2) ; pl.xlim(self.xlim)

        if kp == 'X':
            self.xlim = self.xrng
            pl.ion()
            pl.subplot(2,1,1) ; pl.xlim(self.xlim)
            pl.subplot(2,1,2) ; pl.xlim(self.xlim)

        if kp == '>':
            ''' Fast forward to next uncalib obj '''
            for i in xrange(self.slitno+1, len(self.solutions)):
                if type(self.solutions[i]) is int:
                    self.slitno = i
                    self.setup()
                    break


        if kp == 'G':
            self.slitno += 1
            if self.slitno > 46: 
                print "end limit"
                self.slitno = 46

            self.setup()

        if kp == 'B':
            self.slitno -= 1
            if self.slitno < 1: 
                print "first limit"
                self.slitno = 1
            self.setup()


        if kp == 'Q':
            pl.ioff()
            pl.close(self.fig)

        if kp == 't':
            self.parguess[3] *= 1.05
            self.redraw()


        if kp == 'R':
            self.MAD = None
            self.solutions[self.slitno-1] = self.slitno
            self.setup()


        if kp == 'z':
            self.shift(x)

        if kp == 'd':
            self.drop_point(x)
            self.redraw()

        if kp == 'f':
            [xs, sxs, sigmas] = find_known_lines(self.linelist, self.ll,
                    self.spec, self.options)

            [deltas, params, perror] = fit_model_to_lines(xs, sxs,
                    self.linelist, self.parguess, self.options, False)

            self.ll = wavelength_model(self.parguess, self.pix)
            self.foundlines = xs
            self.foundlinesig = sxs

            ok = np.isfinite(deltas)
            self.STD = np.std(deltas[ok])
            self.MAD = np.median(np.abs(deltas[ok]))

            print "STD: %1.2f MAD: %1.2f" % (self.STD, self.MAD)

            self.parguess = params

            self.solutions[self.slitno-1] = {"params": params, "linelist":
                    self.linelist, "MAD": self.MAD, "foundlines":
                    self.foundlines, "foundlinesig": self.foundlinesig,
                    "sol_1d": [deltas, params, perror, sigmas], "STD":
                    self.STD, "slitno": self.slitno}

            print 'Stored: ', self.solutions[self.slitno-1]['slitno']

            self.redraw()





def fit_wavelength_solution(data, parguess, lines, options, 
        slitno, search_num=145, fixed=False):
    '''Tweaks the guessed parameter values and provides 1d lambda solution
    
    '''

    pix = np.arange(2048.)
    MAD = np.inf

    y0 = parguess[5]
    spec = np.median(data[y0-1:y0+1, :], 
        axis=0) # axis = 0 is spatial direction

    d = 0.1
    dsinbetas = np.sort(np.abs(np.linspace(-d/2., d/2., search_num)))

    sinbetadirection = 1.0

    iteration = 0

    DRAW = False
    if DRAW:
        pl.ion()
        pl.figure(2, figsize=(16,5))
        pl.xlim([2.03,2.3])

    #print "iter  dsb      MAD"
    for dsinbeta in dsinbetas:
        dsinbeta *= sinbetadirection
        sinbetadirection *= -1


        pars = parguess
        pars[1] = parguess[1] + dsinbeta
        ll = wavelength_model(parguess, pix)
        [xs, sxs, sigmas] = find_known_lines(lines, ll, spec, options)
        [deltas, params, perror] = fit_model_to_lines(xs, sxs, lines, 
                pars, options, fixed)

        if DRAW:
            pl.figure(2)
            pl.xlim([1.94,2.1])
            ll2 = wavelength_model(params, pix)
            pl.plot(ll2, spec)
            for line in lines:
                pl.axvline(line ,color='red')
            pl.draw()

        MAD = np.median(deltas)
        iteration += 1
        #print "%2.2i] %3.0i %1.4f %1.4f" % (slitno, iteration, dsinbeta, MAD)


        if MAD > MADLIMIT:
            continue
        else: 
            #print "%i] found: %3i %+1.5f %3.6f" % (slitno, iteration, dsinbeta, MAD)
            break


    if MAD <= MADLIMIT:
        #print("%3i: %3.5f %4.3f %3.3e %4.1f %1.4f" % (slitno, params[0], params[1],
            #params[2], params[3], MAD))

        return [deltas, params, perror, sigmas]
    else:
        print("%3i: Could not find parameters" % slitno)
        return [[], parguess, None, []]

#
# Two dimensional wavelength fitting
#
def fit_outwards_xcor(data, sol_1d, lines, options, slitno):
    lags = np.arange(-5,5)
    pix = np.arange(2048.)

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

            [deltas, params, perror, sigmas] = fit_wavelength_solution( data,
                    params, lines, options, slitno, search_num=25, 
                    fixed=False)

            success = True
            if (len(deltas) < 2) or (np.median(deltas) > .4):
                success = False

            
            if True:
                print ("%2i @ %4i: Success: %i - %3.5f %4.3f %3.3e %2.5f = %2.2f" %
                        (slitno, params[5], success, params[0], params[1],
                            params[2], params[3], np.median(np.abs(deltas))))

            ret.append([params, perror, np.median(deltas), success])

        return ret

    pix = np.arange(2048.)

    params = sweep(xrange(0,10,1))
    params_down = sweep(xrange(-1,-10,-1))

    params.extend(params_down)

    return params

def merge_solutions(lamout, slitno, order, bs, sol_2d, options):
    ff = filter_2d_solutions(sol_2d)
    ar = np.array(map(lambda x: x[0], ff))

    if lamout == None:
        pdb.set_trace()

    if len(ar) == 0:
        print("%3i: no slit level solution has been found." % slitno)
        return lamout


    alphas = ar[:,0]
    betas = ar[:,1]
    gammas = ar[:,2]
    deltas = ar[:,3]
    pixels = ar[:,5]

    alphamodel = np.poly1d(np.polyfit(pixels, alphas, 2))
    betamodel  = np.poly1d(np.polyfit(pixels, betas, 2))
    gammamodel = np.poly1d(np.polyfit(pixels, gammas, 2))
    deltamodel = np.poly1d(np.polyfit(pixels, deltas, 2))

    print "Alpha scatter: %3.3e" % np.std(alphas-alphamodel(pixels))
    print " Beta scatter: %3.3e" % np.std(betas-betamodel(pixels))

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
    ''' Returns wavelength [um] as function of pixel (x)
    
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

def plot_mask_fits(maskname, fname, options):

    from matplotlib.backends.backend_pdf import PdfPages
    path = os.path.join(options["outdir"], maskname)
    print path
    if not os.path.exists(path):
        raise Exception("Output directory '%s' does not exist. This " 
            "directory should exist." % path)

    fp = os.path.join(options['indir'], fname)
    mfits = IO.readmosfits(fp)
    header, data, bs = mfits

    fname = fname.rstrip(".fits")
    solname = os.path.join(path, "mask_solution_%s.npy" % fname)
    maskfit = np.load(solname)[0]

    outname = os.path.join(path, "mask_fit_%s.pdf" % fname)
    pp = PdfPages(outname)


    def plotfun(px, ys, lsf, i):
        pl.plot(px, ys, '.')
        params = np.copy(lsf.params[0:5])
        params[4] = lsf.params[4+i]
        my = mask_model(params, [px])
        pl.plot(px, my)

    alpha_lsf = maskfit["alpha_lsf"]
    beta_lsf = maskfit["beta_lsf"]
    gamma_lsf = maskfit["gamma_lsf"]
    delta_lsf = maskfit["delta_lsf"]

    pixel_set = maskfit["pixels"]
    alpha_set = maskfit["alphas"]
    beta_set = maskfit["betas"]
    gamma_set = maskfit["gammas"]
    delta_set = maskfit["deltas"]

    pixels, alphas, betas, gammas, deltas = map(np.concatenate, [pixel_set,
        alpha_set, beta_set, gamma_set, delta_set])

    pl.clf()
    pl.plot(pixels, alphas, '.')
    pl.ylim([0.9955, 0.9980])
    for i in xrange(len(pixel_set)):
        plotfun(pixel_set[i], alpha_set[i], alpha_lsf, i)
    pp.savefig()

    pl.clf()
    for i in xrange(len(pixel_set)):
        plotfun(pixel_set[i], beta_set[i], beta_lsf, i)
    pp.savefig()

    pl.clf()
    pl.plot(pixels, np.concatenate(gamma_set), '.')
    for i in xrange(len(pixel_set)):
        plotfun(pixel_set[i], gamma_set[i], gamma_lsf, i)
    pp.savefig()

    pl.clf()
    pl.plot(pixels, np.concatenate(delta_set), '.')
    for i in xrange(len(pixel_set)):
        plotfun(pixel_set[i], delta_set[i], delta_lsf, i)
    pp.savefig()

    pp.close()


def plot_sky_spectra(maskname, fname, options):

    from matplotlib.backends.backend_pdf import PdfPages
    path = os.path.join(options["outdir"], maskname)
    if not os.path.exists(path):
        raise Exception("Output directory '%s' does not exist. This " 
            "directory should exist." % path)

    fp = os.path.join(options['indir'], fname)
    mfits = IO.readmosfits(fp)
    header, data, bs = mfits

    
    fname = fname.rstrip(".fits")
    solutions = IO.load_lambdadata(fname, maskname, band, options)

    outname = os.path.join(path, "sky_spectra_%s.pdf" % fname)

    pp = PdfPages(outname)
    band = header['filter'].rstrip()

    # determine region to cutoff spectra for xlims
    linelist = pick_linelist(header)
    hpps = Filters.hpp[band]

    # Pick top 95% of flux for ylims
    sdata = np.sort(data, None)
    ymax = sdata[-15000]


    pix = np.arange(2048)
    for solution in solutions:
        slitno = solution["slitno"]
        parameters = solution["center_sol"][0]
        print "Slit: {0}".format(slitno)

        parguess = guess_wavelength_solution(slitno, header, bs)
        y0 = parguess[-2]

        ll = wavelength_model(parameters, pix)
        measured = data[y0, :]

        pl.clf()
        pl.title("Slit {0}".format(solution["slitno"]))
        pl.plot(ll, measured, linewidth=.2)
        pl.xlim(hpps)
        pl.ylim(-30, ymax)

        for line in linelist:
            pl.axvline(line, color='red', linewidth=.1)

        pp.savefig()

    pp.close()



def plot_data_quality(maskname, fname, options):

    from matplotlib.backends.backend_pdf import PdfPages
    path = os.path.join(options["indir"])
    if not os.path.exists(path):
        raise Exception("Output directory '%s' does not exist. This " 
            "directory should exist." % path)

    fp = os.path.join(path, fname)
    mfits = IO.readmosfits(fp)
    header, data, bs = mfits

    
    fname = fname.rstrip(".fits")
    path = os.path.join(options["outdir"],maskname)
    solname = os.path.join(path, "lambda_coeffs_%s.npy" % fname)
    solutions = np.load(solname)
    solname = os.path.join(path, "mask_solution_%s.npy" % fname)
    masksol = np.load(solname)[0]


    outname = os.path.join(path, "wavelength_fits_%s.pdf" % fname)

    pp = PdfPages(outname)

    filter_fun = (lambda x:
            (x[1] is not None) and
            (x[1][0] < 1e-5) and
            (x[2] < .2) and
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

        if len(ar) == 0: continue

        pixels = ar[:,5]

        alphas = ar[:,0]
        betas = ar[:,1]
        gammas = ar[:,2]
        deltas = ar[:,3]
        sds = ar[:,4]


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
        
        pl.ylim([.993,1/.993])
        pl.xticks(rotation=90)
        pl.ylabel(r'$\alpha$')

        pl.subplot(2,2,2)
        pl.scatter(pixels, betas)
        pl.plot(pixels, betamodel(pixels))
        pl.xticks(rotation=90)
        pl.ylabel(r'$\beta$')

        pl.subplot(2,2,3)
        pl.scatter(pixels, gammas)
        pl.plot(pixels, gammamodel(pixels))
        pl.ylim([0,1e-12])
        pl.xticks(rotation=90)
        pl.ylabel(r'$\gamma$')

        pl.subplot(2,2,4)
        pl.scatter(pixels, deltas)
        pl.plot(pixels, deltamodel(pixels))
        pl.xticks(rotation=90)
        pl.ylabel(r'$\delta$')


        pp.savefig()

    band = header['filter'].rstrip()
    [alpha_pixel, sinbeta_position, sinbeta_pixel, gamma_pixel, 
            delta_pixel] = param_guess_functions(band)
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

    tester()

    if False:
        tester()

        ff = filter(lambda x: (x[1] is not None) and (x[1][0] < 3e-5) and 
                (x[2] < .1) and (x[3] == True), sol_2d)

        ar = np.array(map(lambda x: x[0], ff))
    elif False:
        plot_data_quality("/Users/npk/desktop/c9_reduce/"
        "npk_calib3_q1700_pa_0/lambda_coeffs_m110323_2718.npy")



