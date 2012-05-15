
import os
import time

import numpy as np
import pylab as pl
import pyfits as pf
from multiprocessing import Pool
import scipy as sp
import scipy.ndimage
from scipy import interpolate as II

import pdb

import MOSFIRE
from MOSFIRE import CSU, Fit, IO, Options, Filters, Detector

#from IPython.Shell import IPShellEmbed
#ipshell = IPShellEmbed()


def handle_background(datalist, lamname, maskname, band_name, extension, 
        options):
    global header, bs, edges, data, lam, sky_sub_out, sky_model_out, band

    band = band_name

    path = os.path.join(options["outdir"], maskname)
    if not os.path.exists(path):
        raise Exception("Output directory '%s' does not exist. This " 
                "directory should exist." % path)
    
    sky_sub_out   = np.zeros((2048, 2048), dtype=np.float)
    sky_model_out = np.zeros((2048, 2048), dtype=np.float)

    edges = IO.load_edges(maskname, band, options)
    lam = IO.load_lambdaslit(lamname, maskname, band, options)

    for fname in datalist:
        fp = os.path.join(path, fname)

        header, data, bs = IO.readmosfits(fp,
                extension=os.path.join(path,extension))
        ivar = 1/(data + 25)

        data[np.logical_not(np.isfinite(data))] = 0.
    
        p = Pool()
        solutions = map(background_subtract_helper, xrange(len(bs.ssl)))
        p.close()

        xroi = slice(0,2048)
        for sol in solutions:
            if not sol["ok"]: continue

            yroi = slice(sol["bottom"], sol["top"])
            sky_sub_out[yroi, xroi] = sol["output"]
        
        print path
        output = os.path.join(path, 'b' + fname.lstrip(','))
        hdu = pf.PrimaryHDU(sky_sub_out)
        try: os.remove(output)
        except: pass
        hdu.writeto(output)


    if True:
        dlam = np.median(np.diff(lam[1][1024,:]))
        hpp = np.array(Filters.hpp[band]) * 1e4
        ll_fid = np.arange(hpp[0], hpp[1], dlam)
        nspec = len(ll_fid)

        f = open("ll_fid.txt", "w")
        for i in xrange(nspec):
            f.write("%4.4i %6.6f\n" % (i, ll_fid[i]))

        f.close()


        rectified = np.zeros((2048, nspec), dtype=np.float32)
        rectified_ivar = np.zeros((2048, nspec), dtype=np.float32)

        from scipy.interpolate import interp1d
        for i in xrange(2048):
            ll = lam[1][i,:]
            ss = sky_sub_out[i,:]
            iv = ivar[i,:]

            f = interp1d(ll, ss, bounds_error=False)
            rectified[i,:] = f(ll_fid)

            f = interp1d(ll, iv, bounds_error=False)
            rectified_ivar[i,:] = f(ll_fid)

        hdu = pf.PrimaryHDU(rectified)
        hdu.header.update("wat0_001", "system=world")
        hdu.header.update("wat1_001", "wtype=linear")
        hdu.header.update("wat2_001", "wtype=linear")
        hdu.header.update("dispaxis", 2)
        hdu.header.update("dclog1", "Transform")
        hdu.header.update("dc-flag", 0)
        hdu.header.update("ctype1", "LINEAR")
        hdu.header.update("ctype2", "LINEAR")
        hdu.header.update("crval1", ll_fid[0])
        hdu.header.update("crval2", 1)
        hdu.header.update("cdelt1", dlam)
        hdu.header.update("cdelt2", 1)

        path = os.path.join(options["outdir"], maskname)
        fn = os.path.join(path, "rectified_test.fits")
        try: os.remove(fn)
        except: pass
        hdu.writeto(fn)

        hdu = pf.PrimaryHDU(ivar)
        fn = os.path.join(path, "ivar_test.fits")
        try: os.remove(fn)
        except: pass
        hdu.writeto(fn)

        hdu = pf.PrimaryHDU(rectified*np.sqrt(rectified_ivar))
        fn = os.path.join(path, "rectified_sn_test.fits")
        try: os.remove(fn)
        except: pass
        hdu.writeto(fn)





def background_subtract_helper(slitno):
    '''

    Background subtraction follows the methods outlined by Kelson (2003). Here
    a background is estimated as a function of wavelength using B-splines and
    subtracted off. The assumption is that background is primarily a function
    of wavelength, and thus by sampling the background across the full 2-d
    spectrum the background is sampled at much higher than the native spectral
    resolution of mosfire. 

    Formally, the assumption that background is only a function of wavelength
    is incorrect, and indeed a "transmission function" is estimated from the 2d
    spectrum. This gives an estimate of the throughput of the slit and divided
    out.

    1. Extract the slit from the 2d image.
    2. Find discrepant pixels and remove them from the set of background
    3. Convert the 2d spectrum into a 1d spectrum
    4. Estimate transmission function

    '''

    global header, bs, edges, data, lam, sky_sub_out, sky_model_out, band
    tick = time.time()

    # 1
    top = np.int(edges[slitno]["top"](1024)) - 5
    bottom = np.int(edges[slitno]["bottom"](1024)) + 6
    print "Background subtracting slit %i [%i,%i]" % (slitno, top, bottom)

    pix = np.arange(2048)
    xroi = slice(0,2048)
    yroi = slice(bottom, top)

    slit = data[yroi, xroi]
    ivar = 1/(slit)
    ivar[np.logical_not(np.isfinite(ivar))] = 0

    lslit = lam[1][yroi,xroi]

    # 2
    median_slit = sp.ndimage.median_filter(slit, size=(3,3))

    deviations = np.abs(slit - median_slit)/np.sqrt(np.abs(median_slit))
    deviations[np.logical_not(np.isfinite(deviations))] = 0


    # 3
    xx = np.arange(slit.shape[1])
    yy = np.arange(slit.shape[0])

    X,Y = np.meshgrid(xx,yy)

    ls = lslit.flatten()
    ss = slit.flatten()
    ys = Y.flatten()

    sort = np.argsort(ls)
    ls = ls[sort]
    ys = ys[sort]

    hpps = np.array(Filters.hpp[band] ) * 1e4

    diff = np.append(np.diff(ls), False)
    OK = (diff > 0.001) & (deviations.flatten() < 15) & (ls > hpps[0]) & \
        (ls < hpps[1]) 

    if len(np.where(OK)[0]) < 1000:
        print "Failed on slit ", slitno
        return {"ok": False}

    # 4
    pp = np.poly1d([1.0])

    ss = (slit / pp(Y)).flatten()
    ss = ss[sort]

    knotstart = max(hpps[0], min(ls[OK])) + 5
    knotend = min(hpps[1], max(ls[OK])) - 5

    for i in range(3):
        try:
            delta = 1.5
            knots = np.arange(knotstart, knotend, delta)
            bspline = II.splrep(ls[OK], ss[OK], k=3, task=-1, t=knots)
        except:
            delta = 3.0
            knots = np.arange(knotstart, knotend, delta)
            bspline = II.splrep(ls[OK], ss[OK], k=3, task=-1, t=knots)
            #ipshell()

        ll = lslit.flatten()
        model = II.splev(ll, bspline)

        oob = np.where((ll < knotstart) | (ll > knotend))
        model[oob] = np.median(ss)
        model = model.reshape(slit.shape)
        model *= pp(Y)

        output = slit - model

        std = np.abs(output)/(np.sqrt(np.abs(model)+1))
        tOK = (std < 30).flatten() & np.isfinite(std).flatten()
        OK = OK & tOK[sort]

    return {"ok": True, "slitno": slitno, "bottom": bottom, "top": top,
            "output": output, "model": model, "bspline": bspline}


if __name__ == "__main__":
    background_subtract()


