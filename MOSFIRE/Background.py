
import os
import sys
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

def imcombine(files, maskname, options, flat, outname=None):

    datas = np.zeros((len(files), 2048, 2048))
    itimes = []
    prevssl = None
    prevmn = None
    patternid = None
    maskname = None

    header = None

    for i in xrange(len(files)):
        fname = files[i]
        thishdr, data, bs = IO.readmosfits(fname, options)
        itimes.append(thishdr["truitime"])
        datas[i,:,:] = data.filled(0)

        ''' Construct Header'''
        if header is None:
            header = thishdr
        header.update("imfno%2.2i" % (i), fname, "-------------------")
        for key in header.keys():
            val = header[key]

            if thishdr.has_key(key):
                if val != thishdr[key]:
                    newkey = "hierarch " + key + ("_img%2.2i" % i)
                    header.update(newkey.rstrip(), thishdr[key])

        ''' Now handle error checking'''

        if thishdr["aborted"]:
            raise Exception("Img '%s' was aborted and should not be used" %
                    fname)

        if prevssl is not None:
            if len(prevssl) != len(bs.ssl):
                # todo Improve these checks
                raise Exception("The stack of input files seems to be of "
                        " different masks")
        prevssl = bs.ssl

        if patternid is not None:
            if patternid != thishdr["frameid"]:
                raise Exception("The stack should be of '%s' frames only, but "
                        "the current image is a '%s' frame." % (patternid, 
                            thishdr["frameid"]))

        patternid = thishdr["frameid"]

        if maskname is not None:
            if maskname != thishdr["maskname"]:
                raise Exception("The stack should be of CSU mask '%s' frames "
                        "only but contains a frame of '%s'." % (maskname,
                        thishdr["maskname"]))

        maskname = thishdr["maskname"]

        if thishdr["BUNIT"] != "ADU per coadd":
            raise Exception("The units of '%s' are not in ADU per coadd and "
                    "this violates an assumption of the DRP. Some code " 
                    "is needed in the DRP to handle this case." % fname)

        ''' Error checking is complete'''
        print "%s %s[%s] %5.1f" % (fname, maskname, patternid, itimes[-1])

    
    itimes = np.array(itimes)
    if np.any(itimes[0] != itimes):
        raise Exception("Nod position integration times are not all the same "
                "current DRP cannot handle this situation properly")
    
    datas = np.array(datas) * Detector.gain
    dd = np.sort(datas,axis=0)
    output = np.zeros((2048, 2048))
    cnts = np.zeros((2048,2048),np.int16) + len(datas)
    mask = np.zeros((2048,2048),np.int16) + data.mask

    if len(files) > 2:
        mn = np.mean(datas,axis=0)
        tots = np.sum(datas,axis=0)
        md = np.median(datas,axis=0)
        sd = np.std(datas,axis=0)
    
        issx, issy = np.where((mn-md)/md > 1)
        for i in xrange(len(issx)):
            x = issx[i] ; y = issy[i]
            tots[x,y] = np.sum(dd[1:-1,x,y])
            cnts[x,y] -= 2
            mask[x,y] += 2

    itime = itimes[0]
    rates = tots / (cnts*itime*flat)
    var = itime*cnts*flat*(1/tots + 1/Detector.RN**2)

    path = os.path.join(options["outdir"], maskname)
    if not os.path.exists(path):
        raise Exception("Output directory '%s' does not exist. This " 
                "directory should exist." % path)

    if outname is not None:
        IO.writefits(rates, maskname, outname, options, header=header,
                overwrite=True)
        IO.writefits(mask, maskname, "mask_" + outname, options, header=header,
                overwrite=True)

    return header, rates, var, mask, bs
    

def handle_background(As, Bs, lamname, maskname, band_name, options): 
    
    global header, bs, edges, data, lam, sky_sub_out, sky_model_out, band

    band = band_name

    
    flatname = os.path.join(maskname, "pixelflat_2d_%s.fits" % band_name)
    hdr, flat = IO.read_drpfits(maskname, 
            "pixelflat_2d_%s.fits" % (band_name), options)

    hdrA, ratesA, varA, maskA, bsA = imcombine(As, maskname, options, 
            flat, outname="A.fits")
    hdrB, ratesB, varB, maskB, bsB = imcombine(Bs, maskname, options, 
            flat, outname="B.fits")

    AmB = ratesA - ratesB
    vAmB = varA + varB

    sky_sub_out   = np.zeros((2048, 2048), dtype=np.float)
    sky_model_out = np.zeros((2048, 2048), dtype=np.float)

    edges = IO.load_edges(maskname, band, options)
    lam = IO.load_lambdaslit(lamname, maskname, band, options)

    header = hdrA
    bs = bsA
    data = AmB

    p = Pool()
    solutions = p.map(background_subtract_helper, xrange(len(bs.ssl)))
    p.close()

    xroi = slice(0,2048)
    for sol in solutions:
        if not sol["ok"]: continue

        yroi = slice(sol["bottom"], sol["top"])
        sky_sub_out[yroi, xroi] = sol["output"]
    
    IO.writefits(sky_sub_out, maskname, "bsub_%s_%s.fits" % (maskname, band),
            options, header=header, overwrite=True)


    if True:
        dlam = np.median(np.diff(lam[1][1024,:]))
        print dlam
        hpp = np.array(Filters.hpp[band]) 
        print hpp
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
            ll = lam[1][i,:]*1e4
            ss = sky_sub_out[i,:]
            iv = 1/vAmB[i,:]

            ok = np.isfinite(ll) & np.isfinite(ss) & (ll < hpp[1]) & (ll >
                    hpp[0])

            if len(np.where(ok)[0]) < 100:
                continue
            f = interp1d(ll[ok], ss[ok], bounds_error=False)
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

        hdu = pf.PrimaryHDU(1/vAmB)
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
    2. Convert the 2d spectrum into a 1d spectrum
    3. Estimate transmission function

    '''

    global header, bs, edges, data, lam, sky_sub_out, sky_model_out, band
    tick = time.time()

    # 1
    top = np.int(edges[slitno]["top"](1024)) - 7
    bottom = np.int(edges[slitno]["bottom"](1024)) + 7
    print "Background subtracting slit %i [%i,%i]" % (slitno, top, bottom)

    pix = np.arange(2048)
    xroi = slice(0,2048)
    yroi = slice(bottom, top)

    slit = data[yroi, xroi]
    ivar = 1/(slit)
    ivar[np.logical_not(np.isfinite(ivar))] = 0

    lslit = lam[1][yroi,xroi]

    # 2
    xx = np.arange(slit.shape[1])
    yy = np.arange(slit.shape[0])

    X,Y = np.meshgrid(xx,yy)

    ls = lslit.flatten().filled(0)
    ss = slit.flatten()
    ys = Y.flatten()

    sort = np.argsort(ls)
    ls = ls[sort]
    ys = ys[sort]

    hpps = np.array(Filters.hpp[band] ) * 1e4

    diff = np.append(np.diff(ls), False)

    OK = (diff > 0.01) & (ls > hpps[0]) & (ls < hpps[1]) & (np.isfinite(ls)) \
            & (np.isfinite(ss[sort]))

    if len(np.where(OK)[0]) < 1000:
        print "Failed on slit ", slitno
        return {"ok": False}

    # 3
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
            delta = 5.0
            knots = np.arange(knotstart, knotend, delta)
            try:
                bspline = II.splrep(ls[OK], ss[OK], k=3, task=-1, t=knots)
            except:
                "Could not construct spline on slit ", slitno
                return {"ok": False}

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


