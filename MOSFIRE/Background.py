
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
from MOSFIRE import CSU, Fit, IO, Options, Filters, Detector, Wavelength

#from IPython.Shell import IPShellEmbed
#ipshell = IPShellEmbed()

def rem_header_key(header, key):

    try:
        del header[key]
    except:
        return False

    return True


def imcombine(files, maskname, options, flat, outname=None):
    '''
    From a list of files it imcombine returns the imcombine of several values.
    The imcombine code also estimates the readnoise ad RN/sqrt(numreads) so
    that the variance per frame is equal to (ADU + RN^2) where RN is computed
    in ADUs.

    header -- fits header
    ADUs -- The mean # of ADUs per frame
    var -- the Variance [in adu] per frame. 
    bs -- Barset
    itimes -- The _total_ integration time in second
    Nframe -- The number of frames in a stack.

    
    Thus the number of electron per second is derived as: 
        e-/sec = (ADUs * Gain) * (Nframe/itimes)

    The total number of electrons is:
        el = ADUs * Gain * Nframe

    '''

    ADUs = np.zeros((len(files), 2048, 2048))
    itimes = np.ones((len(files), 2048, 2048))
    prevssl = None
    prevmn = None
    patternid = None
    maskname = None

    header = None

    for i in xrange(len(files)):
        fname = files[i]
        thishdr, data, bs = IO.readmosfits(fname, options)
        itimes[i,:,:] *= thishdr["truitime"]
        ADUs[i,:,:] = data.filled(0.0) 

        ''' Construct Header'''
        if header is None:
            header = thishdr
        header.update("imfno%2.2i" % (i), fname, "-------------------")

        map(lambda x: rem_header_key(header, x), ["CTYPE1", "CTYPE2", "WCSDIM",
            "CD1_1", "CD1_2", "CD2_1", "CD2_2", "LTM1_1", "LTM2_2", "WAT0_001",
            "WAT1_001", "WAT2_001", "CRVAL1", "CRVAL2", "CRPIX1", "CRPIX2",
            "RADECSYS"])

        for key in header.keys():
            if key == '': continue
            val = header[key]

            if thishdr.has_key(key):
                if val != thishdr[key]:
                    newkey = "hierarch " + key + ("_img%2.2i" % i)
                    try: header.update(newkey.rstrip(), thishdr[key])
                    except: pass

        ''' Now handle error checking'''

        if maskname is not None:
            if thishdr["maskname"] != maskname:
                raise Exception("File %s uses mask '%s' but the stack is of '%s'" %
                    (fname, thishdr["maskname"], maskname))

        maskname = thishdr["maskname"]
            
        if thishdr["aborted"]:
            raise Exception("Img '%s' was aborted and should not be used" %
                    fname)

        if prevssl is not None:
            if len(prevssl) != len(bs.ssl):
                # todo Improve these checks
                raise Exception("The stack of input files seems to be of "
                        "different masks")
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
                    "this violates an assumption of the DRP. Some new code " 
                    "is needed in the DRP to handle the new units of "
                    "'%s'." % (fname, thishdr["BUNIT"]))

        ''' Error checking is complete'''
        print "%s %s[%s]/%s %5.1f" % (fname, maskname, patternid,
            header['filter'], np.mean(itimes[i]))

    ADUs = np.array(ADUs) / flat
    ADUs_per_sec = ADUs / itimes

    output = np.zeros((2048, 2048))
    exptime = np.zeros((2048, 2048))


    if len(files) > 5:
        print "Drop min/max CRR"
        srt = np.argsort(ADUs_per_sec,axis=0)
        shp = ADUs_per_sec.shape
        sti = np.ogrid[0:shp[0], 0:shp[1], 0:shp[2]]

        ADUs_per_sec = ADUs_per_sec[srt, sti[1], sti[2]]
        ADUs = ADUs[srt, sti[1], sti[2]]
        itimes = itimes[srt, sti[1], sti[2]]

        ADUs = np.mean(ADUs[1:-1,:,:], axis=0)
        ADUs_per_sec = np.mean(ADUs_per_sec[1:-1,:,:], axis=0)
        itimes = np.sum(itimes[1:-1,:,:], axis=0)
        Nframe = len(files) - 2

    else:
        ADUs = np.mean(ADUs, axis=0)
        itimes = np.sum(itimes, axis=0)
        Nframe = len(files) 


    ''' Now handle variance '''
    numreads = header["READS0"]
    RN_adu = Detector.RN / np.sqrt(numreads) / Detector.gain

    var = (ADUs + RN_adu**2) 

    if header.has_key("RN"): raise Exception("RN already populated in header")
    header.update("RN", "%1.3f ADU" % RN_adu)
    header.update("NUMFRM", Nframe)


    if outname is not None:
        header.update("object", "{0}: ADU per frame".format(maskname))
        header.update("bunit", "ADU")

        IO.writefits(ADUs, maskname, "adu_%s" % (outname),
                options, header=header, overwrite=True)

        header.update("object", "{0}: adu^2 var per frame".format(maskname))
        header.update("bunit", "ADU^2")

        IO.writefits(var, maskname, "var_%s" % (outname),
                options, header=header, overwrite=True, lossy_compress=True)

        header.update("object", "{0}: itime s total".format(maskname))
        header.update("bunit", "SECOND")

        IO.writefits(np.float32(itimes), maskname, "itimes_%s" % (outname),
                options, header=header, overwrite=True, lossy_compress=True)

    return header, ADUs, var, bs, itimes, Nframe

def merge_headers(h1, h2):
    """Merge headers h1 and h2 such that h2 has the nod position name
    appended"""
    
    h = h1.copy()

    patternid = h2["frameid"]

    for key in h2.keys():
        val = h2[key]

        if h.has_key(key):
            if val != h[key]:
                newkey = "hierarch " + key + ("_pos_%s" % patternid)
                try: h.update(newkey.rstrip(), val)
                except: pass

    return h


def handle_background(As, Bs, wavenames, maskname, band_name, options): 
    
    global header, bs, edges, data, ivar, lam, sky_sub_out, sky_model_out, band

    band = band_name

    lamname = Wavelength.filelist_to_wavename(wavenames, band_name, maskname,
            options).rstrip(".fits")
    
    suffix = lamname.lstrip("wave_stack_%s_" % band_name)
    
    flatname = os.path.join(maskname, "pixelflat_2d_%s.fits" % band_name)
    hdr, flat = IO.read_drpfits(maskname, 
            "pixelflat_2d_%s.fits" % (band_name), options)

    if np.abs(np.median(flat) - 1) > 0.1:
        raise Exception("Flat seems poorly behaved.")

    hdrA, avgaduA, varaduA, bsA, timeA, NframeA = imcombine(As, maskname,
            options, flat, outname="%s_%s_A.fits" % (band_name, suffix))

    hdrB, avgaduB, varaduB, bsB, timeB, NframeB = imcombine(Bs, maskname,
            options, flat, outname="%s_%s_B.fits" % (band_name, suffix))

    header = merge_headers(hdrA, hdrB)

    elpersecA = avgaduA * Detector.gain * (NframeA/timeA)
    elpersecB = avgaduB * Detector.gain * (NframeB/timeB)
    varA = varaduA * Detector.gain * (NframeA/timeA)**2
    varB = varaduB * Detector.gain * (NframeB/timeA)**2

    AmB = elpersecA - elpersecB
    vAmB = varA + varB
    ivar = 1/vAmB

    sky_sub_out   = np.zeros((2048, 2048), dtype=np.float)
    sky_model_out = np.zeros((2048, 2048), dtype=np.float)

    edges, meta = IO.load_edges(maskname, band, options)
    lam = IO.load_lambdaslit(lamname, maskname, band, options)

    bs = bsA
    data = AmB

    p = Pool()
    solutions = p.map(background_subtract_helper, xrange(len(bs.ssl)))
    p.close()

    xroi = slice(0,2048)
    for sol in solutions:
        if not sol["ok"]: 
            continue

        yroi = slice(sol["bottom"], sol["top"])
        sky_sub_out[yroi, xroi] = sol["output"]
        sky_model_out[yroi, xroi] = sol["model"]
    
    header.update("object", "{0}: electron".format(maskname))
    header.update("bunit", "ELECTRONS")

    IO.writefits(data, maskname, "sub_%s_%s_%s.fits" % (maskname, band,
        suffix), options, header=header, overwrite=True, lossy_compress=True)

    header.update("object", "{0}: electron/frame".format(maskname))
    header.update("bunit", "ELECTRONS")
    IO.writefits(sky_sub_out, maskname, "bsub_%s_%s_%s.fits" % (maskname, band,
        suffix), options, header=header, overwrite=True)

    header.update("object", "{0}: (1/electron)^2".format(maskname))
    header.update("bunit", "(1/ELECTRONS)^2")
    IO.writefits(ivar, maskname, "bsub_ivar_%s_%s_%s.fits" % (maskname, band,
        suffix), options, header=header, overwrite=True, lossy_compress=True)

    header.update("object", "{0}: electron".format(maskname))
    header.update("bunit", "ELECTRONS")

    IO.writefits(sky_model_out, maskname, "bmod_%s_%s_%s.fits" % (maskname,
        band, suffix), options, header=header, overwrite=True,
        lossy_compress=True)

    '''Now create rectified solutions'''
    dlam = Wavelength.grating_results(band)
    hpp = np.array(Filters.hpp[band]) 
    ll_fid = np.arange(hpp[0], hpp[1], dlam)
    nspec = len(ll_fid)


    rectified = np.zeros((2048, nspec), dtype=np.float32)
    rectified_ivar = np.zeros((2048, nspec), dtype=np.float32)

    from scipy.interpolate import interp1d
    for i in xrange(2048):
        ll = lam[1][i,:]
        ss = sky_sub_out[i,:]

        ok = np.isfinite(ll) & np.isfinite(ss) & (ll < hpp[1]) & (ll >
                hpp[0])

        if len(np.where(ok)[0]) < 100:
            continue
        f = interp1d(ll[ok], ss[ok], bounds_error=False)
        rectified[i,:] = f(ll_fid)

        f = interp1d(ll, ivar[i,:], bounds_error=False)
        rectified_ivar[i,:] = f(ll_fid)

    header.update("wat0_001", "system=world")
    header.update("wat1_001", "wtype=linear")
    header.update("wat2_001", "wtype=linear")
    header.update("dispaxis", 1)
    header.update("dclog1", "Transform")
    header.update("dc-flag", 0)
    header.update("ctype1", "AWAV")
    header.update("cunit1", "Angstrom")
    header.update("crval1", ll_fid[0])
    header.update("crval2", 0)
    header.update("crpix1", 1)
    header.update("crpix2", 1)
    header.update("cdelt1", 1)
    header.update("cdelt2", 1)
    header.update("cname1", "angstrom")
    header.update("cname2", "pixel")
    header.update("cd1_1", dlam)
    header.update("cd1_2", 0)
    header.update("cd2_1", 0)
    header.update("cd2_2", 1)


    header.update("object", "rectified [eps]")
    IO.writefits(rectified, maskname, "%s_pair_%s_%s.fits" % (maskname, band_name,
        suffix), options, header=header, overwrite=True, lossy_compress=True)

    header.update("object", "rectified ivar [1/eps^2]")
    IO.writefits(rectified_ivar, maskname, "%s_ivar_pair_%s_%s.fits" %
            (maskname, band_name, suffix), options, header=header, overwrite=True,
            lossy_compress=True)

    header.update("object", "rectified snr")

    IO.writefits(rectified*np.sqrt(rectified_ivar), maskname,
            "%s_sn_pair_%s_%s.fits" % (maskname, band_name, suffix), options,
            header=header, overwrite=True, lossy_compress=True)


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

    global header, bs, edges, data, ivar, lam, sky_sub_out, sky_model_out, band
    tick = time.time()

    # 1
    top = np.int(edges[slitno]["top"](1024))  
    bottom = np.int(edges[slitno]["bottom"](1024)) 
    print "Background subtracting slit %i [%i,%i]" % (slitno, top, bottom)

    pix = np.arange(2048)
    xroi = slice(0,2048)
    yroi = slice(bottom, top)

    slit = data[yroi, xroi]
    ivar[np.logical_not(np.isfinite(ivar))] = 0

    lslit = lam[1][yroi,xroi]

    # 2
    xx = np.arange(slit.shape[1])
    yy = np.arange(slit.shape[0])

    X,Y = np.meshgrid(xx,yy)

    train_roi = slice(5,-5)
    ls = lslit[train_roi].flatten().filled(0)
    ss = slit[train_roi].flatten()
    ys = Y[train_roi].flatten()

    dl = np.median(np.diff(lslit[lslit.shape[0]/2,:]))
    if dl == 0:
        return {"ok": False}

    sort = np.argsort(ls)
    ls = ls[sort]
    ys = ys[sort]

    hpps = np.array(Filters.hpp[band] ) 

    diff = np.append(np.diff(ls), False)

    OK = (diff > 0.001) & (ls > hpps[0]) & (ls < hpps[1]) & (np.isfinite(ls)) \
            & (np.isfinite(ss[sort]))

    if len(np.where(OK)[0]) < 1000:
        print "Failed on slit ", slitno
        return {"ok": False}

    # 3
    pp = np.poly1d([1.0])
    ss = (slit[train_roi] / pp(Y[train_roi])).flatten()
    ss = ss[sort]

    knotstart = max(hpps[0], min(ls[OK])) + 5
    knotend = min(hpps[1], max(ls[OK])) - 5


    for i in range(3):
        try:
            delta = dl*0.9
            knots = np.arange(knotstart, knotend, delta)
            bspline = II.splrep(ls[OK], ss[OK], k=5, task=-1, t=knots)
        except:
            delta = dl*1.4
            knots = np.arange(knotstart, knotend, delta)
            try:
                bspline = II.splrep(ls[OK], ss[OK], k=5, task=-1, t=knots)
            except:
                print "Could not construct spline on slit ", slitno
                return {"ok": False}

        ll = lslit.flatten()
        model = II.splev(ll, bspline)

        oob = np.where((ll < knotstart) | (ll > knotend))
        model[oob] = np.median(ss)
        model = model.reshape(slit.shape)

        output = slit - model

        std = np.abs(output)/(np.sqrt(np.abs(model)+1))

        tOK = (std[train_roi] < 10).flatten() & \
                np.isfinite(std[train_roi]).flatten()  
        OK = OK & tOK[sort]

    return {"ok": True, "slitno": slitno, "bottom": bottom, "top": top,
            "output": output, "model": model, "bspline": bspline}


if __name__ == "__main__":
    background_subtract()


