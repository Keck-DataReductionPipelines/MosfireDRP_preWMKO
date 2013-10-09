
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



def handle_rectification(maskname, nod_posns, wavenames, band_pass, options,
        commissioning_shift=3.0):
    global edges, dats, ivars, shifts, lambdas, band, fidl


    '''
    '''

    band = band_pass
    lname = Wavelength.filelist_to_wavename(wavenames, band_pass, maskname,
            options).rstrip(".fits")

    suffix = lname.lstrip("wave_stack_%s_" % band_pass)

    dlambda = Wavelength.grating_results(band)

    hpp = Filters.hpp[band]
    fidl = np.arange(hpp[0], hpp[1], dlambda)

    lambdas = IO.load_lambdaslit(lname, maskname, band, options)
    edges, meta = IO.load_edges(maskname, band, options)

    shifts = []
    for pos in nod_posns:

        II = IO.read_drpfits(maskname, "adu_{0}_{1}_{2}.fits".format(band,
            suffix, pos), options)

        off = np.array((II[0]["decoff"], II[0]["raoff"]),dtype=np.float64)
        if II[0].has_key("yoffset"):
            off = -II[0]["yoffset"]
        else:
            # Deal with data taken during commissioning
            if II[0]["frameid"] == 'A': off = 0.0
            else: off = commissioning_shift

        try: off0
        except: off0 = off

        shift = off - off0

        shifts.append(shift)
        print "Position {0} shift: {1:2.2f} as".format(pos, shift)
    

    theBPM = IO.badpixelmask()

    fname = "bsub_{0}_{1}_{2}.fits".format(maskname, band, suffix)
    EPS = IO.read_drpfits(maskname, fname, options)
    EPS[1] = np.ma.masked_array(EPS[1], theBPM, fill_value=0)

    fname = "bsub_ivar_{0}_{1}_{2}.fits".format(maskname, band, suffix)
    IVAR = IO.read_drpfits(maskname, fname, options)
    IVAR[1] = np.ma.masked_array(IVAR[1], theBPM, fill_value=0)

    dats = EPS
    ivars = IVAR
    EPS[0].update("ORIGFILE", fname)

    tock = time.time()
    sols = range(len(edges)-1,-1,-1)

    p = Pool()
    solutions = p.map(handle_rectification_helper, sols)
    p.close()
    tick = time.time()
    print "-----> Mask took %i. Writing to disk." % (tick-tock)

    output = np.zeros((1, len(fidl)))
    snrs = np.zeros((1, len(fidl)))
    ivarout= np.zeros((1, len(fidl)))


    # the barset [bs] is used for determining object position
    drop, drop, bs = IO.readmosfits(wavenames[0], options)


    for i in xrange(len(solutions)):
        solution = solutions[i]
        header = EPS[0].copy()

        pixel_dist = np.float(bs.ssl[-(i+1)]['Target_to_center_of_slit_distance'])/0.18

        header.update("OBJECT", "{0}".format(solution["Target_Name"]))

        ll = solution["lambda"]

        header.update("wat0_001", "system=world")
        header.update("wat1_001", "wtype=linear")
        header.update("wat2_001", "wtype=linear")
        header.update("dispaxis", 1)
        header.update("dclog1", "Transform")
        header.update("dc-flag", 0)
        header.update("ctype1", "AWAV")
        header.update("cunit1", "Angstrom")
        header.update("crval1", ll[0])
        header.update("crval2", -solution["eps_img"].shape[0]/2 - pixel_dist)
        header.update("crpix1", 1)
        header.update("crpix2", 1)
        header.update("cdelt1", 1)
        header.update("cdelt2", 1)
        header.update("cname1", "angstrom")
        header.update("cname2", "pixel")
        header.update("cd1_1", ll[1]-ll[0])
        header.update("cd1_2", 0)
        header.update("cd2_1", 0)
        header.update("cd2_2", 1)

        img = solution["eps_img"]
        ivar = solution["iv_img"]

        S = output.shape
        output  = np.append(output, img, 0)
        output = np.append(output, np.nan*np.zeros((3,S[1])), 0)
        snrs = np.append(snrs, img*np.sqrt(ivar), 0)
        snrs = np.append(snrs, np.nan*np.zeros((3,S[1])), 0)
        ivarout = np.append(ivarout, ivar, 0)
        ivarout = np.append(ivarout, np.nan*np.zeros((3,S[1])), 0)
        

        if True:
            IO.writefits(solution["eps_img"], maskname,
                "eps_{0}_{1}_S{2:02g}.fits".format(band, suffix, i+1), options,
                overwrite=True, header=header, lossy_compress=False)

            IO.writefits(solution["iv_img"], maskname,
                "ivar_{0}_{1}_S{2:02g}.fits".format(band, suffix, i+1), options,
                overwrite=True, header=header, lossy_compress=False)

    header = EPS[0].copy()
    header.update("OBJECT", "{0}/{1}".format(maskname, band))
    header.update("wat0_001", "system=world")
    header.update("wat1_001", "wtype=linear")
    header.update("wat2_001", "wtype=linear")
    header.update("dispaxis", 1)
    header.update("dclog1", "Transform")
    header.update("dc-flag", 0)
    header.update("ctype1", "AWAV")
    header.update("cunit1", "Angstrom")
    header.update("crval1", ll[0])
    header.update("crval2", 1)
    header.update("crpix1", 1)
    header.update("crpix2", 1)
    header.update("cdelt1", 1)
    header.update("cdelt2", 1)
    header.update("cname1", "angstrom")
    header.update("cname2", "pixel")
    header.update("cd1_1", ll[1]-ll[0])
    header.update("cd1_2", 0)
    header.update("cd2_1", 0)
    header.update("cd2_2", 1)

    IO.writefits(output, maskname, "eps_{0}_{1}_{2}.fits".format(maskname,
        suffix, band), options, overwrite=True, header=header,
        lossy_compress=False)

    IO.writefits(snrs, maskname, "snrs_{0}_{1}_{2}.fits".format(maskname,
        suffix, band), options, overwrite=True, header=header,
        lossy_compress=False)

    IO.writefits(ivarout, maskname, "ivars_{0}_{1}_{2}.fits".format(maskname,
        suffix, band), options, overwrite=True, header=header,
        lossy_compress=False)


def r_interpol(ls, ss, lfid, tops, top, shift_pix=0, pad=[0,0], fill_value=0.0):
    '''
    Interpolate the data ss(ls, fs) onto a fiducial wavelength vector.
    ls[n_spatial, n_lam] - wavelength array
    ss[n_spatial, n_lam] - corresponding data array
    lfid[n_lam] - wavelength fiducial to interpolate onto
    shift_pix - # of pixels to shift in spatial direction
    pad - # of pixels to pad in spatial direction
    fill_value - passed through to interp1d
    '''

    S = ss.shape

    output = np.zeros((np.int(S[0]+pad[0]+pad[1]), len(lfid)))

    L = np.double(len(lfid))
    
    # First interpolate onto a common wavelength grid
    for i in xrange(S[0]):

        ll = ls[i,:] ; sp = ss[i,:]
        ok = np.where(ll>1000)[0]

        if len(ok) < 100: continue

        f = II.interp1d(ll[ok], sp[ok], bounds_error=False, fill_value = fill_value)
            

        output[i+pad[0],:] = f(lfid)

    # Now rectify in spatial
    vert_shift = tops-top-shift_pix

    f = II.interp1d(ls[10, :], vert_shift, bounds_error=False, fill_value = fill_value)

    for i in xrange(output.shape[1]):
        to_shift = f(fidl[i])
        x = np.arange(output.shape[0])
        y = II.interp1d(x, output[:, i], bounds_error=False, fill_value=0.0)

        output[:,i] = y(x + to_shift)


            
    return output


def handle_rectification_helper(edgeno):
    global edges, dats, ivars, shifts, lambdas, band, fidl

    pix = np.arange(2048)
    
    edge = edges[edgeno]

    print "Handling edge: ", edge["Target_Name"]

    tops = edge["top"](pix)
    bots = edge["bottom"](pix)

    # Length of the slit in arcsecond
    lenas = (tops[1024] - bots[1024]) * 0.18
    mxshift = np.abs(np.int(np.ceil(np.max(shifts)/0.18)))
    mnshift = np.abs(np.int(np.floor(np.min(shifts)/0.18)))

    top = min(np.floor(np.min(tops)), 2048)
    bot = max(np.ceil(np.max(bots)), 0)

    ll = lambdas[1].data[bot:top, :]
    eps = dats[1][bot:top, :].filled(0.0)
    ivs = ivars[1][bot:top, :].filled(0.0)

    lmid = ll[ll.shape[0]/2,:]
    hpp = Filters.hpp[band]

    minl = lmid[0] if lmid[0]>hpp[0] else hpp[0]
    maxl = lmid[-1] if lmid[-1]<hpp[1] else hpp[1]


    epss = []
    vss = []
    sign = -1
    for shift in shifts:
        output = r_interpol(ll, eps, fidl, tops, top, shift_pix=shift/0.18,
            pad=[mnshift, mxshift])
        epss.append(sign * output)

        var = 1/ivs
        output = r_interpol(ll, var, fidl, tops, top, shift_pix=shift/0.18,
            pad=[mnshift, mxshift], fill_value=np.inf) 
        vss.append(output)

        sign *= -1

    eps_img = np.sum(epss, axis=0)


    # Remove any NaNs or infs from the variance array
    ivs = []
    for var in vss:
        THIS_IVAR = 1/var
        bad = np.where(np.isfinite(THIS_IVAR) == 0)
        THIS_IVAR[bad] = 0
        ivs.append(THIS_IVAR)

    iv_img = np.sum(np.array(ivs), axis=0)


    return {"eps_img": eps_img, "iv_img": iv_img, "lambda": fidl,
            "Target_Name": edge["Target_Name"], "slitno": edgeno+1}

