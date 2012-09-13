
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

        II = IO.read_drpfits(maskname, "eps_{0}_{1}_{2}.fits".format(band,
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
    
    fname = "bsub_{0}_{1}_{2}.fits".format(maskname, band, suffix)
    EPS = IO.read_drpfits(maskname, fname, options)

    fname = "bsub_ivar_{0}_{1}_{2}.fits".format(maskname, band, suffix)
    IVAR = IO.read_drpfits(maskname, fname, options)

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
    for i in xrange(len(solutions)):
        solution = solutions[i]
        header = EPS[0].copy()
        
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
        header.update("crval2", 0)
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

        IO.writefits(solution["eps_img"], maskname,
                "eps_{0}_{1}_S{2:02g}.fits".format(band, suffix, i+1), options,
                overwrite=True, header=header, lossy_compress=True)

        
        IO.writefits(solution["iv_img"], maskname,
                "ivar_{0}_{1}_S{2:02g}.fits".format(band, suffix, i+1), options,
                overwrite=True, header=header, lossy_compress=True)

    header.update("OBJECT", "{0}/{1}".format(maskname, band))

    IO.writefits(output, maskname, "eps_{0}_{1}_{2}.fits".format(maskname,
        suffix, band), options, overwrite=True, header=header,
        lossy_compress=True)

    IO.writefits(snrs, maskname, "snrs_{0}_{1}_{2}.fits".format(maskname,
        suffix, band), options, overwrite=True, header=header,
        lossy_compress=True)

def r_interpol(ls, ss, lfid, shift_pix=0, pad=[0,0]):
    '''
    Interpolate the data ss(ls, fs) onto a fiducial wavelength vector.
    ls[n_spatial, n_lam] - wavelength array
    ss[n_spatial, n_lam] - corresponding data array
    lfid[n_lam] - wavelength fiducial to interpolate onto
    shift_pix - # of pixels to shift in spatial direction
    pad - # of pixels to pad in spatial direction
    '''

    S = ss.shape

    output = np.zeros((np.int(S[0]+pad[0]+pad[1]), len(lfid)))

    L = np.double(len(lfid))
    
    # First interpolate onto a common wavelength grid
    for i in xrange(S[0]):

        ll = ls[i,:] ; sp = ss[i,:]
        ok = np.where(ll>1000)[0]

        if len(ok) < 100: continue

        f = II.interp1d(ll[ok], sp[ok], bounds_error=False, fill_value = 0.0)
        output[i+pad[0],:] = f(lfid)

    # Now shift in the spatial direciton
    if np.abs(shift_pix) > 1e-4:
        if np.abs(shift_pix - np.int(shift_pix)) < 1e-4:
            output = np.roll(output, -shift, axis=0)
        else:
            y = np.arange(output.shape[0])
            for i in xrange(output.shape[1]):

                f = II.interp1d(y, output[:, i], bounds_error=False, fill_value
                        = 0.0)

                output[:,i] = f(y-shift_pix)
            
    return output


def handle_rectification_helper(edgeno):
    global edges, dats, ivars, shifts, lambdas, band, fidl

    pix = np.arange(2048)
    
    edge = edges[edgeno]

    print "Handling edge: ", edge["Target_Name"]

    tops = edge["top"](pix)
    bots = edge["bottom"](pix)

    lenas = (tops[1024] - bots[1024]) * 0.18
    mxshift = np.int(np.ceil(np.max(shifts)/0.18))
    mnshift = np.int(np.floor(np.min(shifts)/0.18))

    top = min(np.floor(np.min(tops)), 2048)
    bot = max(np.ceil(np.max(bots)), 0)

    ll = lambdas[1].data[bot:top, :]
    eps = dats[1][bot:top, :]
    ivs = ivars[1][bot:top, :]

    lmid = ll[ll.shape[0]/2,:]
    hpp = Filters.hpp[band]

    minl = lmid[0] if lmid[0]>hpp[0] else hpp[0]
    maxl = lmid[-1] if lmid[-1]<hpp[1] else hpp[1]


    epss = []
    ivss = []
    sign = -1
    for shift in shifts:

        output = r_interpol(ll, eps, fidl, shift_pix=shift/0.18, pad=[mnshift,
            mxshift])
        epss.append(sign * output)

        output = r_interpol(ll, ivs, 
            fidl, shift_pix=shift/0.18, pad=[mnshift, mxshift])
        ivss.append(output)

        sign *= -1

    eps_img = np.sum(epss, axis=0)
    iv_img = 1/np.sum(1/np.array(ivss), axis=0)

    return {"eps_img": eps_img, "iv_img": iv_img, "lambda": fidl,
            "Target_Name": edge["Target_Name"], "slitno": edgeno+1}

