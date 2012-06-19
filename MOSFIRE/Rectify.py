
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

    orders = {"Y": 6, "J": 5, "H": 4, "K": 3}
    order = orders[band]
    d = 1e3/110.5 # Groove spacing in micron
    pixelsize, focal_length = 18.0, 250e3 # micron
    scale = pixelsize/focal_length
    dlambda = scale * d / order * 10000

    hpp = Filters.hpp[band]
    fidl = np.arange(hpp[0], hpp[1], dlambda)

    lambdas = IO.load_lambdaslit(lname, maskname, band, options)
    edges, meta = IO.load_edges(maskname, band, options)

    shifts = []
    for pos in nod_posns:
        II = IO.read_drpfits(maskname, "eps_{0}_{1}.fits".format(band, pos),
                options)
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
    
    fname = "bsub_{0}_{1}.fits".format(maskname, band)
    EPS = IO.read_drpfits(maskname, fname, options)
    fname = "bsub_ivar_{0}_{1}.fits".format(maskname, band)
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
        ff = solution["frac"]

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


        # TODO: 15:-1 is a hack. there is a deeper problem that needs fixing.
        img = solution["eps_img"][15:-1]
        ivar = solution["iv_img"][15:-1]

        output  = np.append(output, img, 0)
        snrs = np.append(snrs, img*np.sqrt(ivar), 0)

        IO.writefits(solution["eps_img"], maskname,
                "eps_{0}_S{1:02g}.fits".format(band, i+1), options,
                overwrite=True, header=header, lossy_compress=True)

        
        IO.writefits(solution["iv_img"], maskname,
                "ivar_{0}_S{1:02g}.fits".format(band, i+1), options,
                overwrite=True, header=header, lossy_compress=True)

    header.update("OBJECT", "{0}/{1}".format(maskname, band))

    IO.writefits(output, maskname, "eps_{0}_{1}.fits".format(maskname, band),
            options, overwrite=True, header=header, lossy_compress=True)

    IO.writefits(snrs, maskname, "snrs_{0}_{1}.fits".format(maskname, band),
            options, overwrite=True, header=header, lossy_compress=True)

def r_interpol(ls, fs, ss, lfid, ffid):
    '''Interpolate the data ss(ls, fs) onto a grid that looks like 
    
    ^ ffid
    |
    o---> lfid
    '''

    S = ss.shape

    output = np.zeros((len(ffid), len(lfid)))
    
    L = np.double(len(lfid))
    lam_to_pix = II.interp1d(lfid, np.arange(L)/L)

    for i in xrange(len(ffid)):
        ll = ls[i,:] ; sp = ss[i,:]
        ok = np.where(ll>1000)[0]

        if len(ok) < 100: continue

        f = II.interp1d(ll[ok], sp[ok], bounds_error=False)
        output[i,:] = f(lfid)
    
    for i in xrange(len(lfid)):

        pix = lam_to_pix([lfid[i]])[0]
        fr = fs[:,np.round(pix*2048.)] 
        sp = output[:,np.round(pix*len(lfid))]

        f = II.interp1d(fr, sp, bounds_error=False)
        output[:,i] = f(ffid)

    return output


def handle_rectification_helper(edgeno):
    global edges, dats, ivars, shifts, lambdas, band, fidl

    pix = np.arange(2048)
    
    edge = edges[edgeno]

    print "Handling edge: ", edge["Target_Name"]

    tops = edge["top"](pix)
    bots = edge["bottom"](pix)

    lenas = (tops[1024] - bots[1024]) * 0.18
    mxshift = np.ceil(np.max(shifts)/0.18)

    top = min(np.floor(np.min(tops)) + mxshift, 2048)
    bot = max(np.ceil(np.max(bots)) - mxshift, 0)

    slopes = 1.0/(tops-bots)
    X,Y = np.mgrid[bot:top, 0:2048]
    fracs = (X.copy() - bots) * slopes
    ll = lambdas[1].data[bot:top, :]
    eps = dats[1][bot:top, :]
    ivs = ivars[1][bot:top, :]

    fidf = fracs[:,1024]
    lmid = ll[ll.shape[0]/2,:]
    hpp = Filters.hpp[band]

    minl = lmid[0] if lmid[0]>hpp[0] else hpp[0]
    maxl = lmid[-1] if lmid[-1]<hpp[1] else hpp[1]


    epss = []
    ivss = []
    sign = 1
    for shift in shifts:
        epss.append(sign * r_interpol(ll, fracs, eps, 
            fidl, fidf - shift/lenas))

        ivss.append(r_interpol(ll, fracs, ivs, 
            fidl, fidf - shift/lenas))

        sign *= -1

    eps_img = np.sum(epss, axis=0)
    iv_img = 1/np.sum(1/np.array(ivss), axis=0)

    return {"eps_img": eps_img, "iv_img": iv_img, "lambda": fidl,
            "Target_Name": edge["Target_Name"], "slitno": edgeno+1, "frac":
            fidf}

