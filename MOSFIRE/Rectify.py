
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



def handle_rectification(maskname, nod_posns, lname, band_pass, options):
    global edges, dats, ivars, shifts, lambdas, band

    band = band_pass


    lambdas = IO.load_lambdaslit(lname, maskname, band, options)
    edges, meta = IO.load_edges(maskname, band, options)

    shifts = []
    for pos in nod_posns:
        II = IO.read_drpfits(maskname, "cnts_{0}_{1}.fits".format(band, pos),
                options)
        off = np.array((II[0]["decoff"], II[0]["raoff"]))

        try: off0
        except: off0 = off

        shift = np.sqrt((off[0]-off0[0])**2 + (off[1]-off0[1])**2) * \
                180./np.pi * 3600.0 

        if shift > 3: shift = 3.0
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
    p = Pool()
    solutions = p.map(handle_rectification_helper, xrange(len(edges)))
    p.close()
    tick = time.time()
    print "-----> Mask took %i. Writing to disk." % (tick-tock)

    for i in xrange(len(solutions)):
        solution = solutions[i]
        header = EPS[0].copy()
        
        header.update("OBJECT", "{0}".format(solution["Target_Name"]))

        ll = solution["lambda"]
        ff = solution["frac"]

        header.update("wat0_001", "system=world")
        header.update("wat1_001", "wtype=linear")
        header.update("wat2_001", "wtype=linear")
        header.update("dispaxis", 2)
        header.update("dclog1", "Transform")
        header.update("dc-flag", 0)
        header.update("ctype1", "LINEAR")
        header.update("ctype2", "LINEAR")
        header.update("crval1", ll[0])
        header.update("crval2", 0)
        header.update("cdelt1", ll[1]-ll[0])
        header.update("cdelt2", ff[1]-ff[0])
        header.update("crpix1", 0)
        header.update("crpix2", 0)
        header.update("radecsys", "")
        header.rename_key("cd1_1", "ol_cd1_1")
        header.rename_key("cd1_2", "ol_cd1_2")
        header.rename_key("cd2_1", "ol_cd2_1")
        header.rename_key("cd2_2", "ol_cd2_2")

        IO.writefits(solution["eps_img"], maskname,
                "eps_{0}_S{1:02g}.fits".format(band, i+1), options,
                overwrite=True, header=header)

        IO.writefits(solution["iv_img"], maskname,
                "ivar_{0}_S{1:02g}.fits".format(band, i+1), options,
                overwrite=True, header=header)


def r_interpol(ls, fs, ss, lfid, ffid):
    '''Interpolate the data ss(ls, fs) onto a grid that looks like 
    
    ^ ffid
    |
    o---> lfid
    '''

    S = ss.shape

    output = np.zeros((len(ffid), len(lfid)))

    for i in xrange(len(ffid)):
        ll = ls[i,:] ; sp = ss[i,:]
        ok = np.where(ll>1000)[0]

        if len(ok) < 100: continue

        f = II.interp1d(ll[ok], sp[ok], bounds_error=False)
        output[i,:] = f(lfid)
    
    for i in xrange(len(lfid)):
        fr = fs[:,i] ; sp = output[:,i]
        f = II.interp1d(fr, sp, bounds_error=False)
        output[:,i] = f(ffid)

    return output


def handle_rectification_helper(edgeno):
    global edges, dats, ivars, shifts, lambdas, band

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
    dl = np.median(np.diff(lmid))
    hpp = Filters.hpp[band]
    fidl = lmid

    minl = lmid[0] if lmid[0]>hpp[0] else hpp[0]
    maxl = lmid[-1] if lmid[-1]<hpp[1] else hpp[1]

    fidl = fidl[np.where((fidl > minl) & (fidl < maxl))]
    fun = np.poly1d(np.polyfit(np.arange(len(fidl)), fidl, 1))
    fidl = fun(pix)

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

