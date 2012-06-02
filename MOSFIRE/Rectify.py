
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



def handle_rectification(maskname, nod_posns, lname, band, options):
    global edges, dats, shifts, lambdas


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
    II = IO.read_drpfits(maskname, fname, options)

    dats = II
    II[0].update("ORIGFILE", fname)

    tock = time.time()
    p = Pool()
    solutions = p.map(handle_rectification_helper, xrange(len(edges)))
    p.close()
    tick = time.time()
    print "-----> Mask took %i. Writing to disk." % (tick-tock)

    for i in xrange(len(solutions)):
        solution = solutions[i]
        header = II[0].copy()
        
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
        header.update("cdelt2", 0)
        header.update("radecsys", "")
        header.rename_key("cd1_1", "ol_cd1_1")
        header.rename_key("cd1_2", "ol_cd1_2")
        header.rename_key("cd2_1", "ol_cd2_1")
        header.rename_key("cd2_2", "ol_cd2_2")

        IO.writefits(solution["img"], maskname,
                "stacked_{0}_S{1:02g}.fits".format(band, i+1), options,
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
    global edges, dats, shifts, lambdas

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
    ss = dats[1][bot:top, :]

    fidf = fracs[:,1024]
    lmid = ll[ll.shape[0]/2,:]
    fl = np.poly1d(np.polyfit(pix,lmid,1))
    fidl = fl(pix)


    interps = []
    sign = 1
    for shift in shifts:
        interps.append(sign * r_interpol(ll, fracs, ss, 
            fidl, fidf - shift/lenas))

        sign *= -1

    np.array(interps)
    img = np.sum(interps, axis=0)

    return {"img": img, "lambda": fidl, "Target_Name": edge["Target_Name"],
            "slitno": edgeno+1, "frac": fidf}

