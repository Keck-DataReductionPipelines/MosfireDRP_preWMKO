
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
from MOSFIRE import Background, CSU, Fit, IO, Options, Filters, Detector, Wavelength



def handle_rectification(maskname, in_files, wavename, band_pass, barset_file, options,
        commissioning_shift=3.0):
    global edges, dats, stds, itimes, shifts, lambdas, band, fidl, all_shifts


    '''
    '''

    band = band_pass

    
    dlambda = Wavelength.grating_results(band)

    hpp = Filters.hpp[band]
    fidl = np.arange(hpp[0], hpp[1], dlambda)

    lambdas = IO.readfits(wavename, options)
    edges, meta = IO.load_edges(maskname, band, options)
    shifts = []

    posnames = []
    
    for file in in_files:

        II = IO.read_drpfits(maskname, file, options)

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
        posnames.append(II[0]["frameid"])
    
        print "Position {0} shift: {1:2.2f} as".format(off, shift)
    

    plans = Background.guess_plan_from_positions(set(posnames))

    theBPM = IO.badpixelmask()

    all_solutions = []
    all_shifts = [[0.3,3], [0,2.7]]
    cntr = 0
    for plan in plans:
        p0 = plan[0].replace("'", "p")
        p1 = plan[1].replace("'", "p")
        suffix = "%s-%s" % (p0,p1)
        print "Handling plan %s" % suffix
        fname = "bsub_{0}_{1}_{2}.fits".format(maskname,band,suffix)

        EPS = IO.read_drpfits(maskname, fname, options)
        EPS[1] = np.ma.masked_array(EPS[1], theBPM, fill_value=0)

        fname = "std_{0}_{1}_{2}.fits".format(maskname, band, suffix)
        STD = IO.read_drpfits(maskname, fname, options)
        STD[1] = np.ma.masked_array(STD[1], theBPM, fill_value=np.inf)

        fname = "itime_{0}_{1}_{2}.fits".format(maskname, band, suffix)
        ITIME = IO.read_drpfits(maskname, fname, options)
        ITIME[1] = np.ma.masked_array(ITIME[1], theBPM, fill_value=0)


        dats = EPS
        stds = STD
        itimes = ITIME

        EPS[0].update("ORIGFILE", fname)

        tock = time.time()
        sols = range(len(edges)-1,-1,-1)

        shifts = all_shifts[cntr]
        cntr += 1
        p = Pool()
        solutions = p.map(handle_rectification_helper, sols)
        p.close()

        all_solutions.append(solutions)

    tick = time.time()
    print "-----> Mask took %i. Writing to disk." % (tick-tock)


    output = np.zeros((1, len(fidl)))
    snrs = np.zeros((1, len(fidl)))
    sdout= np.zeros((1, len(fidl)))
    itout= np.zeros((1, len(fidl)))


    # the barset [bs] is used for determining object position
    x, x, bs = IO.readmosfits(barset_file, options)
    

    for i_slit in xrange(len(solutions)):
        solution = all_solutions[0][i_slit]
        header = EPS[0].copy()

        target_name = bs.ssl[-(i_slit+1)]['Target_Name']
        pixel_dist = np.float(bs.ssl[-(i_slit+1)]['Target_to_center_of_slit_distance'])/0.18

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


        S = output.shape

        img = solution["eps_img"]
        std = solution["sd_img"]
        tms = solution["itime_img"]


        for i_solution in xrange(1,len(all_solutions)):
            print "Combining solution %i" %i_solution
            solution = all_solutions[i_solution][i_slit]
            img += solution["eps_img"]
            std += solution["sd_img"]
            tms += solution["itime_img"]

        output = np.append(output, img, 0)
        output = np.append(output, np.nan*np.zeros((3,S[1])), 0)
        snrs = np.append(snrs, img*tms/std, 0)
        snrs = np.append(snrs, np.nan*np.zeros((3,S[1])), 0)
        sdout = np.append(sdout, std, 0)
        sdout = np.append(sdout, np.nan*np.zeros((3,S[1])), 0)
        itout = np.append(itout, tms, 0)
        itout = np.append(itout, np.nan*np.zeros((3,S[1])), 0)

        IO.writefits(img, maskname,
            "eps_{0}_S{1:02g}_S{2}.fits".format(band, i_slit+1,target_name), options,
            overwrite=True, header=header, lossy_compress=False)

        IO.writefits(std, maskname,
            "sd_{0}_S{1:02g}_S{2}.fits".format(band, i_slit+1,target_name), options,
            overwrite=True, header=header, lossy_compress=False)

        IO.writefits(tms, maskname,
            "itime_{0}_S{1:02g}_S{2}.fits".format(band, i_slit+1, target_name), options,
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

    IO.writefits(output, maskname, "eps_{0}_{1}.fits".format(maskname,
        band), options, overwrite=True, header=header,
        lossy_compress=False)

    IO.writefits(snrs, maskname, "snrs_{0}_{1}.fits".format(maskname,
        band), options, overwrite=True, header=header,
        lossy_compress=False)

    IO.writefits(sdout, maskname, "sd_{0}_{1}.fits".format(maskname,
        band), options, overwrite=True, header=header,
        lossy_compress=False)

    IO.writefits(itout, maskname, "itime_{0}_{1}.fits".format(maskname,
        band), options, overwrite=True, header=header,
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
    global edges, dats, stds, itimes, shifts, lambdas, band, fidl,all_shifts

    pix = np.arange(2048)
    
    edge = edges[edgeno]

    print "Handling edge: ", edge["Target_Name"]

    tops = edge["top"](pix)
    bots = edge["bottom"](pix)

    # Length of the slit in arcsecond
    lenas = (tops[1024] - bots[1024]) * 0.18
    mxshift = np.abs(np.int(np.ceil(np.max(all_shifts)/0.18)))
    mnshift = np.abs(np.int(np.floor(np.min(all_shifts)/0.18)))

    top = min(np.floor(np.min(tops)), 2048)
    bot = max(np.ceil(np.max(bots)), 0)

    ll = lambdas[1].data[bot:top, :]
    eps = dats[1][bot:top, :].filled(0.0)
    sds = stds[1][bot:top, :].filled(np.inf)
    it  = itimes[1][bot:top, :].filled(0.0)

    lmid = ll[ll.shape[0]/2,:]
    hpp = Filters.hpp[band]

    minl = lmid[0] if lmid[0]>hpp[0] else hpp[0]
    maxl = lmid[-1] if lmid[-1]<hpp[1] else hpp[1]

    epss = []
    ivss = []
    itss = []
    sign = -1
    for shift in shifts:
        output = r_interpol(ll, eps, fidl, tops, top, shift_pix=shift/0.18,
            pad=[mnshift, mxshift])
        epss.append(sign * output)

        var = (sds)**2
        ivar = 1/var
        bad = np.where(np.isfinite(ivar) ==0)
        ivar[bad] = 0.0
        output = r_interpol(ll, ivar, fidl, tops, top, shift_pix=shift/0.18,
            pad=[mnshift, mxshift], fill_value=0.0) 
        ivss.append(output)

        output = r_interpol(ll, it, fidl, tops, top, shift_pix=shift/0.18,
            pad=[mnshift, mxshift], fill_value=0.0) 
        itss.append(output)

        sign *= -1

    eps_img = np.sum(epss, axis=0)
    it_img = np.sum(np.array(itss), axis=0)


    # Remove any NaNs or infs from the variance array
    ivar_img = []
    for ivar in ivss:
        bad = np.where(np.isfinite(ivar) == 0)
        ivar[bad] = 0.0

        ivar_img.append(ivar)
    ivar_img = np.sum(np.array(ivar_img), axis=0)
    sd_img = 1/np.sqrt(ivar_img)

    return {"eps_img": eps_img, "sd_img": sd_img, "itime_img": it_img, 
            "lambda": fidl, "Target_Name": edge["Target_Name"], 
            "slitno": edgeno+1}

