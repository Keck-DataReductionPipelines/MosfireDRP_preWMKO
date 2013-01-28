
# MOSFIRE Longslit Reductions
# 5 Aug 2012
# Nick Konidaris

import os
import pdb
import numpy as np
import scipy

from MOSFIRE import Detector, IO, Filters, Wavelength

def rectify(path, dname, lamdat, A, B, maskname, band, wavoptions, 
        longoptions):

    header, data = IO.readfits(os.path.join(path, dname))
    raw_img = data * Detector.gain / header['TRUITIME']

    dlam = Wavelength.grating_results(band)
    hpp = np.array(Filters.hpp[band]) 
    ll_fid = np.arange(hpp[0], hpp[1], dlam)

    rectified = np.zeros((2048, len(ll_fid)))

    from scipy.interpolate import interp1d

    for i in xrange(2048):
        ll = lamdat[1][i,:]
        ss = raw_img[i,:]
        ok = np.isfinite(ll) & np.isfinite(ss) & (ll < hpp[1]) & (ll >
                hpp[0])

        if len(np.where(ok)[0]) < 30:
            continue

        f = interp1d(ll[ok], ss[ok], bounds_error=False)
        rectified[i,:] = f(ll_fid)

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
    IO.writefits(rectified, maskname, "rectified_%s" % (dname), 
        wavoptions, header=header, overwrite=True, lossy_compress=True)


def imdiff(A, B, maskname, band, options):
    s = "[0]"

    targname = A[1]["targname"].rstrip(" ")
    if targname == "":
        objname = A[1]["object"].replace(" ", "_")
    else:
        objname = targname.replace(" ", "_")

    outpath = os.path.join(options["outdir"], maskname)

    try: operand1 = os.path.join(IO.fname_to_path(A[0], options), A[0]) + s
    except: operand1 = A[0]

    try: operand2 = os.path.join(IO.fname_to_path(B[0], options), B[0]) + s
    except: operand2 = B[0]

    dname = "{0}_{1}_{2}_{3}-{4}.fits".format(maskname, objname, band,
        A[1]["frameid"], B[1]["frameid"])
    IO.imarith(operand1, '-', operand2, os.path.join(outpath, dname))
    print dname

    return outpath, dname

def go(maskname,
        band,
        filenames,
        wavoptions,
        longoptions):

    wavename = Wavelength.filelist_to_wavename(filenames, band, maskname,
            wavoptions).rstrip(".fits")

    lamdat = IO.load_lambdaslit(wavename, maskname, band, wavoptions)

    print("Wavelength solution {0}".format(wavename))
    print("{0:18s} {1:30s} {2:2s} {3:5s}".format("filename", "object", "pos",
    "offset"))
    positions = []
    objname = None
    for file in filenames:
        header, data, bs = IO.readmosfits(file, wavoptions)

        if objname is None:
            objname = header["object"]

        if objname != header["object"]:
            raise Exception("Trying to combine longslit stack of object {0} " 
                    "with object {1}".format(objname, header["object"]))

        print("{0:18s} {1:30s} {2:2s} {3:4.1f}".format(file, header["object"],
            header["frameid"], header["yoffset"]))

        positions.append([file, header, data, bs])

    print("{0:2g} nod positions found. Producing stacked difference" \
           " image.".format(len(positions)))

    for i in xrange(len(positions)-1):
        A = positions[i]
        B = positions[i+1]

        path, dname = imdiff(A, B, maskname, band, wavoptions)
        rectify(path, dname, lamdat, A, B, maskname, band, wavoptions,
                longoptions)
        print dname

    path, dname = imdiff(B, A, maskname, band, wavoptions)
    print dname
    rectify(path, dname, lamdat, B, A, maskname, band, wavoptions,
            longoptions)
    
    fname = os.path.join(path, wavename + ".fits")
    B = IO.readfits(fname)
    B = [fname, B[0], B[1]]
    for i in xrange(len(positions)):
        A = positions[i]
        imdiff(A, B, maskname, band, wavoptions)
        rectify(path, dname, lamdat, A, B, maskname, band, wavoptions,
            longoptions)
    imdiff(B, A, maskname, band, wavoptions)
    rectify(path, dname, lamdat, B, A, maskname, band, wavoptions,
            longoptions)


