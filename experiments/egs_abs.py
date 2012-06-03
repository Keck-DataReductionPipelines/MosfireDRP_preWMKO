
'''
    npk 7 may 2012

'''


import os, time

import MOSFIRE
import warnings

from MOSFIRE import Flats, Options, IO, Wavelength, Background, Detector, Rectify 
import numpy as np, pylab as pl, pyfits as pf


reload(Flats)
reload(IO)
reload(Wavelength)
reload(Background)

# Check these
flatops = Options.flat
flatops["outdir"] = "/scr2/mosfire/secondlight/"
flatops["indir"] = "/users/npk/desktop/"
wavlops = Options.wavelength
wavlops["outdir"] = flatops["outdir"]
wavlops["indir"] = flatops["indir"]

maskname = 'egs_abs'
band = 'J'

# Flat Names
flatnames = []
for flat in range(14,24): flatnames.append("m120507_%4.4i.fits" % flat)

# Use the file "lname" for wavelength calibration
lname = "m120507_0230"


# Create A/B positions
As = ["m120507_%4.4i.fits" % i for i in range(229,249,2)]
Bs = ["m120507_%4.4i.fits" % i for i in range(230,249,2)]
As.extend(["m120509_%4.4i.fits" % i for i in range(315,369,2)])
Bs.extend(["m120509_%4.4i.fits" % i for i in range(316,369,2)])

if False:
    Flats.handle_flats(flatnames, maskname, band, flatops)

if False:
    fs = [lname + ".fits"]
    for fname in fs:
        mfits = IO.readmosfits(fname, wavlops)
        header, data, bs = mfits

        #Wavelength.fit_lambda_interactively(mfits, fname, maskname, options)
        #Wavelength.fit_lambda(mfits, fname, maskname, options)
        Wavelength.apply_lambda_simple(mfits, fname, maskname, wavlops)

if True:
    Background.handle_background(As, Bs, lname, maskname, band, wavlops)

if True:
    Rectify.handle_rectification(maskname, ["A", "B"], lname, band, wavlops)


