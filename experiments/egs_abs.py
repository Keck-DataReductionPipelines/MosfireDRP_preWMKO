
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


maskname = 'egs_abs'
band = 'J'


if False:
    flatlist = range(14,24)
    fs = []
    for flat in flatlist:
        fs.append("m120507_%4.4i.fits" % flat)

    options = Options.flat
    options["outdir"] = "/scr2/mosfire/secondlight/"
    options["indir"] = "/users/npk/desktop/"
    Flats.handle_flats(fs, maskname, band,  options)


options = Options.wavelength
options["outdir"] = "/scr2/mosfire/secondlight/"
options["indir"] = "/users/npk/desktop/"


lname = "m120507_0230"
fs = [lname + ".fits"]
np.set_printoptions(precision=2)
if False:
    for fname in fs:
        mfits = IO.readmosfits(fname, options)
        header, data, bs = mfits

        #Wavelength.fit_lambda_interactively(mfits, fname, maskname, options)
        #Wavelength.fit_lambda(mfits, fname, maskname, options)
        Wavelength.apply_lambda_simple(mfits, fname, maskname, options)

if True:

    As = ["m120507_%4.4i.fits" % i for i in range(229,249,2)]
    Bs = ["m120507_%4.4i.fits" % i for i in range(230,249,2)]

    As.extend(["m120509_%4.4i.fits" % i for i in range(315,369,2)])
    Bs.extend(["m120509_%4.4i.fits" % i for i in range(316,369,2)])

    Background.handle_background(As, Bs, lname, maskname,
            band, options)


if True:
    Rectify.handle_rectification(maskname, ["A", "B"], lname, band, options)


