
'''
    npk 7 may 2012

'''


import os, time

import MOSFIRE
from MOSFIRE import Flats, Options, IO, Wavelength, Background, Detector
import numpy as np, pylab as pl, pyfits as pf


reload(Flats)
reload(IO)
reload(Wavelength)
reload(Background)


maskname = 'egs_abs'
band = 'J'


if True:
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


ind = options['indir']
fs = ['m120507_0230.fits']
np.set_printoptions(precision=2)
if False:
    for fname in fs:
        mfits = IO.readmosfits(fname, options)
        header, data, bs = mfits

        #Wavelength.fit_lambda_interactively(mfits, fname, maskname, options)
        #Wavelength.fit_lambda(mfits, fname, maskname, options)
        #Wavelength.apply_lambda_simple(mfits, fname, maskname, options)
        #Wavelength.plot_mask_fits(maskname, fname, options)

if False:

    As = [ind + "m120507_%4.4i.fits" % i for i in range(229,249,2)]
    Bs = [ind + "m120507_%4.4i.fits" % i for i in range(230,249,2)]

    Background.handle_background(As, Bs, 'm120507_0230', maskname,
            band, options)


