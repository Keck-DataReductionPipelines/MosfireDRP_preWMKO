
'''
    npk 23 apr 2012

'''


import os, time

import MOSFIRE
from MOSFIRE import Flats, Options, IO, Wavelength
import numpy as np, pylab as pl, pyfits

fs = ['m120406_0291.fits']

maskname = 'NGC5053'
options = Options.wavelength

path = os.path.join(options["outdir"], maskname)
if not os.path.exists(path):
    raise Exception("Output directory '%s' does not exist. This " 
            "directory should exist." % path)

if False:
    for fname in fs:
        fp = os.path.join(path, fname)

        mfits = IO.readmosfits(fp)
        header, data, bs = mfits

        Wavelength.plot_mask_solution_ds9(fname, maskname, options)
        Wavelength.fit_lambda(mfits, fname, maskname, options)
        Wavelength.apply_lambda(mfits, fname, maskname, options)
        Wavelength.plot_data_quality(maskname, fname, options)
        Wavelength.plot_sky_spectra(maskname, fname, options)
        Wavelength.plot_mask_fits(maskname, fname, options)

if True:
    for fname in fs:
        pass



