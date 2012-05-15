
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


if False:
    flatlist = range(14,24)
    fs = []
    for flat in flatlist:
        fn = os.path.join("/users/npk/desktop/7may2012/m120507_%4.4i.fits" % flat)
        fs.append(fn)

    options = Options.flat
    options["outdir"] = "/scr2/mosfire/secondlight/"
    Flats.handle_flats(fs, maskname, band,  options)


options = Options.wavelength
options["outdir"] = "/scr2/mosfire/secondlight/"
options["indir"] = "/users/npk/desktop/7may2012/"
path = os.path.join(options["outdir"], maskname)
if not os.path.exists(path):
    raise Exception("Output directory '%s' does not exist. This " 
            "directory should exist." % path)

ind = options['indir']

if False:
    As = [ind + "m120507_%4.4i" % i for i in range(229,249,2)]
    Bs = [ind + "m120507_%4.4i" % i for i in range(230,249,2)]

    IO.imcombine(As, "A", reject='none')
    IO.imcombine(Bs, "B", reject='none')

    header, datA = IO.readfits("A.fits")
    header, datB = IO.readfits("A.fits")
    ffp = os.path.join(options['outdir'], '%s/pixelflat_2d_%s.fits' %
            (maskname, band))
    header, flat = IO.readfits(ffp)

    AmB = (datA/flat) - (datB/flat)
    AmB *= Detector.gain

    outname = os.path.join(options['outdir'], 'AmB.fits') 
    hdu = pf.PrimaryHDU(AmB)
    try: os.remove(outname)
    except: pass
    hdu.writeto(outname)



inpath = options['indir']
fs = ['m120507_0230.fits']
np.set_printoptions(precision=2)
if True:
    for fname in fs:
        fp = os.path.join(inpath, fname)

        mfits = IO.readmosfits(fp)
        header, data, bs = mfits

        #Wavelength.plot_mask_solution_ds9(fname, maskname, options)
        Wavelength.fit_lambda_interactively(mfits, fname, maskname, options)
        #Wavelength.fit_lambda(mfits, fname, maskname, options)
        #Wavelength.apply_lambda(mfits, fname, maskname, options)
        #Wavelength.plot_data_quality(maskname, fname, options)
        #Wavelength.plot_sky_spectra(maskname, fname, options)
        #Wavelength.plot_mask_fits(maskname, fname, options)

if False:

    Background.handle_background(['AmB.fits'], 'm120507_0230', maskname,
            band, "120506_164909_egs_abs.fits", options)


