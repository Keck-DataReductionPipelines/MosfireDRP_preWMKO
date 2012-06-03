
'''

    MOSFIRE pre production template for slitmask reductions

    Instructions

    1. Check lines under the comment "Check paths" and modify the input and
    output dirs. 

    * The input directory should be a top-level directory on ramekin looks
      something like:

    ramekin% ls -d 2012*
    2012jun02/  2012may05/  2012may06a/  2012may08/
    2012jun03/  2012may06/  2012may07/   2012may09/

    * The output directory is a top level directory where a directory with mask
      names will be put. e.g.

    ramekin% cd /scr2/npk/mosfire_redux/
    /scr2/npk/mosfire_redux
    ramekin% ls
    egs_abs/  egs_abs.py  Q1700.spike/  q1700_spike.py

    * I recommend you copy the template file to the output directory, rename
      the template file to the mask name (one template file per band) and run
      the pipeline from the outdir

    2. Change the variable maskname and band to the mask and band names. e.g.
    maskname = 'egs_abs' and bandname = 'J'

    3. Change the flat names and range. If you took ten flats on 3 Jun 2012,
    then you would use the following line:
    flatnames = ["m120603_%4.4i.fits" % i for i in range(58, 68)]

    4. Change the variable lname to the file you want to perform the wavelength
    calibration on. I generally pick the first file in a set of A/Bs

    5. Change the A/B positions in the same way. If you have files from
    multiple nights, then you will need to use the python "extend" function.
    
    6. Uncomment the step you would like to perform.
        * Each step has an if False, which means it won't be executed.
        * You will have to perform each step sequentially!!
        * Some steps take a long time, one step is interactive, hence the need
          to if True the appropriate step


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

# Check paths
flatops = Options.flat
flatops["outdir"] = "/scr2/mosfire/secondlight/"
flatops["indir"] = "/users/npk/desktop/"
wavlops = Options.wavelength
wavlops["outdir"] = flatops["outdir"]
wavlops["indir"] = flatops["indir"]

maskname = 'CHANGE'
band = 'CHANGE'

# Flat Names
flatnames = ["m120603_%4.4i.fits" % i for i in range(58, 68)]

# Use the file "lname" for wavelength calibration
lname = "m120603_0310"

# Create A/B positions
As = ["m120603_%4.4i.fits" % i for i in range(310,330,2)]
Bs = ["m120603_%4.4i.fits" % i for i in range(311,330,2)]

fs = lname + ".fits"
mfits = IO.readmosfits(fname, wavlops)
header, data, bs = mfits

# Change if False to if True when you want to execute that step
if False: Flats.handle_flats(flatnames, maskname, band, flatops)
if False: Wavelength.fit_lambda_interactively(mfits, fname, maskname, wavlops)
if False: Wavelength.fit_lambda(mfits, fname, maskname, wavlops)
if False: Wavelength.apply_lambda_simple(mfits, fname, maskname, wavlops)
if False: Background.handle_background(As, Bs, lname, maskname, band, wavlops)
if False: Rectify.handle_rectification(maskname, ["A", "B"], lname, band,
        wavlops)


