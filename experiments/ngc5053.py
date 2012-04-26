
'''
    npk 23 apr 2012

'''


import MOSFIRE
import time
from MOSFIRE import Flats, Options, IO, Wavelength
import numpy as np, pylab as pl, pyfits


reload(Flats)
reload(IO)
reload(Wavelength)

fs = ['m120406_0292.fits']

MOSFIRE.Wavelength.handle_lambdas(fs, 'NGC5053', Options.wavelength)

