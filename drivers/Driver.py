# Help, bugs to: http://mosfire.googlecode.com

import os, time
import MOSFIRE

from MOSFIRE import Background, Combine, Detector, Flats, IO, Options, \
    Rectify
from MOSFIRE import Wavelength

import numpy as np, pylab as pl, pyfits as pf

np.seterr(all="ignore")

maskname = ''
band = ''

flatops = Options.flat
waveops = Options.wavelength

#Flats.handle_flats('Flat.txt', maskname, band, flatops)
#Wavelength.imcombine('Offset_1.5.txt', maskname, band, waveops)
#Wavelength.fit_lambda_interactively(maskname, band, 'Offset_1.5.txt', 
    #waveops)
#Wavelength.fit_lambda(maskname, band, 'Offset_1.5.txt', 'Offset_1.5.txt',
    #waveops)

#Wavelength.apply_lambda_simple(maskname, band, 'Offset_1.5.txt', waveops)
#Background.handle_background(['Offset_1.5.txt', 'Offset_-1.5.txt',
    #'Offset_1.2.txt', 'Offset_-1.2.txt'],
    #'fill in after apply_lambda_simple'
    #maskname, band, waveops)

#Rectify.handle_rectification(maskname, ['eps_Offset_1.5.txt.fits',
#    'eps_Offset_-1.5.txt.fits','eps_Offset_1.2.txt.fits',
#    'eps_Offset_-1.2.txt.fits',],
#    'fill in after handle_background step'
#    band, 
#    "/path/to/first/file in Offset_1.5.txt"
#    waveops)
#
