# Help, bugs to: http://mosfire.googlecode.com

import os, time
import MOSFIRE

from MOSFIRE import Background, Combine, Detector, Flats, IO, Options, \
    Rectify
from MOSFIRE import Wavelength

import numpy as np, pylab as pl, pyfits as pf

np.seterr(all="ignore")

maskname = ''
band = 'K'

flatops = Options.flat
waveops = Options.wavelength

#Flats.handle_flats('Flat.txt', maskname, band, flatops)

#Wavelength.imcombine('Offset_1.5.txt', maskname, band, waveops)
#Wavelength.imcombine('Ne.txt', maskname, band, waveops)

#Wavelength.fit_lambda_interactively(maskname, band, 'Offset_1.5.txt', 
    #waveops)

#Wavelength.apply_interactive(maskname, band, waveops, 
    #apply='Offset_1.5.txt', to='Ne.txt', neon=True)

#Wavelength.fit_lambda(maskname, band, 'Offset_1.5.txt', 'Offset_1.5.txt',
    #waveops)
#Wavelength.fit_lambda(maskname, band, 'Ne.txt', 'Ne.txt',
    #waveops)
#LROI = [[21000, 22800]] * 32
#LROIs = Wavelength.check_wavelength_roi(maskname, band, 'Offset_1.5.txt', 
    #'Ne.txt', LROI, waveops)

#Wavelength.apply_lambda_simple(maskname, band, 'Offset_1.5.txt', waveops)
#Wavelength.apply_lambda_simple(maskname, band, 'Ne.txt', waveops)
#Wavelength.apply_lambda_sky_and_arc(maskname, band, 'Offset_1.5.txt', 
    #'Ne.txt', LROIs, waveops, neon=True)
#Background.handle_background(['Offset_1.5.txt', 'Offset_-1.5.txt'],
#    'fill this in after apply_lambda_sky_and_arc'
#    maskname,
#    band,
#    waveops)
#Rectify.handle_rectification(maskname, ['eps_Offset_1.5.txt.fits', 'eps_Offset_-1.5.txt.fits'],
#    'fill this in after handle_background'
#   band, 
#   "Copy the first file in Offset_1.5.txt here"
#   waveops)
#

