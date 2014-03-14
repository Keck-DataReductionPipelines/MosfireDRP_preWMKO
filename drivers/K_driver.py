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

obsfiles = ['Offset_1.5.txt', 'Offset_-1.5.txt']
#Wavelength.imcombine(obsfiles, maskname, band, waveops)
#Wavelength.imcombine('Ne.txt', maskname, band, waveops)

#Wavelength.fit_lambda_interactively(maskname, band, obsfiles, 
    #waveops)

#Wavelength.apply_interactive(maskname, band, waveops, 
    #apply=obsfiles, to='Ne.txt', neon=True)

#Wavelength.fit_lambda(maskname, band, obsfiles, obsfiles,
    #waveops)
#Wavelength.fit_lambda(maskname, band, 'Ne.txt', 'Ne.txt',
    #waveops)
#LROI = [[21000, 22800]] * 32
#LROIs = Wavelength.check_wavelength_roi(maskname, band, obsfiles,
    #'Ne.txt', LROI, waveops)

#Wavelength.apply_lambda_simple(maskname, band, obsfiles, waveops)
#Wavelength.apply_lambda_simple(maskname, band, 'Ne.txt', waveops)
#Wavelength.apply_lambda_sky_and_arc(maskname, band, obsfiles,
    #'Ne.txt', LROIs, waveops, neon=True)
#Background.handle_background(obsfiles,
#    'Fill in with appropriate merged_lambda_solution_..fits
#    maskname,
#    band,
#    waveops)

#redfiles = ["eps_" + file + ".fits" for file in obsfiles]
#Rectify.handle_rectification(maskname, ['eps_Offset_1.5.txt.fits', 'eps_Offset_-1.5.txt.fits'],
#    'Fill in with appropriate merged_lambda_solution_..fits
#   band, 
#   "Copy the first file in Offset_1.5.txt here"
#   waveops)
#

