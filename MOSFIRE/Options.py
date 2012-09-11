'''

===================
MOSFIRE Options
===================




npk April 18th 2011

'''

import getpass
import os

__version__ = 0.1


npix = 2048

indir = "/scr2/mosfire"
outdir = "/scr2/{0}/mosfire_redux".format(getpass.getuser())
path_bpm = "/scr2/mosfire/badpixels/badpix_10sep2012.fits"

flat = {
        "indir": indir,
        "outdir": outdir,
        "version": 1, 
        "edge-order": 4, # Polynomial order for edge of slit 
        "edge-fit-width": 20,
        "flat-field-order": 7 # Order of polynomial for fitting the flat field profile
}


wavelength = {
        "indir": indir,
        "outdir": outdir,
        "datadir" : os.path.join(os.environ["MOSPATH"], "code", "data"),
        "version": 1,
        #"fractional-wavelength-search": 0.9988, # used in determining oned wavelength solutions
        "fractional-wavelength-search": 0.99935, # used in determining oned wavelength solutions
        "chebyshev-degree": 5, # polynomial order for fitting wavelengths

}


