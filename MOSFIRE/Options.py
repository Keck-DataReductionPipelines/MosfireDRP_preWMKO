'''

===================
MOSFIRE Options
===================




npk April 18th 2011

'''


__version__ = "19_04_2011"

import MOSFIRE

npix = 2048

flat = {"outdir": "/scr2/mosfire/c9_npk",
        "indir": "/scr2/mosfire/110326",
        "version": 1, 
        "edge-order": 5, # Polynomial order for edge of slit 
        "edge-fit-width": 15,
        "flat-field-order": 7 # Order of polynomial for fitting the flat field profile
}
