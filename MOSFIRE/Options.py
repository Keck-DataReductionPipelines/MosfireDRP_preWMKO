'''

===================
MOSFIRE Options
===================




npk April 18th 2011

'''


import MOSFIRE

npix = 2048

flat = {"outdir": "/users/npk/desktop/c9_reduce/", 
        "indir": "/users/npk/desktop/c9/", 
        "version": 1, 
        "edge-order": 5, # Polynomial order for edge of slit 
        "edge-fit-width": 15,
        "flat-field-order": 7 # Order of polynomial for fitting the flat field profile
}
