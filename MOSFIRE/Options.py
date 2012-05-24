'''

===================
MOSFIRE Options
===================




npk April 18th 2011

'''


__version__ = 0.1


npix = 2048

flat = {"outdir": "/scr2/mosfire/secondlight/",
        "indir": "/users/npk/desktop",
        "version": 1, 
        "edge-order": 5, # Polynomial order for edge of slit 
        "edge-fit-width": 20,
        "flat-field-order": 7 # Order of polynomial for fitting the flat field profile
}


wavelength = {"outdir": flat["outdir"],
        "indir": flat["indir"],
        "datadir" : "/Users/npk/Dropbox/MOSFIRE/code/data",
        "version": 1,
        #"fractional-wavelength-search": 0.9988, # used in determining oned wavelength solutions
        "fractional-wavelength-search": 0.99935, # used in determining oned wavelength solutions
        "chebyshev-degree": 5, # polynomial order for fitting wavelengths

}


path_bpm = "/users/npk/desktop/8may/badpix_18may2012.fits"
