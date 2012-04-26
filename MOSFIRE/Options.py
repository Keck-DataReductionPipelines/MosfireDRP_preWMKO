'''

===================
MOSFIRE Options
===================




npk April 18th 2011

'''


__version__ = "19_04_2011"


npix = 2048

flat = {"outdir": "/scr2/mosfire/firstlight/",
        "indir": "/scr2/mosfire/firstlight/",
        "version": 1, 
        "edge-order": 5, # Polynomial order for edge of slit 
        "edge-fit-width": 20,
        "flat-field-order": 7 # Order of polynomial for fitting the flat field profile
}


wavelength = {"outdir": flat["outdir"],
        "indir": flat["indir"],
        "version": 1,
        "fractional-wavelength-search": 0.9988, # used in determining oned wavelength solutions

}
