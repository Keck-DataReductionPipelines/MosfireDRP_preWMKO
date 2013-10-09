'''







===================



MOSFIRE Options



===================



















npk April 18th 2011







'''







import getpass



import os







__version__ = 2.0











npix = 2048







indir = '/Users/npk/mosdrp/data'
outdir = '/Users/npk/mosdrp/output'
path_bpm = '/Users/npk/mosdrp/badpixels/badpix_10sep2012.fits'




indir = '/Users/npk/mosdrp/data'




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



        "version": 2,



        #"fractional-wavelength-search": 0.9988, # used in determining oned wavelength solutions



        "fractional-wavelength-search": 0.99935, # used in determining oned wavelength solutions



        "chebyshev-degree": 5, # polynomial order for fitting wavelengths







}











