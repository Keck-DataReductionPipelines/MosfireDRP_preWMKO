
import numpy as np
import pyfits
import sys, os

fname = sys.argv[1]

# Right - just put in 

hdulist = pyfits.open (fname)

im = hdulist[0].data

# Flip the image in the x axis ...

# Flat variance image
varIm = np.ones_like (im)

imNew = np.fliplr (im)

imHdu = pyfits.PrimaryHDU (imNew)
varHdu = pyfits.ImageHDU (varIm)

newHduList = pyfits.HDUList ([imHdu, varHdu] + hdulist[1:])

newname = fname + ".new.fits"

if os.path.exists (newname):
	os.remove (newname)

newHduList.writeto (newname)
