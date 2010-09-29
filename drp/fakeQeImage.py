
import numpy as np
import pyfits
import sys, os

def saveAsFits (arr, fname):
	if os.path.exists (fname):
		os.remove (fname)
	hdu = pyfits.PrimaryHDU (arr)
	hdu.writeto (fname)

fname = sys.argv[1]

arr = np.ones ((2048,2048))

# set first 4 rows and last 4 rows to zero
arr[0:4] = 0.0
arr[2048-5:2048] = 0.0

# set first and last 4 columns to zero
arr[...,0:4] = 0.0
arr[...,2048-5:2048] = 0.0

saveAsFits (arr, fname)
