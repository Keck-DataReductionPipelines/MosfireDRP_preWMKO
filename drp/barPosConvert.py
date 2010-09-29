# Convert B01POS etc. into the required parameters for the DRP
# Robert Lasenby 2010

from drpUtils import *

from mosfireRaytraceModel import *

import pyfits

import os
import sys
import getopt

import pickle

# values given are in millimetres in focal plane
def getBarPosValues (hdulist):
	names = ["B%02iPOS" % (i+1) for i in xrange(92)]
	values = [hdulist[0].header[s] for s in names]
	return grouper (2, values)

def slitYPos (inst, slit):
	"""Returns (bottom of slit, top of slit, middle of slit),
	where bottom and top refer to the bottom and top of the lower
	and upper gaps adjacent to the slit"""
	y0 = inst.barPitch * slit
	y1 = y0 + inst.barPitch + inst.barGap
	yC = y0 + 0.5 * (inst.barPitch + inst.barGap)
	return (y0/inst.fieldAngle - 1.0, y1/inst.fieldAngle - 1.0, yC/inst.fieldAngle - 1.0)

def computeStuff (inst, n, (b0, b1)):
	slitWidth = (b1-b0)/inst.focalSurfaceScale # in arcsec
	# Right - how to get the slit length?
	slitLength = inst.barPitch - inst.barGap # nominal slit length in arcsec
	# Check that no rigid transformation is needed below
	slitX = 0.5*(b0+b1)/(inst.focalSurfaceScale*inst.fieldAngle) - 1.0 # field angle
	_, _, slitY = slitYPos (inst, n)
	return (n, slitWidth, slitLength, slitX, slitY)

def allStuff (inst, bp):
	return [computeStuff (inst, n, b) for n, b in zip (xrange(46), bp)]

def dumpPickle (hdulist, fname):
	bp = getBarPosValues (hdulist)
	inst = mosfire
	stuff = allStuff (inst, bp)
	d = [[stuff[j][i] for j in xrange(46)] for i in xrange(5)]
	f = open (fname, "w")
	pickle.dump (d, f)
	f.close ()

def main ():
	opts, args = getopt.getopt (sys.argv[1:],"f:")
	imFile = None
	for o,a in opts:
		if o == "-f":
			imFile = a
	# Right - ideally we'd set the band via all this as well ...
	hdulist = pyfits.open (imFile)
	dumpPickle (hdulist, imFile + ".config")
	return

if __name__ == "__main__":
	main()
