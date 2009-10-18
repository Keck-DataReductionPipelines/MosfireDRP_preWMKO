# Wrapper program performing calibration on a specific slit configuration
# Robert Lasenby 2009

from mosfireRaytraceModel import *
from flatFieldEdgeTracing import *
from coarseAbsoluteWavelengthCalibration import *
from absoluteWavelengthCalibration import *
from relativeWavelengthCalibration import *

import pyfits

import os
import sys
import getopt

import pdb

def saveAsFits (arr, fname):
	if os.path.exists (fname):
		os.remove (fname)
	hdu = pyfits.PrimaryHDU (arr)
	hdu.writeto (fname)

helpstring = """Wrapper program performing calibration on a specific slit configuration
Options:
-d [dir]  : set output directory to [dir]
"""
# utility function to write ds9 reg file
def writeRegLines (f, xA, yA, colour):
	for i in xrange(len(xA)-1):
		# Off-by-one correction needed for ds9 pixel coords
		px0 = xA[i]*1024 + 1024 + 1
		px1 = xA[i+1]*1024 + 1024 + 1
		py0 = yA[i]*1024 + 1024 + 1
		py1 = yA[i+1]*1024 + 1024 + 1
		f.write ("line " + str(px0) + " " + str(py0) + " " + str(px1) + " " + str(py1) + " # color=" + colour + "\n")

def drawSlitEdges (inst, f, slit, pyF):
	# TODO - refactor this pattern?
	pbX = np.linspace (-1.0, 1.0, 50)
	y0, y1 = inst.slitEdges (slit)
	pbY = pyF (pbX, np.repeat (y0, len(pbX)))
	ptY = pyF (pbX, np.repeat (y1, len(pbX)))
	writeRegLines (f, pbX, pbY, "green")
	writeRegLines (f, pbX, ptY, "red")

# Okay - what do want to do as regard diagnostics?
# Need to output a set of sufficiently pretty pictures so as to be convincing ...
# Well, obvious things we need to do:
# - probably want to display the results of the edge tracing in some sensible way ...
# 	- make it so that we can output a .reg file consisting of a *sensible* # of points
# - produce p_x, y image to demonstrate that our p_x, p_y -> y map is in order
# - produce \eta, y image to show that our p_x, p_y -> \eta map has worked
#	- produce a pixel values -> 1D thing here, allowing for closer inspection
#		- basically, just dump stuff out to a file, from where it can be plotted if necessary
# - output a basic pixel values -> 1D lambda thing here, to allow for easy comparison with
# 	model spectrum
# - output a "rectified 1D" reduction, i.e. a sample of the actual "end-product"
#	- think about this once we're sorted out the rest of the stuff ...
#
# Should also do target extraction? Hmmmm ... would need background
# subtraction for that ...

def main ():
	opts, args = getopt.getopt (sys.argv[1:],"hd:")
	#if opts == [] or opts[0][0] == "-h":
	if opts != [] and opts[0][0] == "-h":
		print helpstring
		return
	outputdir = "calibOutput1"
	for o,a in opts:
		if o == "-d":
			outputdir = a
	if not os.path.exists (outputdir):
		os.mkdir (outputdir)
	# Okay - so, we need to modify the simulator so that we're including
	# the data as to band, exposure time and all that in the headers - probably
	# the course of least resistance for the moment ... (?)
	# Also need to sort out the 0.5 pixel mismatch ...
	#
	# Will be used for edge tracing ...
	ffFile = '/home/robert/Projects/mosfire/drp/mosfire/spectroscopicSimulator/session5/flatField.fits'
	#
	sciFile = '/home/robert/Projects/mosfire/drp/mosfire/spectroscopicSimulator/session5/0.fits'
	ffHduList = pyfits.open (ffFile)
	sciHduList = pyfits.open (sciFile)
	slitConfig = slitConfigFromFits (ffHduList)
	qe0 = pyfits.getdata (detectorQEFile)
	# TODO - what do the -ve qe values mean?
	qe = np.clip (qe0, 0.0, 1.0)
	#
	# TODO - should get this stuff from fits ..
	bandName = "K"
	inst = mosfire
	optData0 = loadRaytraceData (raytraceFile)
	optData1 = cleanOpticalData (optData0)
	processOpticalData (inst, optData1)
	band = bandFromRaytrace (inst, bandName, optData1)
	getTransferFn (band)
	# Which slits do we want to deal with?
	# TODO - take this as input
	slitNums = [35]
	# Flat field edge tracing, i.e.
	# inferring p_x, y -> p_y
	# TODO - would be nice to be able to only process certain of the slits ...
	# make this happen?
	slits = [slitConfig[n] for n in slitNums]
	pyFA = calibFlat (inst, band, ffHduList, qe, shape1, slits)
	# TODO - take this stuff as input
	outputSlitEdges = True 
	# should probably put all the stuff in a suitable folder ...
	slitRegFile = os.path.join (outputdir, "slitCalibTest.reg")
	if outputSlitEdges:
		freg = open (slitRegFile, "w")
		for slit, pyF in zip (slits, pyFA):
			drawSlitEdges (inst, freg, slit, pyF)
		freg.close ()
	#return
	# Preparing the sky image for use (masking out object, flat-fielding, adding in Poisson-variance stuff, etc.)
	# Okay ... 
	iSky, wSky = qeCompensate (sciHduList[0].data, sciHduList[1].data, qe)
	# Load up the sky model ...
	skyBg1 = loadSkyBg (skyBgFile)
	lASky, iASky = zip(*[(l, i) for (l, i) in skyBg1 if (l >= band.minL and l <= band.maxL)])
	lASky = np.array (lASky)
	iASky = np.array (iASky)
	lineListFile = "list_v2.0.dat"
	lineList = loadLineList (lineListFile)
	outputIntermediateImages = True
	# Now, do things a slit at a time from here ...
	for n in xrange(len(slits)):
		slit = slits[n]
		# Relative wavelength calibration, i.e.
		# inferring p_x, p_y -> \eta
		# Have h : \eta -> \lambda
		# TODO - want to pass pyF in here ...
		pyF = pyFA[n]
		py0, slitIm, h, hInv = relativeWavelengthCalibrateSlit (inst, band, slit, iSky, wSky, pyF)
		# Coarse absolute wavelength calibration, i.e.
		# inferring kAb = \kappa (nAb)
		nA0, kA0, kStep = coarseAbsoluteWavelengthCalibrateSlit (inst, band, slitIm, iASky, lASky, h, hInv)
		# Fine absolute wavelength calibration, i.e.
		# inferring g : \kappa -> \eta
		# Okay - we need to pass throught the results of the coarse calibration here ...
		g = absoluteWavelengthCalibrateSlit (inst, band, slitIm, hInv, lineList, nA0, kA0, kStep)
		# Okay ... what to do now?
		iAA, wAA, nAA, cond = slitIm
		nMin = np.min (nAA)
		nMax = np.max (nAA)
		slitName = os.path.join (outputdir, "slit" + str(slitNums[n]))
		if outputIntermediateImages:
			y0, y1 = inst.slitEdges (slit)
			pxA = inst.px2pxn (np.arange (inst.nPx))
			yC = np.linspace (y0, y1, 40)
			pyAA = np.transpose([inst.pyn2py (pyF (np.ones_like(yC)*px, yC)) for px in pxA])
			pyAA -= py0
			pxCoords = np.tile (inst.pxn2px (pxA),  (len(yC),1))
			coords = np.array ([pyAA, pxCoords])
			pxyIm = scipy.ndimage.interpolation.map_coordinates (iAA, coords)
			saveAsFits (pxyIm, slitName + "pxy.fits")
			# Okay ... now the \eta, y image, I suppose ...
			# make sure that stuff behaves properly at the edges, I spose ...
			nIm = scipy.ndimage.interpolation.map_coordinates (nAA, coords)
			nAAc = np.arange (nMin, nMax, 2.0/float(inst.nPx))
			nyIm = np.array ([ resample1D (nIm[i], pxyIm[i], nAAc) for i in xrange(len(nIm))])
			saveAsFits (nyIm, slitName + "ny.fits")
		# How about a comparison with the model spectrum?
		# Definitely need to allow this in some way ... maybe it's best to just dump the data,
		# leave it at that? Okay - do that ftm ...
		# Easiest thing at first, I spose, is to dump something like ...
		# iAA, wAA, nAA, lAA, cond
		# Then can build a separate little visualisation thing ...
		saveAsFits (iAA, slitName + "iAA.fits")
		saveAsFits (wAA, slitName + "wAA.fits")
		saveAsFits (nAA, slitName + "nAA.fits")
		# FIXME - pyfits apparently doesn't like saving boolean images
		#saveAsFits (cond, slitName + "cond.fits")
		saveAsFits (cond.astype(np.int), slitName + "cond.fits")
		# Now - how to deal with the lambda stuff, exactly ... ?
		# TODO - find k0 st g(k0) = nMin - epsilon etc.
		#k0 = scipy.optimise.brentq (lambda k : g(k) - (nMin - 0.01), nMin-0.1, nMin+0.1)
		#k1 = scipy.optimise.brentq (lambda k : g(k) - (nMax + 0.01), nMax - 0.1, nMax + 0.1)
		# FIXME - magic numbers
		k0 = nMin - 0.1
		k1 = nMax + 0.1
		kA1 = np.linspace (k0, k1, 50)
		nA1 = g (kA1)
		gInv = scipy.interpolate.InterpolatedUnivariateSpline (nA1, kA1)
		# TODO - should make sure that everything is in the right range, etc. etc.
		#lAA = h (gInv (nAA))
		lAA = np.apply_along_axis (h, 1, np.apply_along_axis (gInv, 1, nAA))
		#lAA = np.apply_along_axis (h, 1, nAA)
		saveAsFits (lAA, slitName + "lAA.fits")

if __name__ == "__main__":
	main()
