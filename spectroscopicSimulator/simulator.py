# Simulation of MOSFIRE output, spectroscopic mode
# Robert Lasenby 2009

# TODO - more use of docstrings

# List of problems
# -> We do not infer dy / dpy from the data (though this should vary very little)
# -> We do not infer the width of the convolution needed from the data -
# should investigate to determine the severity of this problem
# -> We do not take into account the pixel width when doing the convolution,
# which will cause problems for thin slits

# Should also document ways in which we've hard-coded mosfire behaviour ...

# Also, want a "feature list"
# e.g.
# -> cosmic rays
# -> mascgen (slit mask configuration files)
# -> read noise

from math import *
import numpy as np
import pyfits
from scipy import interpolate
import pyfits

# TODO - refactor
def unzip (l):
	"""Transposes the first two levels of a list of list (resp tuples etc.)"""
	c = zip(*l)
	return [list(a) for a in c]

class Instrument():
	"""Stores the data pertaining to the properties of the instrument"""
	pass

def slitYPos (inst, slit):
	"""Returns (bottom of slit, top of slit, middle of slit),
	where bottom and top refer to the bottom and top of the lower
	and upper gaps adjacent to the slit"""
	y0 = inst.barPitch * slit
	y1 = y0 + inst.barPitch + inst.barGap
	yC = y0 + 0.5 * (inst.barPitch + inst.barGap)
	return (y0/inst.fieldAngle - 1.0, y1/inst.fieldAngle - 1.0, yC/inst.fieldAngle - 1.0)

def slitXInterval (inst, slit, up):
	if slit < 0 or slit >= inst.nBars: return None
	xW = inst.slitWidth[slit]
	xC = inst.slitX[slit]
	#offset = 0.5 * inst.barPitch * inst.sinSlitTilt * up / inst.fieldAngle
	x0 = xC - 0.5*xW + inst.barOffset * up
	x1 = x0 + xW
	return (x0, x1)

def intersectIntervals ((x0,x1),(y0,y1)):
	z0 = max (x0, y0)
	z1 = min (x1, y1)
	if z0 <= z1: return (z0,z1)
	else: return None

def intervalLen ((x,y)) : return y - x

def computeSlitOverlaps (inst):
	def overlaps (slit1, slit2):
		i1 = slitXInterval (inst, slit1, 1) 
		i2 = slitXInterval (inst, slit2, -1)
		if i1 == None or i2 == None:
			return None
		else:
			iC = intersectIntervals (i1, i2)
			if iC == None: return None
			else:
				iLen = intervalLen (iC)
				trans1 = iLen / inst.slitWidth[slit1]
				trans2 = iLen / inst.slitWidth[slit2]
				rtrans1 = trans1 / (trans1 + trans2)
				rtrans2 = trans2 / (trans1 + trans2)
				return ((0.5*trans1, rtrans1), (0.5*trans2, rtrans2))
	inst.slitOverlaps = [overlaps (slit, slit+1) for slit in xrange (-1,inst.nBars)]

# group together slits that share an overlap
def computeSlitGroups (inst):
	groups = []
	currentGroup = [0]
	for i in xrange (inst.nBars):
		if inst.slitOverlaps[i+1] == None:
			groups += [currentGroup]
			currentGroup = [i+1]
		else:
			currentGroup += [i+1]
	inst.slitGroups = groups

def slitYSamp (inst, slit, ySpacing):
	y0 = inst.barPitch * slit
	y0d = y0 + inst.barGap
	y1d = y0 + inst.barPitch
	y1 = y1d + inst.barGap
	yA = np.arange (y0, y1, ySpacing * inst.fieldAngle)
	# note - need epsilon offsets to cope with floating-point error
	lowerUpperMiddle = [yA - y0 - 1e-10 <= inst.barGap, y1 - yA - 1e-10 <= inst.barGap, True]
	# exponential falloff
	eA1 = np.array( [exp (- (((y - y0d) / inst.slitFalloffScale)**2) / 2.0) for y in yA] )
	eA2 = np.array( [exp (- (((y - y1d) / inst.slitFalloffScale)**2) / 2.0) for y in yA] )
	oLower = inst.slitOverlaps[slit]
	oUpper = inst.slitOverlaps[slit+1]
	transLower, rtransLower = oLower[1] if oLower != None else (0.0, 1.0)
	transUpper, rtransUpper = oUpper[0] if oUpper != None else (0.0, 1.0)
	eA = np.select (lowerUpperMiddle, 
			[transLower + (0.5 - transLower)*eA1, 
			transUpper + (0.5 - transUpper)*eA2, 1.0])
	eAx = np.select (lowerUpperMiddle, [rtransLower, rtransUpper, 1.0])
	# need to convert yA to field angle
	return (yA / inst.fieldAngle - 1.0), eA, eAx

class Band:
	pass

def loadSkyBg (fname):
	"""Load the sky background data from a file, returning
	[(um, ph/sec/arcsec^2/nm/m^2)]"""
	f = open(fname, 'r')
	def readLn(l) : return (float(l[0])/1000.0, float(l[1]))
	d = [readLn(l.split()) for l in f]
	f.close()
	return d

def loadTransferFn (fname):
	"""Load the instrument transfer function from a file, returning
	[(um, fractional efficiency)]"""
	f = open(fname, 'r')
	def readLn(l) : return (float(l[0]), float(l[1]))
	d = [readLn(l.split()) for l in f]
	f.close()
	return d

def loadRaytraceData (fname):
	"""Load the raytracing data from a file, returning a record array"""
	f = open (fname, 'r')
	# check that the data is in the correct order
	l1 = f.readline()
	assert l1.split() == ['config','xin','yin','lambda(um)','xfp(mm)','yfp(mm)','xout(mm)','yout(mm)','slit#']
	def readLn(l) : return [l[0]] + [float(n) for n in l[1:8]] + [int(float(l[8]))]
	d = [readLn(l.split()) for l in f if l.strip() != "done"]
	dt = np.dtype([('band', 'a1'), ('xin', 'f8'), ('yin', 'f8'),
		('lambda', 'f8'), ('xfp', 'f8'), ('yfp', 'f8'),
		('px', 'f8'), ('py', 'f8'), ('slit', 'i4')])
	da = np.rec.array (d, dtype = dt)
	f.close()
	return da

###########
# Overview of process:
# -> Produce p_x, y image using p_x, x, y -> lambda map, finely sampled in p_x
# -> Convolve in p_x direction to take into account width of slit and pixels,
# then resample down to detector resolution in p_x
# -> Use p_x, x, y -> p_y map to transform to p_x, p_y, then resample
# down to detector resolution in p_y

# We get input intensity in
# I0 = photons / (sec arcsec^2 um m^2)
# (actually in nm, but we convert to um immediately)
# In the first stage, we convert to
# I1 = photons / (sec arcsec px) = I0 * area * slit width in arcsec * (um/px)
# In the second stage, we convert to
# I2 = photons / (sec px py) = I1 * (arsec/py)

def drawSlits (inst, band):
	im = np.zeros ((inst.nPy, inst.nPx))
	# Aim for ~1 sample per pixel
	ySpacing = inst.barPitch / (ceil (inst.barPitch / inst.pixScale))
	nYslit = int (floor ((inst.barPitch + inst.barGap) / inst.pixScale)) + 1
	ySpacing /= inst.fieldAngle
	for group in inst.slitGroups:
		print group
		nSlits = len(group)
		if nSlits == 1:
			slit = group[0]
			xA, yA, iAA = slitImage (inst, band, slit, ySpacing)
			pyA, iAAd = distortY (inst, band, xA, yA, iAA)
			im[pyA] += iAAd
		else:
			y0, _, _ = slitYPos (inst, group[0])
			_, y1, _ = slitYPos (inst, group[-1])
			nYg = int (round ((y1 - y0) / ySpacing)) + 1
			iAAg = np.zeros ((nYg, inst.nPx))
			yAg = np.linspace (y0, y1, nYg)
			xAg = np.zeros (nYg)
			for slit in group:
				print slit
				y0s, _, _ = slitYPos (inst, slit)
				# xA, iAA already modified for falloff / mixing
				xA, yA, iAA = slitImage (inst, band, slit, ySpacing)
				offset = int (round ((y0s - y0) / ySpacing))
				xAg[offset:offset+nYslit] += xA
				iAAg[offset:offset+nYslit] += iAA
			pyA, iAAd = distortY (inst, band, xAg, yAg, iAAg)
			im[pyA] += iAAd
	return im

# use p_x, x, y -> p_y map to do a series of 1D resamplings
def distortY (inst, band, xA, yA, iAA):
	pxA = np.arange (inst.nPx)
	pyAA = np.squeeze ([band.pyF (pxA, x, y) for (x,y) in zip (xA, yA)])
	# for the moment, use fixed scale rather than inferring dp_y / dy
	iAA *= inst.pixScale
	pyMin = max (int (floor (np.min (pyAA))), 0)
	pyMax = min (int (ceil (np.max (pyAA))), inst.nPy-1)
	pyA = np.arange (pyMin, pyMax+1)
	def resampleCol (c):
		# y and py are reversed wrt each other
		# again, we should really do sorting here,
		# but since we know about the reversal we can just do that
		# (need pyC monotonically increasing for the interpolation)
		pyC = np.flipud (pyAA[...,c])
		iC = np.flipud (iAA[...,c])
		py0 = pyC[0]
		py1 = pyC[-1]
		#iF = interpolate.InterpolatedUnivariateSpline (pyC, iC)
		# go with linear interpolation at the moment to prevent ringing
		# artifacts where we have rapid falloff
		iF = interpolate.InterpolatedUnivariateSpline (pyC, iC, k=1)
		iA = iF (pyA)
		cond = np.logical_and (pyA >= py0, pyA <= py1)
		return np.where (cond, iA, 0.0)
	return pyA, np.transpose ([resampleCol(c) for c in pxA])

# note - slit inaccuracy in that we're not convolving over pixel width
# - we should really add that effect in if we're going to deal with
# narrow slits
# Of course, our convolution strategy isn't really accurate anyway ...
# the anamorphic factor stuff we're using is just an approximation,
# and we could get the "correct" answer from the raytracing data ...
def slitImage (inst, band, slit, ySpacing):
	# transfer function has already been applied
	lA, iA = band.spectralIntensity
	iF = interpolate.InterpolatedUnivariateSpline (lA, iA)
	y0, y1, yC = slitYPos (inst, slit)
	yA, eA, eAx = slitYSamp (inst, slit, ySpacing)
	# Need to think about proper sampling here as well
	perPix = int ( ceil (len (lA) / float(inst.nPx) ) )
	pxA, pixSpacing = np.linspace (0, inst.nPx, perPix * inst.nPx, retstep = True)
	slitWidth = inst.slitWidth[slit] * inst.fieldAngle
	convolveS = 0.5 * (slitWidth / (band.anamorphicFactor * inst.pixScale))
	kSigma = convolveS / pixSpacing
	kern1 = np.array([exp(-n*n / (2*kSigma*kSigma)) 
		for n in range(int (floor(-6*kSigma)), int (ceil (6*kSigma)) + 1)])
	kern = kern1 / kern1.sum()
	def processRow (x, y):
		p0 = band.pxF (band.minL, x, y)
		p1 = band.pxF (band.maxL, x, y)
		pMin = max (int (floor (perPix * min (p0, p1))), 0)
		pMax = min (int (ceil (perPix * max (p0, p1)))+1, perPix * inst.nPx)
		pxC = pxA[pMin:pMax]
		lAr = np.squeeze (band.lF (pxC, x, y))
		iAr = np.zeros (pxA.size)
		iAr[pMin:pMax] = np.abs(slitWidth * inst.intensityConv * iF(lAr) * np.squeeze ( band.lF (pxC, x, y, dx=1) ) )
		# TODO - check that this is handling the edges properly
		iAconv = np.convolve (iAr, kern, 'same')
		# now downsample to the real pixel grid
		iAp = iAconv[perPix//2 :: perPix]
		return iAp
	offset = inst.slitX[slit]
	iAA = np.array ([processRow (offset + (y-yC) * inst.sinSlitTilt, y) for y in yA])
	xA = offset + (yA - yC) * inst.sinSlitTilt
	return xA*eAx, yA, iAA*eA[:,np.newaxis]

# does linear interpolation between 2d fns
class FInterp3D ():
	def __init__ (self, z0, zStep, zA, fA):
		self.z0 = z0
		self.zStep = zStep
		self.zA = zA
		self.fA = fA

	# does linear interpolation between 2d B-spline functions
	# TODO - think about how best to deal with z values
	# that don't fall within the range we've been given ...
	def __call__ (self, x, y, z, dx = 0, dy = 0, dz = 0):
		# TODO - check if x or y are outside allowed range
		#if abs(x) > 1.0000001 or abs(y) > 1.000001: print "FInter3D warning!"
		# a "proper" version should probably use a binary search
		# to locate the z stuff ... for now, we assume
		# that the z spacing is uniform and specified
		i0 = int (floor ((z-self.z0) / self.zStep))
		zN = len (self.zA)
		i0 = max (0, min (zN - 2, i0))
		i1 = i0 + 1
		z0 = self.zA[i0]
		z1 = self.zA[i1]
		v0 = interpolate.bisplev (x, y, self.fA[i0], dx, dy)
		v1 = interpolate.bisplev (x, y, self.fA[i1], dx, dy)
		# TODO - handle dz
		l = (z1 - z) / (z1 - z0)
		return l * v0 + (1 - l) * v1

def getBand (inst, bandName):
	band = Band ()
	band.name = bandName
	band.anamorphicFactor = inst.anamorphicFactors[bandName]
	condBand = inst.optData['band'] == band.name
	bandL = np.extract (condBand, inst.optData['lambda'])
	band.minL = np.min (bandL)
	band.maxL = np.max (bandL)
	#####
	yA = np.unique (inst.optData['yin'])
	def spl (y):
		cond = np.logical_and (condBand, inst.optData['yin'] == y)
		l = np.extract (cond, inst.optData['lambda'])
		x = np.extract (cond, inst.optData['xin'])
		px = np.extract (cond, inst.optData['px'])
		py = np.extract (cond, inst.optData['py'])
		# note - for my version of scipy, this seems to be
		# a "integer argument expected, got float" warning
		# This seems to be a feature of running 
		# scipy.interpolate.bisplrep with arrays of length >= 72
		# i.e.
		# scipy.interpolate.bisplrep (range(72), range(72), range(72)) 
		# gives the error, but
		# scipy.interpolate.bisplrep (range(71), range(71), range(71)) 
		# does not
		# TODO - check whether this happens on a more recent version of scipy
		lSpl = interpolate.bisplrep (px, x, l)
		ySpl = interpolate.bisplrep (px, x, py)
		xSpl = interpolate.bisplrep (l, x, px, xb = band.minL, xe = band.maxL)
		return lSpl, ySpl, xSpl
	lSpl, ySpl, xSpl = unzip([spl(y) for y in yA])
	band.lF = FInterp3D (yA[0], yA[1] - yA[0], yA, lSpl)
	band.pyF = FInterp3D (yA[0], yA[1] - yA[0], yA, ySpl)
	band.pxF = FInterp3D (yA[0], yA[1] - yA[0], yA, xSpl)
	return band

def processOpticalData (inst, optData):
	# convert from mm to pixels
	optData['px'] /= 0.001 * inst.pixSize
	optData['px'] += inst.nPx/2
	optData['py'] /= 0.001 * inst.pixSize
	optData['py'] += inst.nPy/2
	inst.optData = optData

# TODO - it's rather unfortunate that all of our data sources
# have rather different ideas about what the bands are
# It would be nice to get better data ...
# In the mean time, compromise as best we can ...
def applyTransferFunction (band, transfer, skyBg):
	lA, tA = unzip (transfer)
	# Note - we need to feed this a monotonically increasing
	# lA thing, otherwise it goes silly and gives NaNs (bit stupid ...)
	# In our case, we know that it's in reverse order, so we just reverse
	# things
	lAd, tAd = (np.flipud(lA), np.flipud(tA)) if (lA[0] > lA[1]) else (np.array(lA), np.array(tA))
	l0 = lAd[0]
	l1 = lAd[-1]
	tF = interpolate.InterpolatedUnivariateSpline (lAd, tAd)
	lA1, iA1 = zip(*[(l, i) for (l, i) in skyBg if (l >= band.minL and l <= band.maxL)])
	iA2 = tF (lA1) * iA1 * 1000.0 # 1000.0 for nm to um conversion
	cond = np.logical_and (lA1 >= l0, lA1 <= l1)
	iA3 = np.where (cond, iA2, 0.0)
	band.spectralIntensity = (lA1, iA3)

def genPoissonVals (im): return np.random.poisson (im)

def saveAsFits (arr, fname):
	hdu = pyfits.PrimaryHDU (arr)
	hdu.writeto (fname)

def raFromRad (rad):
	rH = rad * 12.0 / pi
	h = int (floor (rH))
	rMin = (rH - h)*60.0
	m = int (floor (rMin))
	s = (rMin - m)*60.0
	return (h % 24, m, s)

def dmsFromRad (rad):
	rDeg = rad * 180.0 / pi
	d = int (floor (rDeg))
	rMin = (rDeg - d)*60.0
	m = int (floor (rMin))
	s = (rMin - m)*60.0
	return (d % 180, m, s)

# Note - mechanically contiguous slits are given a single entry in
# the table, according to MSDN18
# Also - how do we indicate that there is no target in a given slit?
# Will that ever arise? Seems silly not to allow it.
def saveAsFitsWithExtensions (inst, im, fname):
	imHDU = pyfits.PrimaryHDU (im)
	####
	# Assume that long slits are laid out sensibly
	slitWidths = np.array ([inst.slitWidth[g[0]]*inst.fieldAngle 
		for g in inst.slitGroups])
	slitLengths = np.array ([inst.barPitch*len(g) - 2*inst.barGap 
		for g in inst.slitGroups])
	# need to calculate RA and DEC stuff ...
	# For the moment, just pretend that x is dec, y is ra
	ra0 = 0.0
	dec0 = 0.0
	arcsec = 2 * pi / (360 * 60 * 60)
	# want ra, dec of slit centre
	slitPos = [(dmsFromRad((dec0 + np.mean(np.array(inst.slitX)[g])*inst.fieldAngle)*arcsec), 
		raFromRad((ra0 + np.mean([slitYPos(inst, s)[2] for s in g])*inst.fieldAngle)*arcsec)) 
		for g in inst.slitGroups]
	slitDec, slitRa = unzip (slitPos)
	slitRaH, slitRaM, slitRaS = unzip (slitRa)
	slitDecD, slitDecM, slitDecS = unzip (slitDec)
	# should slit # be zero-based?
	nGroups = len(inst.slitGroups)
	c1 = pyfits.Column (name='slit number', format='J', array=np.arange(nGroups))
	c2 = pyfits.Column (name='slit RA hours', format='J', array=slitRaH)
	c3 = pyfits.Column (name='slit RA minutes', format='J', array=slitRaM)
	c4 = pyfits.Column (name='slit RA seconds', format='D', array=slitRaS)
	c5 = pyfits.Column (name='slit DEC degrees', format='J', array=slitDecD)
	c6 = pyfits.Column (name='slit DEC minutes', format='J', array=slitDecM)
	c7 = pyfits.Column (name='slit DEC seconds', format='D', array=slitDecS)
	c8 = pyfits.Column (name='slit width', format='D', array=slitWidths)
	c9 = pyfits.Column (name='slit length', format='D', array=slitLengths)
	# Don't deal with targets for the moment
	tId = ["No target"]*nGroups
	c10 = pyfits.Column (name='target id', format='10A', array=tId )
	c11 = pyfits.Column (name='target priority', format='J', array=[0]*nGroups)
	c12 = pyfits.Column (name='target location', format='D', array=[0.0]*nGroups)
	slitHDU = pyfits.new_table([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12])
	hdulist = pyfits.HDUList ([imHDU, slitHDU])
	hdulist.writeto (fname)
