# Absolute wavelength calibration
# Robert Lasenby 2009

from instrumentModel import *
import scipy.optimize

import matplotlib.pyplot as plt
import pdb

# TODO - refactor
def loadLineList (fname):
	f = open (fname, 'r')
	def readLn (l) : return (float(l[0])/10000.0, float(l[1]))
	d = [readLn(l.split()) for l in f if l[0] != '#']
	f.close()
	return d

class lineGroup():
	def __init__ (self, k, y, k0, k1):
		self.k = k
		self.y = y
		self.k0 = k0
		self.k1 = k1
		self.mightMerge = False

# TODO - document the grouping procedure
def groupLines (inst, band, h, lA1, iA1):
	# so, the necessary parameters ...
	readNoiseSigma = 10.0 # read noise of 10 counts --- too high, I spose ...
	readNoiseVar = readNoiseSigma * readNoiseSigma
	# theoretically, the number of samples
	# per pixel will basically be the number of rows in the slit ...
	# however, this assumes a totally reliable relative wavelength soln ...
	# let's try this for now ...
	nSampPerPixel = 10 
	# we need to get this from the simulator to start with ...
	# Number in pixels ...
	# Should, of course, base this on slit width (and similarly
	# for the simulator)
	# TODO - refactor all these constants
	# TODO - might also want to allow this to vary with wavelength? Check ...
	# or at least document how one would do this
	peakSpreadSigma = 1.4328993695242775 # K and H bands
	#peakSpreadSigma = 1.4565126924677487 # J and Y bands
	# we're also going to want some kind of characterisation
	# of how the theoretical amplitudes can vary ...
	# anyway, let's go for this :)
	integrationTime = 30.0 # 
	lineCountMax = 800.0 * integrationTime # guess based on the H band long-slit image ...
	#
	condBand = np.logical_and (lA1 > band.minL, lA1 < band.maxL)
	lA1b = lA1[condBand]
	iA1b = iA1[condBand]
	pxA1 = h(lA1b) * (0.5 * inst.nPx)
	#
	iSort = np.argsort (pxA1)
	lA = lA1b[iSort]
	pxA = pxA1[iSort]
	iA = iA1b[iSort]
	iMax = np.max (iA)
	# TODO - compensate for maximum of transfer function as well ...
	# think about exactly what's needed here ...
	transferVals = band.transfer (lA)
	iA *= transferVals
	iA *= lineCountMax / iMax
	# Now, time to do the groups stuff ...
	# let's set the multiplier sensibly ...
	ePeak = 2.0 # complete guess ftm ...
	def inI (a, (b,c)): return (a > b and a < c)
	def combineQ (k0, y0, k1, y1):
		d = (k1-k0)/peakSpreadSigma
		if d < 2:
			return 2 # always
		else:
			dis = sqrt(d*d - 4)
			t = 0.5 * (d*(d + dis)-2) * exp (-0.5 * d * dis)
			c0 = inI (ePeak*y0, (t*y1, y1)) or inI (y1, (t*ePeak*y0, ePeak*y0))
			c1 = inI (y0, (t*ePeak*y1, ePeak*y1)) or inI (ePeak*y1, (t*y0, y0))
			c2 = inI (y0, (y1, ePeak*y1)) or inI (y1, (y0, ePeak*y0))
			if c0 and c1:
				return 0 # never
			elif c2 or c0 or c1:
				return 1 # sometimes
			else:
				return 2 # always
	def approxPeak (x0, y0, x1, y1):
		# weighted average as first guess... then, we want to solve ...
		x = (x0*y0 + x1*y1)/(y0 + y1)
		# well, let's just go with the weighted average ftm ... see how much
		# difference it makes to do it properly later ...
		return x
	def addLine (g, k, y):
		kd = approxPeak (g.k, g.y, k, y)
		k0d = approxPeak (g.k0, ePeak * g.y, k, y/ePeak)
		k1d = approxPeak (g.k1, g.y/ePeak, k, ePeak*y)
		yd = g.y * exp (- (kd-g.k)*(kd-g.k) / (2.0 * peakSpreadSigma * peakSpreadSigma))
		yd += y * exp (- (kd-k)*(kd-k) / (2.0 * peakSpreadSigma * peakSpreadSigma))
		g.k = kd
		g.k0 = k0d
		g.k1 = k1d
		g.y = yd
	groups = []
	nPx = len(pxA)
	g0 = None
	for i in xrange(nPx):
		l = lA[i]
		k = pxA[i]
		y = iA[i]
		if g0 == None:
			g0 = lineGroup (k, y, k, k)
		else:
			c0 = combineQ (g0.k1, g0.y, k, y)
			c1 = combineQ (g0.k0, g0.y, k, y)
			if c0 == 2 and c1 == 2: # always, always
				addLine (g0, k, y)
			elif c0 == 0 and c1 == 0: # never, never
				groups += [g0]
				g0 = lineGroup (k, y, k, k)
			else:
				g0.mightMerge = True
				groups += [g0]
				g0 = lineGroup (k, y, k, k)
				g0.mightMerge = True
	# TODO - refactor
	cy = 0.125 # for Fowler-8
	# Now .. we want to get the measurement error for each peak ...
	for g in groups:
		g.measurementW = sqrt(pi) * g.y * g.y * nSampPerPixel / (2 * peakSpreadSigma * (readNoiseVar + cy * g.y))
		g.pW = g.measurementW + (g.k1 - g.k0)/2
	return groups
	#return pxA, iA, groups

def fitPeak (inst, (nA, iA, wA, nu2i), nu0):
	# So, our first trick is to "segment out"
	# the correct part of the arrays
	windowR = 5.0 * (2.0/inst.nPx) # TODO - refactor
	# TODO - might want to reconsider the nu2i approach ...
	nuMin = nu0 - windowR
	nuMax = nu0 + windowR
	iMin = max (0, int (floor( nu2i (nuMin))))
	# TODO - make sure that we don't go out of range here
	nNu = len(nA)
	while nA[iMin] < nuMin: iMin += 1
	while iMin > 0 and nA[iMin] > nuMin: iMin -= 1
	iMax = min (len(nA)-1, int (ceil(  nu2i (nuMin))))
	while nA[iMax] > nuMax: iMax -= 1
	while iMax < nNu-1 and nA[iMax] < nuMax: iMax += 1
	# Windowed versions of the arrays
	nAw = nA[iMin:iMax+1]
	iAw = iA[iMin:iMax+1]
	wAw = wA[iMin:iMax+1]
	# TODO - refactor
	peakS = 1.4328993695242775 * (2.0/inst.nPx) # K and H bands
	#peakSpreadSigma = 1.4565126924677487 * (2.0/inst.nPx) # J and Y bands
	tw = np.sum (wAw)
	iwA = iAw * wAw
	def f(nu):
		pA = np.exp (- (nAw - nu)**2 / (2.0 * peakS * peakS))
		tpw = np.sum (pA * wAw)
		pA -= tpw / tw
		dp = np.dot (iwA, pA)
		# think about whether to return squared version ...
		return - np.dot (iwA, pA) / np.sqrt (np.sum (wAw * np.square (pA)))
	# Now ... the trick is, we need to start with a bracketing
	# interval. Well ... shall we just assume we've
	# got one, throw an error otherwise?
	# Bear in mind that we might want to play around with
	# number of iterations allowed, etc ...
	searchR = 2.0 * (2.0/inst.nPx)
	#
	nuA, nuB, nuC = nu0 - searchR, nu0, nu0 + searchR
	# Check if we have a bracketing interval
	# If we don't, we should probably do some kind of
	# search ... at the moment, just give up :)
	fA = f(nuA); fB = f(nuB); fC = f(nuC)
	if fB > fA or fB > fC:
		return nuB, 0.0
	nuArgMin, fMin, _, _ = scipy.optimize.brent (f, brack=(nuA, nuB, nuC), full_output=True)
	# Okay ... we want to return a "real" weight here ...
	# Basically, numerically evaluate second derivative at peak ...
	dNu = 1e-6
	# remember that actual log-posterior is f^2/2 ...
	fl = f(nuArgMin - dNu); fh = f(nuArgMin + dNu)
	wNu = 1e12 * (fMin*fMin - 0.5*(fl*fl + fh*fh))
	#
	#return nuArgMin, 3.0 * (inst.nPx * inst.nPx * 0.25) # TODO - proper weight
	return nuArgMin, wNu # TODO - proper weight

# where spec = (nA, iA, wA, nu2i)
# g : \kappa -> \nu
def fitPeaks (inst, groups, g, spec):
	fitList = []
	# TODO - sensible number here
	# Also, think about using height instead ...
	wThreshold = 4.0 * (inst.nPx * inst.nPx * 0.25) # \sigma < 0.5
	pWThreshold = 4.0 # i.e. \sigma < 0.5
	nuMin = np.min (spec[0])
	nuMax = np.max (spec[0])
	for group in groups:
		if group.pW < pWThreshold:
			break # since the groups are ordered by this field
		if group.mightMerge:
			continue
		# Scale change ... look at what the errors mean as well ...
		gk = 2.0 * group.k  / inst.nPx
		if gk <= nuMin or gk >= nuMax:
			continue
		nu1 = g (gk)
		nu, wnu = fitPeak (inst, spec, nu1)
		if wnu > wThreshold:
			fitList += [(gk, nu, wnu)]
			g.updatePoint (gk, nu, wnu)
	return fitList

# where spec = (nA, iA, wA, nu2i)
# Should probaby pass through the "exposure parameters" here ...
def absoluteWavelengthCalibrateSpectrum (inst, band, spec, h, lineList):
	g = glm1D (getLegendre1D (5))
	# Should do coarse alignment - wait on implementation of this
	# until have the relative alignment sorted
	# Okay ... how to init g?
	# FIXME - look at this hack
	w = (1.0/5.0) * inst.nPx
	nuInit = np.linspace (-1.0, 1.0, 10)
	g.updateFull (nuInit, nuInit, np.ones (10) * w)
	lA1, iA1 = unzip (lineList)
	lA1 = np.array (lA1); iA1 = np.array (iA1)
	# Now, want to get our groupings ...
	# So, we need to have h ...
	groups = groupLines (inst, band, h, lA1, iA1)
	#pdb.set_trace()
	# Sort groups by predicted error (descending weight)
	groups.sort (lambda g1, g2 : cmp (g2.pW, g1.pW))
	# Then do the peak-fitting update
	fitList = fitPeaks (inst, groups, g, spec)
	# TODO - fix this update stuff ...
	fitK, fitNu, fitW = unzip (fitList)
	fitK = np.array (fitK); fitNu = np.array (fitNu); fitW = np.array (fitW)
	g.updateFull (fitK, fitNu, fitW)
	#pdb.set_trace ()
	return g

# Alright ... we also need to flat-field the science image ...
# TODO - refactor
def flatFieldCompensate (sciIm, sciW, ffIm, ffW):
	sciImFF = sciIm / ffIm
	# TODO - check this formula (comes from Wikipedia, which claims it's
	# from a series expansion)
	sciImFFw = np.square (ffIm) * ffW * sciW / (ffW + np.square (sciImFF) * sciW)
	infF = np.isinf(sciImFF)
	nanF = np.isnan(sciImFF)
	sciImFFw[np.isnan(sciImFFw)] = 0.0
	sciImFFw[np.isinf(sciImFFw)] = 0.0
	sciImFFw[infF] = 0.0
	sciImFFw[nanF] = 0.0
	sciImFF[infF] = 0.0
	sciImFF[nanF] = 0.0
	return sciImFF, sciImFFw

# Okay ... so, we need to obtain nu2i, h in some way ...
def resample1D (oldx, y, newx):
	iF = scipy.interpolate.InterpolatedUnivariateSpline (oldx, y)
	return iF (newx)

def get_nu2i (nA, nSteps = 20):
	n = len(nA)
	step = n/nSteps
	iAsamp = np.arange (0, n, step)
	return scipy.interpolate.InterpolatedUnivariateSpline (nA[iAsamp], iAsamp)

# So, first task to is prepare our spectrum for processing ...
# how best to do this? Let's spose, at first, that we've
# got ... well, how about the lF : p_x, y -> \lambda function?
# And ... maybe the pyF : p_x, y -> p_y function?
def prepareSpectrum (inst, band, slit, im, imW):
	# TODO - refactor this pattern?
	y0, y1 = inst.slitEdges (slit)
	slitX = slit['slitX']
	slitY = slit['slitY']
	slitX0 = slitX - inst.tanSlitTilt * slitY
	# TODO - probably don't want to use the "band" values ... want
	# to pass things through to replace these ...
	pyF = band.pyF3.xyF (slitX0)
	lF = band.lF3.xyF (slitX0)
	#
	pxCn = inst.px2pxn (np.arange (inst.nPx))
	yC = np.linspace (y0, y1, 40) # FIXME - magic number
	pyAA = np.transpose([inst.pyn2py (pyF (np.ones_like(yC)*px, yC)) for px in pxCn])
	# TODO - think about whether we have any off-by-one problems
	# here and in similar situations
	py0 = max (0, int (floor (np.min (pyAA))))
	py1 = min (inst.nPy-1, int (ceil (np.max (pyAA))))
	pyC = np.arange (py0, py1+1)
	nPyW = len(pyC)
	#
	yAA = np.transpose ([ resample1D (np.flipud(pyAA[:,i]), np.flipud (yC), pyC) for i in xrange(inst.nPx)])
	lAA = lF(np.tile (pxCn, nPyW), np.ravel(yAA)).reshape ((nPyW, inst.nPx))
	#
	pxCns = np.linspace (-1.0, 1.0, 20)
	yCentre = 0.5 * (y0 + y1)
	lAs = lF (pxCns, np.ones_like (pxCns) * yCentre)
	# TODO - make sure that we always need to flip here ...
	# and probably refactor this "decision" in some sensible way?
	f0Inv = scipy.interpolate.InterpolatedUnivariateSpline (np.flipud (lAs), np.flipud (pxCns))
	#
	nAA = f0Inv(np.ravel(lAA)).reshape ((len(pyC), inst.nPx))
	#
	# Now, we look at the actual image ...
	nA0 = np.ravel(nAA)
	iA0 = np.ravel(im[py0:py1+1])
	wA0 = np.ravel(imW[py0:py1+1])
	# cut stuff that's not "suitable"
	cond1 = wA0 > 0.2 # FIXME - magic number
	cond1 = np.logical_and (cond1, np.ravel(yAA) > y0 + 2 * (2.0 / inst.nPx))
	cond1 = np.logical_and (cond1, np.ravel(yAA) < y1 - 2 * (2.0 / inst.nPx))
	#
	nA2 = nA0[cond1]
	iSort = np.argsort (nA2)
	nA1 = nA2[iSort]
	iA1 = (iA0[cond1])[iSort]
	wA1 = (wA0[cond1])[iSort]
	# Okay ... that gives us our spectrum to work with ...
	nu2i = get_nu2i (nA1)
	return (nA1, iA1, wA1, nu2i), f0Inv

# Now, the overall trick ...
def absoluteWavelengthCalibrateSlit (inst, band, slit, im, imW, lineList):
	spec, h = prepareSpectrum (inst, band, slit, im, imW)
	return absoluteWavelengthCalibrateSpectrum (inst, band, spec, h, lineList)