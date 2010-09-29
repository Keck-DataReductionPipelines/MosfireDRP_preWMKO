# Tracing slit edges in flat field images
# Robert Lasenby 2009

from instrumentModel import *
import scipy.ndimage
import scipy.optimize

# Do we want to include qe uncertainty? Would give more realistic values ...
# Don't bother for the moment, I spose ...
# Gives approximate variance for "true" flux (not taking into account
# binomial qe noise ftm)
def qeCompensate (ff, ffW, qe, poissComp=True):
	# Seem to remember that pixels below a certain qe aren't
	# likely to be working properly - ask about this ...
	# For the moment, assume that pixels of qe <= 2.0 are "broken"
	qe0 = 0.2
	qeCond = qe > qe0
	ffWd = np.where (qeCond, ffW * np.square (qe), 0.0)
	# Is it important that we don't have NaNs here? Probably an idea ...
	ffd = np.where (qeCond, ff/qe, ff)
	if poissComp:
		# clipping the -ve stuff ... all a bit of a hack atm
		ffWd /= (1.0 + np.maximum (ffd, 0.0) * ffWd)
	return ffd, ffWd

# Hmmm ... if we're going to be doing this for multiple things, we want to
# save out the stuff that doesn't depend on the kernel, I spose ...
# Put in that optimisation later
def edgeDetectWeighted (ff, ffW, ker, origin=0):
	fyfw = ff * ffW # y_i / s_i^2
	boxcar = np.ones_like (ker)
	fyfwBox = scipy.ndimage.filters.convolve1d (fyfw, boxcar, axis=0, mode='constant', origin=origin)
	wBox = scipy.ndimage.filters.convolve1d (ffW, boxcar, axis=0, mode='constant', origin=origin)
	wP = scipy.ndimage.filters.convolve1d (ffW, ker, axis=0, mode='constant', origin=origin)
	fyfwP = scipy.ndimage.filters.convolve1d (fyfw, ker, axis=0, mode='constant', origin=origin)
	pW = wP / wBox
	ypW = pW * fyfwBox
	#
	zI = scipy.ndimage.filters.convolve1d (ffW, np.square(ker), axis=0, mode='constant', origin=origin)
	zI -= pW * wP
	#
	ffConvW = (fyfwP - ypW) / np.sqrt(zI)
	ffConvW[np.isinf(ffConvW)] = 0.0
	ffConvW[np.isnan(ffConvW)] = 0.0
	return ffConvW

# Assumes that pImB, pImT have had l -> l^2 / 2 operation carried out
def iterativeOpt (inst, pyF, bhStart, y0, y1, pImB, pImT):
	bh0 = np.array(pyF.bh)
	# TODO - refactor
	# far-too-simple prior ... think about sensible value ...
	mb = len(bh0)
	wB = np.ones_like (bh0) * mb * mb # \sigma = 1/m_\beta
	wB[::mb] = 1.0 / (2.0 * 2.0) # wide prior on those terms x^i y^0
	#wB[1::mb] = 1.0 / (2.0 * 2.0) # wide prior on those terms x^i y^1
	wB[0] = 1.0 / (5.0*5.0) # even wider prior on the constant term
	# Okay ... now we need to convert these into a suitably
	# normalised form ...
	wB *= 0.25 * inst.nPx * inst.nPx
	# design matrix doesn't change throughout ...
	pxA = inst.px2pxn (np.arange (inst.nPx))
	yA0 = np.tile (y0, inst.nPx)
	yA1 = np.tile (y1, inst.nPx)
	gm0 = pyF.gm (pxA, yA0)
	gm1 = pyF.gm (pxA, yA1)
	gm0T = np.transpose (gm0)
	gm1T = np.transpose (gm1)
	# Get the necessary interpolations
	py0_0 = inst.pyn2py (np.dot (gm0, bh0))
	py1_0 = inst.pyn2py (np.dot (gm1, bh0))
	windowPx = 10
	pyB0 = max (0, floor (np.min (py0_0) - windowPx))
	pyB1 = min (inst.nPy-1, ceil (np.max (py0_0) + windowPx))
	pyT0 = max (0, floor (np.min (py1_0) - windowPx))
	pyT1 = min (inst.nPy-1, ceil (np.max (py1_0) + windowPx))
	blf = [scipy.interpolate.InterpolatedUnivariateSpline (np.arange(pyB0, pyB1+1),
		pImB[pyB0:pyB1+1, px])
		for px in xrange(inst.nPx)]
	tlf = [scipy.interpolate.InterpolatedUnivariateSpline (np.arange(pyT0, pyT1+1),
		pImT[pyT0:pyT1+1, px])
		for px in xrange(inst.nPx)]
	#
	def getL (b):
		# now, how to evaluate the prior probability here? Well ... 
		priorL = - 0.5 * np.sum (wB * np.square (b - bh0))
		# Now, the posterior ...
		# now ... evaluate the py positions ...
		# Need to convert these in "real" pixel things ...
		# Also, need to make sure that we treat things
		# that've gone outside the interpolation area sensibly ...
		py0 = inst.pyn2py (np.dot (gm0, b))
		py1 = inst.pyn2py (np.dot (gm1, b))
		# This loop will probably be very slow ...
		# how to speed it up?
		#postL = np.sum(blf[px](py0[px]) + tlf[px](py1[px]) for px in xrange(inst.nPx))
		#return priorL + postL
		for px in xrange(inst.nPx):
			py0px = py0[px]
			if py0px >= pyB0 and py0px <= pyB1:
				priorL += blf[px](py0px)[0]
			py1px = py1[px]
			if py1px >= pyT0 and py1px <= pyT1:
				priorL += tlf[px](py1px)[0]
			#priorL += blf[px](py0[px])[0] + tlf[px](py1[px])[0]
		return -priorL
	# gradient of L
	def getLp (b):
		# Prior contribution
		priorP = - wB * (b - bh0)
		# Now for the posterior ...
		py0 = inst.pyn2py (np.dot (gm0, b))
		py1 = inst.pyn2py (np.dot (gm1, b))
		condB = np.logical_or (py0 < pyB0, py0 > pyB1)
		condT = np.logical_or (py1 < pyT0, py1 > pyT1)
		bp = np.array([blf[px](py0[px],nu=1)[0] * 0.5 * inst.nPy for px in xrange (inst.nPx)])
		tp = np.array([tlf[px](py1[px],nu=1)[0] * 0.5 * inst.nPy for px in xrange (inst.nPx)])
		bp[condB] = 0.0
		tp[condT] = 0.0
		# Okay ... now ... what's the best way to do this?
		return -(priorP + np.dot (gm0T, bp) + np.dot (gm1T, tp))
	# Hessian of L
	def getLh (b):
		priorH = np.diag (wB)
		py0 = inst.pyn2py (np.dot (gm0, b))
		py1 = inst.pyn2py (np.dot (gm1, b))
		condB = np.logical_or (py0 < pyB0, py0 > pyB1)
		condT = np.logical_or (py1 < pyT0, py1 > pyT1)
		bpp = np.array([blf[px](py0[px],nu=2)[0]*0.25*inst.nPy*inst.nPy for px in xrange (inst.nPx)])
		tpp = np.array([tlf[px](py1[px],nu=2)[0]*0.25*inst.nPy*inst.nPy for px in xrange (inst.nPx)])
		bpp[condB] = 0.0
		tpp[condT] = 0.0
		priorH -= np.dot (gm0T * bpp, gm0)
		priorH -= np.dot (gm1T * tpp, gm1)
		return priorH
	bhStart = np.array (bhStart)
	bOpt = scipy.optimize.fmin_ncg (getL, bhStart, fprime = getLp, fhess = getLh, avextol = 1e-6)
	return bOpt

# The "SVD" method as well ...

# returns both argmax and width
def quadraticMaxEven ((fa, fb, fc)):
	denom = 2*fb - fa - fc
	return 0.5 * (fc - fa) / denom, denom

# find maximum value in range, perform quadratic interpolation
# using surrounding pixels to guess the subpixel maximum
# Also ... where do we want to do the "rejection" of insufficiently-strong
# peaks?
def getColMaxima (py0, py1, pIm):
	nPy = len(pIm)
	nPx = len(pIm[0])
	pxA = np.arange (nPx)
	pyMin = max(0, int(floor(np.min(py0))))
	# TODO - think about the off-by-one stuff here
	pyMax = min(nPy, int(ceil(np.max(py1))))
	pWindow = pIm[pyMin:pyMax,:]
	argmaxPy = pyMin + np.argmax (pWindow, axis=0)
	# Want to do the filtering stage here?
	maxVals = pIm [argmaxPy, pxA]
	# So ... reject things which are on the edge of the range
	# ... or should we just make sure they're not off the edge 
	# of the image?
	c1 = np.logical_and (argmaxPy > pyMin, argmaxPy < pyMax-1)
	# TODO - factor out the magic level, or at least document it
	c1 = np.logical_and (c1, maxVals > 2.0)
	px1 = pxA[c1]
	# and quadratic-interpolate at the maximum
	# TODO - this idiom again ... sort it out ...
	argmaxSubpix, weightSubpix = unzip ([
		(lambda py : py + quadraticMaxEven(pIm[py-1:py+2,px]))(argmaxPy[px])
		#(lambda py : quadraticMaxEven(0.5*np.square(pIm[py-1:py+2,px])))(argmaxPy[px])
		for px in px1 ])
	# Now ... return this stuff?
	#return px1, argmaxPy[px1] + np.array(argmaxSubpix), np.array (weightSubpix)
	return px1, np.array(argmaxSubpix), np.array (weightSubpix)

# Update from the column maxima found ...
def updateFromMaxima (inst, pyF, y0, y1, pImB, pImT):
	pxAn = inst.px2pxn (np.arange(inst.nPx))
	pyBottomN = pyF (pxAn, np.ones_like(pxAn)*y0)
	pyTopN = pyF (pxAn, np.ones_like(pxAn)*y1)
	#
	pyBottom = inst.pyn2py (pyBottomN)
	pyTop = inst.pyn2py (pyTopN)
	pbX, pbY, pbW = getColMaxima (pyBottom - 5.0, pyBottom + 5.0, pImB)
	ptX, ptY, ptW = getColMaxima (pyTop - 5.0, pyTop + 5.0, pImT)
	# So, let's get our design matrix for the update ...
	pxA = inst.px2pxn (np.concatenate ((pbX, ptX)))
	yA = np.concatenate ((np.ones_like(pbX)*y0, np.ones_like(ptX)*y1))
	pyA = inst.py2pyn (np.concatenate ((pbY, ptY)))
	wA = (0.5 * inst.nPy) * np.sqrt(np.concatenate ((pbW, ptW)))
	# Okay ... how're we going to choose the noise level for our coefficients?
	# Hmmm ... thoughts? Well ... basically, we decide on the "noise"
	# level we want, I spose ... hmmm ... okay ... thoughts? Problem
	# is that we're not properly centred .... 
	# Well, let's decide ... say we want to have ~2 pix stddev ...
	# in the best or the worse case? Well, let's say the best case, i.e.
	# x^T x \approx 2, then with noise increasing from there on in ...
	wb0 = 2.0 / 4.0
	wb0 *= 0.25 * inst.nPy * inst.nPy # convert into correct values ...
	gm1 = pyF.gm (pxA, yA)
	pyF.reg1Bh (gm1, pyA, wA, wb0)
	return

########

# More "practical" side of things ...

# Okay ... thoughts ... do we want to do the SVD solution first, or not?
# Hmmm ... well, can't really hurt, can it? :)

# Provide a "prototype" shape
# In reality, we'd want to derive this from real instrument data
# TODO - probably want to refactor this
def shape1 (y):
	return np.where (y <= 0, np.exp (- (y**2 / 2.0)), 1.0)

def calibFlat (inst, band, ffHduList, qe, shape, slits):
	ff, ffW = qeCompensate (ffHduList[0].data, ffHduList[1].data, qe)
	#
	ker = shape (np.arange (-3, 8) - inst.edgeOffset)
	ker -= np.mean (ker)
	ker /= np.linalg.norm (ker)
	#
	# TODO - should we just flip the kernel?
	pImB = edgeDetectWeighted (ff, ffW, ker, -2)
	pImT = edgeDetectWeighted (ff, ffW, np.flipud (ker), 2)
	pImB *= 0.5 * pImB * np.sign (pImB)
	pImT *= 0.5 * pImT * np.sign (pImT)
	#
	# Now call all the appropriate functions ...
	def getPyF (slit):
		# TODO - refactor this pattern?
		slitX = slit['slitX']
		slitY = slit['slitY']
		slitX0 = slitX - inst.tanSlitTilt * slitY		
		pyF = band.pyF3.xyF (slitX0)
		# HACK - return the first guess
		return pyF
		y0, y1 = inst.slitEdges (slit)
		#
		bh0 = pyF.bh
		updateFromMaxima (inst, pyF, y0, y1, pImB, pImT)
		bhStart = pyF.bh
		pyF.bh = bh0
		bh1 = iterativeOpt (inst, pyF, bhStart, y0, y1, pImB, pImT)
		pyF.bh = bh1
		return pyF
	#
	return [getPyF (slit) for slit in slits]
