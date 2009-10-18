# Coarse Absolute Wavelength Calibration
# Robert Lasenby 2009

from instrumentModel import *
import scipy.special
import scipy.optimize
import scipy.interpolate
from scipy import weave
from scipy.weave import converters

import pdb

# Alright - the idea is that we start with iAA, nAA (for the image)
# and iA, kA (for the "model sky"), and we want to find
# g : k -> n such that things match up (we also need to do the convolution
# stuff, of course ...)

# slitIm = (iAA, wAA, nAA, cond)
# h : \kappa -> \lambda
# hInv : \lambda -> \kappa
# Do we guarantee that these things agree? Do we need to? Think
# a bit about this issue ..
def coarseAbsoluteWavelengthCalibrateSlit (inst, band, slitIm, iASky, lASky, h, hInv):
	# where to get the convolution width from?
	# TODO - refactor
	peakS = 1.4328993695242775 * (2.0/inst.nPx) # K and H bands (in normalised units)
	# deltaTarget in units of kappa
	deltaTarget = 0.5 * (2.0 / inst.nPx)
	skyInfo = prepareSkySpectrum (iASky, lASky, h, hInv, peakS, deltaTarget)
	# need to choose fAb and fR, then ...
	fR = peakS / 2.0 # is this justified? Think about ...
	# Okay ... now, what's this fAb thing going to be?
	nAb, imageInfo = prepareImage (inst, slitIm, fR)
	nA0, kA0, kStep = optimiseG (nAb, imageInfo, skyInfo)
	#optimiseG (nAb, imageInfo, skyInfo, slitIm)
	return nA0, kA0, kStep

# Okay ... given that we've precompute the appropriate things, we want to
# write out the optimisation process ...

def optimiseG (nAb, (xtiw, xtwV, xtwx0, xtwx1), (k0, k1, iF)):
#def optimiseG (nAb, (xtiw, xtwV, xtwx0, xtwx1), (k0, k1, iF), (iAA, wAA, nAA, cond)):
	# Firstly, we do a full-scale translational thing over a fairly wide
	# search range ...
	# Then, how exactly is our lerp sequence going to go?
	# Well, first two stages will deal with the whole b vector,
	# then we'll progress to modifying different parts of it ...
	# Structure ... ?
	# Special-case the first translational thing ...
	# want some rather extended nB here, I spose ...
	# check which direction we should shift in ...
	# Okay - how should b and bShift be related?
	def f0 (b):
		# see how well numpy does here ... may be a job for C?
		return np.dot (xtiw, b) / sqrt (np.dot (b * xtwx0, b) + 
				2.0 * np.dot (b[1:] * xtwx1, b[:-1]) -
				np.dot (xtwV, b))
	def f (kA):
		b = iF (kA)
		b[np.logical_or (kA < k0, kA > k1)] = 0.0
		# see how well numpy does here ... may be a job for C?
		return np.dot (xtiw, b) / sqrt (np.dot (b * xtwx0, b) + 
				2.0 * np.dot (b[1:] * xtwx1, b[:-1]) -
				np.dot (xtwV, b))
	# Go up to here, make sure that everything's working as it should
	# So, search the space of translations at first ... then ...
	# at first, shall we implement without the incremental stuff,
	# just for simplicity? Note that it would be possible to make
	# cut down the work in a localising manner, but don't actually do
	# it just yet ...
	# So, how to represent and work with our n -> k function?
	# Basically, want to add control pts at each stage, and mess
	# with different of them, in some sense ... hmmm ... could just operate
	# directly, I spose ... might be easiest ...
	#
	# TODO - more sensible choice of step size ...
	kStep = 0.25 / 1024.0 # atm, 1/4 of a pixel.
	kA0 = np.copy (nAb)
	#
	# Firstly, then, search over rigid translations
	nSteps = 64 # i.e. a 16-pixel search ... seems reasonable, can always change later
	bestL = -1e10
	bestShift0 = 0.0
	for kS in xrange (-nSteps, nSteps+1):
		kShift = kS * kStep
		ll = f (kA0 + kShift)
		if ll > bestL:
			bestL = ll
			bestShift0 = kShift
	kA0 += bestShift0
	#pdb.set_trace()
	#
	nlevels = 4
	# Okay ... hmmm ... bestL is a little annoying then ...
	# how're we sposed to modify it properly 
	def optInteriorPoint (i0, i1, i2):
		k0 = kA0[i0]; k1 = kA0[i1]; k2 = kA0[i2]
		# So ... search over what range? What should be the gap?
		# Hmmm ... range thing is a bigger issue ... the problem, 
		# I spose, is that we could a rather skewed representation
		# of our function if we allow things to squish up too much
		# - i.e. if we leave some areas undersampled. Hmmm ...
		# well ... hopefully that won't be too much of a problem ... anyway ...
		# halve the search range each time or something?
		# Seems like a vaguely sensible hack for the moment ...
		# To check - these should be power-of-two things?
		n1 = i1-i0+1
		n2 = i2-i1+1
		bestShift = 0.0
		newBest = bestL
		for kS in xrange (-nSteps, nSteps+1):
			kShift = kS * kStep
			k1d = k1 + kShift
			kA0[i0:i1+1] = np.linspace (k0, k1d, n1)
			kA0[i1:i2+1] = np.linspace (k1d, k2, n2)
			ll = f (kA0)
			if ll > newBest:
				newBest = ll
				bestShift = kShift
		k1d = k1 + bestShift
		kA0[i0:i1+1] = np.linspace (k0, k1d, n1)
		kA0[i1:i2+1] = np.linspace (k1d, k2, n2)
		return newBest
	# Rather inelegant way of doing this - how to improve?
	def optLeftEnd (i1):
		i0 = 0
		k0 = kA0[i0]; k1 = kA0[i1]
		n1 = i1 - i0 + 1
		bestShift = 0.0
		newBest = bestL
		for kS in xrange (-nSteps, nSteps+1):
			kShift = kS * kStep
			k0d = k0 + kShift
			kA0[i0:i1+1] = np.linspace (k0d, k1, n1)
			ll = f (kA0)
			if ll > newBest:
				newBest = ll
				bestShift = kShift
		k0d = k0 + bestShift
		kA0[i0:i1+1] = np.linspace (k0d, k1, n1)
		return newBest
	nbf = len (nAb)
	def optRightEnd (i1):
		i2 = nbf-1
		k1 = kA0[i1]; k2 = kA0[i2]
		n2 = i2-i1+1
		bestShift = 0.0
		newBest = bestL
		for kS in xrange (-nSteps, nSteps+1):
			kShift = kS * kStep
			k2d = k2 + kShift
			kA0[i1:i2+1] = np.linspace (k1, k2d, n2)
			ll = f (kA0)
			if ll > newBest:
				newBest = ll
				bestShift = kShift
		k2d = k2 + bestShift
		kA0[i1:i2+1] = np.linspace (k1, k2d, n2)
		return newBest
	for k in xrange (nlevels):
		nControl = 2**k + 1
		controlI = np.floor (np.linspace (0, nbf-1, nControl)).astype (int)
		# Firstly, optimise the odd-index control pts
		# Special-case k=0, I spose ....
		if k > 0:
			for i in xrange (2**(k-1)):
				i0 = controlI[2*i]
				i1 = controlI[2*i+1]
				i2 = controlI[2*i+2]
				bestL = optInteriorPoint (i0, i1, i2)
		# then the even-index interior ones
		if k > 1:
			for i in xrange (2**(k-1)-1):
				i0 = controlI[2*i+1]
				i1 = controlI[2*i+2]
				i2 = controlI[2*i+3]
				bestL = optInteriorPoint (i0, i1, i2)
		# now, the end points
		bestL = optLeftEnd (controlI[1])
		bestL = optRightEnd (controlI[-2])
		nSteps /= 2
	nControl0 = 2**(nlevels-1) + 1
	controlI0 = np.floor (np.linspace (0, nbf-1, nControl0)).astype (int)
	return nAb[controlI0], kA0[controlI0], kStep

# So, our first task is to prepare the sky spectrum - that is, we want
# to produce an interpolated version of the convolved thing. So ...

# nA the chosen locations of the basis functions ...
# actually, should we just choose the correct part of the band
# and be done with it? Probably ...
# h : \kappa -> \lambda (?)
def prepareSkySpectrum (iASky, lASky, h, hInv, convWidth, deltaTarget):
	# Can this go weird with strong peaks etc? Check ... might want to be a bit
	# careful here ... if it looks dodgy, just lerp, I spose ...
	iF = scipy.interpolate.InterpolatedUnivariateSpline (lASky, iASky)
	# Figure out what our sampling rate in kappa should be
	# Note the h' < 0 (i.e. \kappa and \lambda anticorrelate)
	k0 = hInv (lASky[-1])[0]
	k1 = hInv (lASky[0])[0]
	epsilonKappa = 1e-5 # TODO - epsilon <= 1e-8 seems to cause problems ... why?
	hp0 = -(h(k0+epsilonKappa)[0] - h(k0)[0])/epsilonKappa
	hp1 = -(h(k1+epsilonKappa)[0] - h(k1)[0])/epsilonKappa
	hpMax = max (hp0, hp1)
	# assumes a constant sampling rate for lASky ...
	deltaKappa = (lASky[1] - lASky[0])/hpMax
	kA = np.arange (k0, k1, deltaKappa)
	lAk = h (kA)
	# Don't bother here about any deriv-based adjustments - we're mostly
	# interested in local rather than global structure
	iAk = iF (lAk)
	# Now, to convolve ... what do we assume convWidth is given in here?
	kSigma = convWidth / deltaKappa
	kern1 = np.array ([exp(-n*n / (2*kSigma*kSigma)) 
		for n in range(int (floor(-6*kSigma)), int (ceil (6*kSigma)) + 1)])
	kern = kern1/kern1.sum()
	iAConv = np.convolve (iAk, kern, 'same')
	# Okay - we've now got a sensible convolved version ... and it's going
	# to be much too finely sampled. So, the sensible thing to do is to
	# obtain an efficient representation ... meaning that we should downsample,
	# then interpolate, I spose.
	# How about ... sample at the 1/2 pixel level, or something similar?
	stride = int (floor (deltaTarget / deltaKappa))
	assert stride != 0
	# TODO - need to make sure that we're within the correct range ...
	iConvF = scipy.interpolate.InterpolatedUnivariateSpline (kA[::stride], iAConv[::stride])
	return (k0, k1, iConvF)

# Now, want to get i^T X thing ... steal the C code here, I spose ...
def prepareImage (inst, (iAA, wAA, nAA, cond), fR):
	# prepare X^T (iw), X^T w, X^T W X arrays here
	# TODO - should probably link this (and fR) explicitly to convolution width
	# so that, if that changes, things behave sensibly
	nMin = np.min (nAA)
	nMax = np.max (nAA)
	bSpacing = 2.0 / inst.nPx 
	fAb = np.arange (nMin, nMax, bSpacing)
	#
	nbf = len (fAb)
	npy, npx = iAA.shape
	xtiw = np.zeros (nbf)
	xtwV = np.zeros (nbf)
	xtwx0 = np.zeros (nbf)
	xtwx1 = np.zeros (nbf-1)
	code = """
	double sum_iw = 0.0;
	double sum_w = 0.0;

	double edenom = 0.5 * fR * fR;

	for (int i = 0; i < npy; i++)
	{
		int offset = i*npx;
		double *iA = iAA + offset;
		double *wA = wAA + offset;
		double *fA = nAA + offset;
		bool *cA = cond + offset; // pixel mask (should probably window out bad pixels ...
								// check that that's being done ...)

		// more intelligent choice here?
		int i0 = 0;
		int i1 = 0;

		for (int j = 0; j < npx; j++)
		{
			if (!cA[j]) continue;

			double ij = iA[j];
			double wj = wA[j];

			sum_iw += ij * wj;
			sum_w += wj;

			double fp = fA[j];

			// sort out the basis functions that we're going to deal with ...
			// Need to get indexing right here ...
			while ((fp - fAb[i0] > fR) && i0 < nbf-1) i0++;
			while ((fAb[i1] - fp < fR) && i1 < nbf-1) i1++;

			// loop over basis functions

			double prevG = 0.0;
			for (int ibf = i0; ibf < i1; ibf++)
			{
				double r = fp - fAb[ibf];
				double gj = exp (- r*r / edenom);
				xtiw[ibf] += ij * wj * gj;
				xtwV[ibf] += wj * gj;

				// now, want tridiagonal approx ...
				xtwx0[ibf] += wj * gj * gj;
				if (ibf != 0)
					xtwx1[ibf-1] += wj * gj * prevG;

				prevG = gj;
			}
		}
	}

	double mval = sum_iw/sum_w;
	for (int i = 0; i < nbf; i++)
	{
		xtiw[i] -= mval * xtwV[i];
		xtwV[i] /= sum_w;
	}

	// do probably want to modify X^T w with the funny factor ...
	// easiest to do it once and for all here ...

	return_val = 0;
	"""
	ret = weave.inline (code, ['iAA', 'wAA', 'nAA', 'cond', 'fAb', 'fR', 'npx', 
		'npy', 'nbf', 'xtiw', 'xtwV', 'xtwx0', 'xtwx1'])
	# TODO - check exactly what xtwvx1 means
	return fAb, (xtiw, xtwV, xtwx0, xtwx1)


