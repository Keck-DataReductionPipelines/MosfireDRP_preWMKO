# Relative wavelength calibration
# Robert Lasenby 2009

from instrumentModel import *
import scipy.special
import scipy.optimize
import scipy.interpolate
from scipy import weave
from scipy.weave import converters

import pdb

# Need to look into how to distribute this code to platforms
# that don't have C compilers installed?

# Right - let's get all the preparation code ready ...
# We're given the image, the inferred p_x, y -> p_y
# function, the "nominal" p_x, y -> \lambda
# function, and our task is to optimise
# the p_x, p_y -> \eta function, where \eta is as before.
# So, let's set all this up ...

def relativeWavelengthCalibrateSlit (inst, band, slit, im, imW, pyF):
	# So ... first task, I spose, is to get the first-guess at
	# p_x, p_y -> \eta
	# Getting \lambda -> \eta is pretty much the same as the prepareSpectrum
	# function in absoluteWavelengthCalibration - might want to refactor that
	# 
	# Right - we're going to set this up so that
	# we have f* : p_x, p_y -> \eta
	# with f* (p_x, p_y) = f_0 (p_x, p_y) + f_\beta (p_x, p_y)
	# where f_0 is some first-guess, and f_\beta is separable
	# polynomial thing.
	# The idea here is that we choose the p_y -wise basis functions
	# for f_\beta such that f_\beta (p_x_0, p_y) = 0 for all \beta,
	# where p_x_0 is some chosen "central" row of pixels.
	#
	# So, our first task is, given p_x, y -> \lambda and
	# p_x, y -> p_y, to get a first guess at p_x, p_y -> \eta
	# Probably easiest to get this in the form of an image ...
	# which frees up our choice of methods somewhat ...
	#
	# Plan, then - obtain a p_x, p_y -> \lambda image,
	# then pick the correct row of this, obtain the
	# \lambda -> \eta function (in a similar way to the prepareSpectrum
	# function in absoluteWavelengthCalibration) and
	# use this to get our p_x, p_y -> \eta image
	py0, iAA, wAA, nAA, cond, h, hInv = prepareImage (inst, band, slit, im, imW, pyF)
	# Don't actually need separate image and weight things (with the approximation
	# that we're using ...)
	# FIXME - why does't iwAA = iAA * wAA work properly? Gives odd results ...
	#iwAA = iAA * wAA
	iwAA = iAA
	# Don't bother optimising stuff atm - just pass right through ...
	shiftAA = optimiseShift (inst, iwAA, nAA, cond)
	#shiftAA = 0.0 # FIXME - pass straight through
	return py0, (iAA, wAA, nAA + shiftAA, cond), h, hInv

# Returns negative of proper value in order to interface properly to optimisation code ...

# How bad is python function call overhead?
# Live with it for the moment ...
# Don't actually need separate weight image quite yet ...
def getLL (iwAA, nAA, cond, fAb, fR, fxA, fyA, b):
	npy, npx = iwAA.shape
	mx = len(fxA); my = len(fyA);
	nbf = len(fAb)
	# How is beta arranged? (x, y) flattened ...
	code = """
	// Question - how inefficient is it to do this allocation every time?
	// Could share the storage by passing as a parameter ... we'll see whether
	// this seems worth it ...
	double *sxv = new double[nbf];
	memset (sxv, 0, nbf * sizeof(double));
	double edenom = 0.5 * fR * fR;

	for (int i = 0; i < npy; i++)
	{
		// calculational nonsense of the correct form here, given
		// that we're windowing the correct section out in some way ...
		int offset = i*npx;
		double *ia = iwAA + offset;
		double *fa = nAA + offset;
		bool *ca = cond + offset; // pixel mask (should probably window out bad pixels ...
								// check that that's being done ...)

		// more intelligent choice here?
		int i0 = 0;
		int i1 = 0;

		for (int j = 0; j < npx; j++)
		{
			if (!ca[j]) continue;

			// evaluate the map
			double fp = fa[j];
			for (int kx = 0; kx < mx; kx++)
			{
				for (int ky = 0; ky < my; ky++)
				{
					fp += b[ky + my*kx] * fyA[ky*npy + i] * fxA[kx*npx + j];
				}
			}

			//if (isnan (fp))
				//printf (\"%i %i %f\\n\", i, j, fa[j]);

			// sort out the basis functions that we're going to deal with ...
			// Need to get indexing right here ...
			while ((fp - fAb[i0] > fR) && i0 < nbf-1) i0++;
			while ((fAb[i1] - fp < fR) && i1 < nbf-1) i1++;

			// loop over basis functions
			for (int ibf = i0; ibf < i1; ibf++)
			{
				bool notAtFirst = !isnan (sxv[ibf]);
				double valFirst = sxv[ibf];

				double r = fp - fAb[ibf];
				double gj = exp (- r*r / edenom);
				sxv[ibf] += ia[j] * gj;
				
				if (isnan(sxv[ibf]) && notAtFirst)
				{
					//printf (\"%i %f %i %f %f %f %f %f\\n\", ibf, valFirst, j, ia[j], gj, r, fAb[ibf], fp);
				}
			}
		}
	}

	/*
	for (int i = 0; i < 20; i++)
		printf (\"%f, ", sxv[i]);
	printf (\"\\n\");
	*/

	// Now, sum of squares and all that ...
	double tot = 0.0;
	for (int i = 0; i < nbf; i++)
	{
		tot += sxv[i]*sxv[i];
	}
	delete[] sxv;
	return_val = -tot;
	"""
	return weave.inline (code, ['iwAA', 'nAA', 'cond', 'fAb', 'fR', 'npx', 
		'b', 'npy', 'nbf', 'fxA', 'fyA', 'mx', 'my'])

def getLLp (iwAA, nAA, cond, fAb, fR, fxA, fyA, b):
	npy, npx = iwAA.shape
	mx = len(fxA); my = len(fyA);
	nbf = len(fAb)
	gradient = np.zeros (mx*my)
	code = """
	// Question - how inefficient is it to do this allocation every time?
	// Could share the storage by passing as a parameter ... we'll see whether
	// this seems worth it ...
	int mb = mx*my;
	double *sxv = new double[nbf];
	memset (sxv, 0, nbf * sizeof(double));

	double *sxbv = new double[nbf*mb]; // holds gradient stuff
	memset (sxbv, 0, nbf * mb * sizeof(double));

	double edenom = 0.5 * fR * fR;

	for (int i = 0; i < npy; i++)
	{
		// calculational nonsense of the correct form here, given
		// that we're windowing the correct section out in some way ...
		int offset = i*npx;
		double *ia = iwAA + offset;
		double *fa = nAA + offset;
		bool *ca = cond + offset; // pixel mask (should probably window out bad pixels ...
								// check that that's being done ...)

		// more intelligent choice here?
		int i0 = 0;
		int i1 = 0;

		for (int j = 0; j < npx; j++)
		{
			if (!ca[j]) continue;

			// evaluate the map
			double fp = fa[j];
			for (int kx = 0; kx < mx; kx++)
			{
				for (int ky = 0; ky < my; ky++)
				{
					fp += b[ky + my*kx] * fyA[ky*npy + i] * fxA[kx*npx + j];
				}
			}

			// sort out the basis functions that we're going to deal with ...
			// Need to get indexing right here ...
			while ((fp - fAb[i0] > fR) && i0 < nbf-1) i0++;
			while ((fAb[i1] - fp < fR) && i1 < nbf-1) i1++;

			// loop over basis functions
			for (int ibf = i0; ibf < i1; ibf++)
			{
				bool notAtFirst = !isnan (sxv[ibf]);
				double valFirst = sxv[ibf];

				double r = fp - fAb[ibf];
				double gj = exp (- r*r / edenom);
				sxv[ibf] += ia[j] * gj;

				// gradient calculation
				// loop through all the coefficients
				for (int ix = 0; ix < mx; ix++)
				{
					for (int iy = 0; iy < my; iy++)
					{
						int ib = ix*my + iy;
						// we accumulate into this basis function's slot for this coefficient
						sxbv[mb*ibf + ib] += - 2.0 * ia[j] * (r/edenom) * 
							gj * fyA[i + iy*npy] * fxA[j + npx*ix];
					}
				}
			}
		}
	}

	// Now, we need to accumulate all the gradient stuff ...
	// act by modifying gradient

	for (int i = 0; i < mb; i++)
	{
		for (int j = 0; j < nbf; j++)
		{
			gradient[i] += sxv[j] * sxbv[j*mb + i];
		}
		gradient[i] *= -2.0;
	}

	delete[] sxv;
	delete[] sxbv;

	return_val = 0;
	"""
 	ret = weave.inline (code, ['iwAA', 'nAA', 'cond', 'fAb', 'fR', 'npx', 
		'b', 'npy', 'nbf', 'fxA', 'fyA', 'mx', 'my', 'gradient'])
	return gradient

# Getting the Hessian (assess how much better this is ... seems to be worth it)
def getLLH (iwAA, nAA, cond, fAb, fR, fxA, fyA, b):
	npy, npx = iwAA.shape
	mx = len(fxA); my = len(fyA);
	nbf = len(fAb)
	hessian = np.zeros (((mx*my),(mx*my)))
	code = """
	// Question - how inefficient is it to do this allocation every time?
	// Could share the storage by passing as a parameter ... we'll see whether
	// this seems worth it ...
	int mb = mx*my;
	double *sxv = new double[nbf];
	memset (sxv, 0, nbf * sizeof(double));

	double *sxbv = new double[nbf*mb]; // holds gradient stuff
	memset (sxbv, 0, nbf * mb * sizeof(double));

	// storage for Hessian stuff ... hmmm ... starts to get a little
	// silly at this stage ... does it? Well, maybe not so much ... let's try :)
	double *sxbhv = new double[nbf*mb*mb];
	memset (sxbhv, 0, nbf*mb*mb*sizeof(double));

	double edenom = 0.5 * fR * fR;

	for (int i = 0; i < npy; i++)
	{
		// calculational nonsense of the correct form here, given
		// that we're windowing the correct section out in some way ...
		int offset = i*npx;
		double *ia = iwAA + offset;
		double *fa = nAA + offset;
		bool *ca = cond + offset; // pixel mask (should probably window out bad pixels ...
								// check that that's being done ...)

		// more intelligent choice here?
		int i0 = 0;
		int i1 = 0;

		for (int j = 0; j < npx; j++)
		{
			if (!ca[j]) continue;

			// evaluate the map
			double fp = fa[j];
			for (int kx = 0; kx < mx; kx++)
			{
				for (int ky = 0; ky < my; ky++)
				{
					fp += b[ky + my*kx] * fyA[ky*npy + i] * fxA[kx*npx + j];
				}
			}

			// sort out the basis functions that we're going to deal with ...
			// Need to get indexing right here ...
			while ((fp - fAb[i0] > fR) && i0 < nbf-1) i0++;
			while ((fAb[i1] - fp < fR) && i1 < nbf-1) i1++;

			// loop over basis functions
			for (int ibf = i0; ibf < i1; ibf++)
			{
				bool notAtFirst = !isnan (sxv[ibf]);
				double valFirst = sxv[ibf];

				double r = fp - fAb[ibf];
				double gj = exp (- r*r / edenom);
				sxv[ibf] += ia[j] * gj;

				// gradient calculation
				// loop through all the coefficients
				for (int ix = 0; ix < mx; ix++)
				{
					for (int iy = 0; iy < my; iy++)
					{
						int ib = ix*my + iy;
						double val = ia[j] * gj * fyA[i + iy*npy] * fxA[j + npx*ix];
						// we accumulate into this basis function's slot for this coefficient
						sxbv[mb*ibf + ib] += - 2.0 * val * (r/edenom);

						// and now we do the whole Hessian thing ...
						for (int ix1 = ix; ix1 < mx; ix1++)
						{
							int iyStart = ix1 == ix ? iy : 0;
							for (int iy1 = iyStart; iy1 < my; iy1++)
							{
								sxbhv[mb*mb*ibf + ib*mb + ix1*my + iy1] += val *
								(4.0 * (r/edenom) * (r/edenom) - 2.0/edenom)
								* fxA [j + npx*ix1] * fyA[i + npy*iy1];
							}
						}
					}
				}
			}
		}
	}

	// TODO - we're leaving a lot of the sxbhv stuff unused, which is
	// rather unnecessary - change that

	for (int i1 = 0; i1 < mb; i1++)
	{
		for (int i2 = i1; i2 < mb; i2++)
		{
			double *hv = hessian + i1*mb + i2;

			for (int j = 0; j < nbf; j++)
			{
				// single derivs ...
				*hv += sxbv[j*mb + i1] * sxbv[j*mb + i2];
				// now we need to put in double derivs ...
				*hv += sxv[j] * sxbhv [j*mb*mb + i1*mb + i2];
			}
			
			*hv *= -2.0;

			// symmetric, so also want to put value in other half of matrix ...
			if (i1 != i2)
				*(hessian + i2*mb + i1) = *hv;
		}
	}

	delete[] sxv;
	delete[] sxbv;
	delete[] sxbhv;

	return_val = 0;
	"""
 	ret = weave.inline (code, ['iwAA', 'nAA', 'cond', 'fAb', 'fR', 'npx', 
		'b', 'npy', 'nbf', 'fxA', 'fyA', 'mx', 'my', 'hessian'])
	return hessian

def optimiseShift (inst, iwAA, nAA, cond):
	# So ... let's see ... our first task, I spose, is to construct
	# our funny separable functions ... alright then ...
	# Actually, we can do this bit once and for all, I spose ...
	mx = 5
	pxPolys = [scipy.special.legendre(i) for i in xrange(mx)]
	# Question - do we want to vary the order here depending on how 
	# big the slit is that we're looking at? Maybe not ... 5 is quite low
	# anyway, I spose ...
	my = 4
	# Don't want the constant term here, since we're basically removing that dof
	pyPolys = [scipy.special.legendre(i) for i in xrange(1,my+1)]
	for i in xrange(my):
		# subtract the value at zero ...
		pyPolys[i] = np.polysub (pyPolys[i], np.poly1d(pyPolys[i](0)))
	# Okay - we've got our basis things, then - so, let's evaluate
	# them :)
	# Let's just do things like this to start with, I spose ...
	# can see how this fits in with everything else later ...
	# Also, do we want arange here or not? Think about that as well ...
	pxA = ((2.0/float(len(iwAA[0]))) * np.arange(len(iwAA[0]))) - 1.0
	# Hmmm ... will this be sufficient? Or, with things not being exactly zero,
	# will we have funny results? Well ... hmmm ... let's see ...
	# intuitively, it should behave pretty much the same ... ?
	pyA = ((2.0/float(len(iwAA))) * np.arange(len(iwAA))) - 1.0
	# Now, evaluate the basis functions here ...
	fpxA = np.array ([f(pxA) for f in pxPolys])
	fpyA = np.array ([f(pyA) for f in pyPolys])
	# Let's get a few things set up here ...
	nMin = np.min (nAA)
	nMax = np.max (nAA)
	pSpacing = 2.0 / inst.nPx
	# For the moment, we assume that our first guess wasn't more than a few
	# pixels off ...
	fAb = np.arange (nMin - 3*pSpacing, nMax + 4*pSpacing, pSpacing)
	#nbf = len(fAb)
	fR = 3.0 * pSpacing # TODO - check that this choice is a good one ...
	# initial param vector - nice and easy :)
	b0 = np.zeros (mx * my)
	#npy, npx = iwAA.shape
	#pdb.set_trace()
	#
	# Okay - the functions themselves ...
	def f (b):
		return getLL (iwAA, nAA, cond, fAb, fR, fpxA, fpyA, b)
	def fprime (b):
		return getLLp (iwAA, nAA, cond, fAb, fR, fpxA, fpyA, b)
	def fhess (b):
		return getLLH (iwAA, nAA, cond, fAb, fR, fpxA, fpyA, b)
	#
	# Do we want to bother with a prior term?
	#bOpt = scipy.optimize.fmin_ncg (f, b0, fprime = fprime, fhess = fhess, avextol = 1e-6)
	bOpt = scipy.optimize.fmin_ncg (f, b0, fprime = fprime, fhess = fhess)
	fShift = glmSep2D (pxPolys, pyPolys, bOpt)
	pxAA, pyAA = np.meshgrid (pxA, pyA)
	shiftAA = fShift (np.ravel (pxAA), np.ravel (pyAA)).reshape (nAA.shape)
	return shiftAA

def resample1D (oldx, y, newx):
	iF = scipy.interpolate.InterpolatedUnivariateSpline (oldx, y)
	return iF (newx)

# TODO - this host lots of similarities to prepareSpectrum function ...
# want to refactor all this out properly ...
def prepareImage (inst, band, slit, iAA, wAA, pyF):
	# TODO - refactor this pattern?
	y0, y1 = inst.slitEdges (slit)
	slitX = slit['slitX']
	slitY = slit['slitY']
	slitX0 = slitX - inst.tanSlitTilt * slitY
	# TODO - probably don't want to use the "band" values ... want
	# to pass things through to replace these ...
	#pyF = band.pyF3.xyF (slitX0)
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
	# Okay ... so, that gives us our \lambda image ...
	# Now, one thing we need to do is make sure that the \lambda -> \eta function
	# behaves sensibly outside the pixel range of the row that we're looking at ...
	# hmmm ... thoughts? Seem to remember that things behave sensibly ...
	# might want to check on that, of course ... well, do this and see what happens
	lA0 = lAA[(py1-py0)/2]
	# could make this more sparse ... useful?
	# Depends how often we'll be using it, I spose ... let's make it
	# dense for now, think about sparsity later ...
	# Remember that \lambda goes in the opposite direction, so need to
	# flip things ...
	fn = scipy.interpolate.InterpolatedUnivariateSpline (pxCn, lA0)
	fInv = scipy.interpolate.InterpolatedUnivariateSpline (np.flipud (lA0), np.flipud(pxCn))
	#
	nAA = fInv(np.ravel(lAA)).reshape ((len(pyC), inst.nPx))
	# Okay - we've got our first-guess, then
	# Now, next step is to get our conditional mask, I spose ...
	# again, we definitely want to share all our logic here with the later stages ...
	iAAw = iAA[py0:py1+1]
	wAAw = wAA[py0:py1+1]
	# Do we want to cond out low-weight pixels for efficiency reasons
	# make sure that we're properly within the slit ...
	cond1 = np.logical_and (yAA > y0 + 2 * (2.0 / inst.nPx), yAA < y1 - 2 * (2.0 / inst.nPx))
	return py0, iAAw, wAAw, nAA, cond1, fn, fInv
