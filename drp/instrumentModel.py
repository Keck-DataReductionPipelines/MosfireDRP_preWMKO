# Model of semi-generic instrument
# Robert Lasenby 2009

from math import *
from drpUtils import *
import numpy as np
import scipy.special
import scipy.interpolate

class Instrument ():
	# Converts between [-1,1] normalised and [0, nPx] pixel coordinates
	def px2pxn (self, px):
		return 2.0 * px / self.nPx - 1.0
	def pxn2px (self, pxn):
		return 0.5 * self.nPx * (pxn + 1.0)
	def py2pyn (self, py):
		return 2.0 * py / self.nPy - 1.0
	def pyn2py (self, pyn):
		return 0.5 * self.nPy * (pyn + 1.0)
	
	# return yMin, yMax in field angle units
	def slitEdges (self, slit):
		y0 = slit['slitY']
		yd = 0.5 * slit['slit length']/self.fieldAngle
		return y0 - yd, y0 + yd

# apologies for the pun :)
class Band ():
	pass

# TODO - document exactly what this assumes
def slitConfigFromFits (hdulist):
	# TODO - need to check that we're looking in the right header ...
	td = hdulist[2].data
	# Also, need to put in the slitX, slitY fields
	d = [td.field ('slit number'),
			td.field ('slit width'), td.field ('slit length'),
			td.field ('slitX'), td.field ('slitY')]
	dt = np.dtype ([('slit number', 'i4'), ('slit width', 'f8'),
		('slit length', 'f8'),  ('slitX', 'f8'), ('slitY', 'f8')])
	da = np.rec.fromarrays (d, dtype=dt)
	return da

# Legendre polynomials, orthonormal on [-1,1]
def getPolyBasis (n, m):
	xPolys = [scipy.special.legendre(i) for i in xrange(n+1)]
	yPolys = [scipy.special.legendre(i) for i in xrange(m+1)]
	return xPolys, yPolys

# This classes "implements" a bivariate function
# f : R^2 -> R
#     x,y -> z
# where f(x,y) = \sum_{i, j} \beta_{i, j} f_i (x) g_j (y)
# with (f_1, \dots, f_{m_x}) and (g_1, \dots, g_{m_y})
# sets of basis functions in x and y
# Thus, f is separable in x and y (e.g. the f_i and g_j are
# polynomials)
# TODO - clean this class up a bit
class glmSep2D ():
	# should probably give it the set of basis functions it's going to use ...
	def __init__ (self, fAx, fAy, bh=None):
		self.fAx = fAx
		self.fAy = fAy
		self.mx = len(fAx)
		self.my = len(fAy)
		self.bh = np.zeros (self.mx * self.my) if bh == None else bh

	# get the design matrix
	# TODO : check inputs?
	def gm (self, xA, yA):
		xV = np.array([fx(xA) for fx in self.fAx])
		yV = np.array([fy(yA) for fy in self.fAy])
		# so, next ... need to take outer product of these, in some sense ...
		return (xV[:,np.newaxis,:] * yV[np.newaxis,:]).reshape(self.mx*self.my, -1).transpose()

	# predictions - returns the single-point variances, rather than
	# the whole covariance matrix ...
	def pred (self, gAA):
		errs = np.array([np.dot (x, np.dot (self.bCov, x)) for x in gAA])
		return np.dot (gAA, self.bh), errs

	# quick hack - probably not generally useful ...
	def pred_noErrs (self, gAA):
		return np.dot (gAA, self.bh)
	def __call__ (self, xA, yA):
		gAA = self.gm (xA, yA)
		return np.dot (gAA, self.bh)

	# version not taking into account any prior stuff ...
	def newBh (self, gAA, zA, wA):
		# TODO - make sure we're being consistent wrt whether wA
		# is variance or standard deviation ...
		return np.dot (np.linalg.pinv (wA[:,np.newaxis] * gAA), wA*zA)

	# regularised update, where we have a diagonal observation-covariance
	# matrix, and a prior covariance matrix as a multiple of the identity
	# wA is an array with the 1/\sigma_i in - i.e. the sqrt of what
	# we'd usually call the weight array ...
	# TODO - think about whether we want the new covariance matrix as well ...
	# wb0 should be the squared version ...
	def reg1Bh (self, gAA, zA, wA, wb0):
		# so, let's do the SVD ...
		wAA = wA[:,np.newaxis] * gAA
		zW = wA * zA
		# do we need full matrices?
		U, s, VT = np.linalg.svd (wAA, full_matrices=False)
		sReg = s / (np.square (s) + wb0)
		zM = np.dot (np.transpose (VT), sReg * np.dot (np.transpose (U), zW))
		bhM = np.dot (np.transpose (VT), (s*sReg) * np.dot (VT, self.bh))
		self.bh += zM - bhM

	# So - initialising from a set of points ... no weighting or anything ...
	# Don't bother with covariance matrix ...
	def initFromPts (self, xA, yA, zA):
		gAA = self.gm (xA, yA)
		self.bh = self.newBh (gAA, zA, np.ones_like (zA))

# Does interpolation between bivariate fits
class TrivariateInterp ():
	# bX, bY the sets of basis functions for x and y
	def __init__ (self, bX, bY):
		self.bX = bX
		self.bY = bY

	# a parameter vector specifies each bivariate fit,
	# and we interpolate between the parameter vectors
	# z0, z1 are endpoints of z param interval
	# TODO - think about interpolating uncertainty in
	# bivariate fit ...
	def interpBivariateFits (self, zA, fpA, z0, z1):
		# splprep only interpolates between 11d vectors and lower
		# for the moment, just interpolate component-wise
		# come back later and think about how bad this is,
		# and how we might improve
		fpAt = np.transpose (fpA)
		self.tckA = [scipy.interpolate.splrep (zA, cA,
			xb = z0, xe = z1) for cA in fpAt]
		#self.tck, _ = scipy.interpolate.splprep (x = fpA, u = zA, ub = z0, ue = z1)

	def xyF (self, zA):
		#return splev (zA, self.tck)
		return glmSep2D (self.bX, self.bY, [scipy.interpolate.splev(zA, tck) for tck in self.tckA])

####

# At the moment, this stuff used in the absolute wavelength calibration
def getLegendre1D (m):
	return [scipy.special.legendre(i) for i in xrange(m+1)]

class glm1D ():
	def __init__ (self, fAx):
		self.fAx = fAx

	# single-point update
	# wy = 1/\sigma_y^2
	def updatePoint (self, x, y, wy):
		xA = [f(x) for f in self.fAx]
		zA = np.dot (self.bCov, xA)
		yh = np.dot (self.bh, xA)
		denom = wy * np.dot (xA, zA) + 1.0
		self.bh += (wy * (y-yh)/denom) * zA
		self.bCov -= ((zA / denom) * zA[:,np.newaxis])

	def gAA (self, xA):
		return np.transpose ([f(xA) for f in self.fAx])

	def __call__ (self, x):
		return np.dot (self.bh, [f(x) for f in self.fAx])

	# update ab-initio
	# wA_i = 1/\sigma_i
	def updateFull (self, xA, yA, wA):
		gAA = self.gAA (xA)
		gAA *= wA[:,np.newaxis]
		gAAp = np.linalg.pinv (gAA)
		self.bh = np.dot (gAAp, wA * yA)
		self.bCov = np.dot (gAAp, np.transpose (gAAp))
