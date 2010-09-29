# Model of MOSFIRE instrument, from Zemax raytracing data
# Robert Lasenby 2009

# TODO - should refactor simulator to use this model

from instrumentModel import *
import os

mosfire = Instrument()
mosfire.nPx = 2048
mosfire.nPy = 2048
mosfire.pixSize = 18 # width (and height) of a pixel in microns
mosfire.fieldAngle = 60 * 3.072 # arcsec per unit field angle
mosfire.slitTilt = 4*pi/180.0 # angle of slits relative to detector y axis
mosfire.tanSlitTilt = tan(mosfire.slitTilt) # angle of slits relative to detector y axis
# This approximation doesn't seem to be working as well as it should - TODO : investigate
mosfire.edgeOffset = 0.70910000000000006

mosfire.focalSurfaceScale = 0.7238 # mm per arcsecond
mosfire.barPitch = 8.0 # arcsec : spatial y distance between centres of bars
mosfire.barGap = 0.7 / mosfire.focalSurfaceScale # arcsec : width of blocked-out gap between bars

if os.environ.has_key('MOSFIRE_DATA'):
	path = os.environ['MOSFIRE_DATA']
else:
	path = '../data'

raytraceFile = os.path.join(path,"raytrace-1.0.txt")
#detectorQEFile = os.path.join(path, "MOSFIRE-5_2000nm_GLS4.fits")
detectorQEFile = os.path.join(path, "qe1.fits")
transferFiles = {'K' : os.path.join(path,"K_tp_tot.dat"),
		'H' : os.path.join(path, "H_tp_tot.dat"),
		'J' : os.path.join(path, "J_tp_tot.dat"),
		'Y' : os.path.join(path, "Y_tp_tot.dat")}
skyBgFile = os.path.join(path, 'nearIR_skybg_16_15_stripped.dat')
lineListFile = os.path.join(path, 'list_v2.0.dat')

def processOpticalData (inst, optData):
	# convert from mm to pixels
	optData['px'] /= 0.001 * inst.pixSize
	optData['px'] += inst.nPx/2
	optData['py'] /= 0.001 * inst.pixSize
	optData['py'] += inst.nPy/2

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

def cleanOpticalData (optData):
	# hack to clean up the bad values from the raytrace
	# this is specific to raytrace-1.0.txt
	cond1 = optData['px'] > -30
	cond2 = np.logical_and (optData['slit'] <= 10, optData['py'] < 0)
	cond3 = np.logical_and (optData['slit'] >= 36, optData['py'] > 0)
	cond = np.logical_and (cond1, np.logical_not (np.logical_or (cond2, cond3)))
	return np.extract (cond, optData)

def loadTransferFn (fname):
	"""Load the instrument transfer function from a file, returning
	[(um, fractional efficiency)]"""
	f = open(fname, 'r')
	def readLn(l) : return (float(l[0]), float(l[1]))
	d = [readLn(l.split()) for l in f]
	f.close()
	return d

# TODO - refactor these functions
# Also, want to document where the various files came from ...
def loadSkyBg (fname):
	"""Load the sky background data from a file, returning
	[(um, ph/sec/arcsec^2/nm/m^2)]"""
	f = open(fname, 'r')
	def readLn(l) : return (float(l[0])/1000.0, float(l[1]))
	d = [readLn(l.split()) for l in f]
	f.close()
	return d

def loadLineList (fname):
	f = open (fname, 'r')
	def readLn (l) : return (float(l[0])/10000.0, float(l[1]))
	d = [readLn(l.split()) for l in f if l[0] != '#']
	f.close()
	return d

def getTransferFn (band):
	transfer = loadTransferFn (transferFiles[band.name])
	lA, tA = unzip (transfer)
	# Note - we need to feed this a monotonically increasing
	# lA thing, otherwise it goes silly and gives NaNs (bit stupid ...)
	# In our case, we know that it's in reverse order, so we just reverse
	# things
	lAd, tAd = (np.flipud(lA), np.flipud(tA)) if (lA[0] > lA[1]) else (np.array(lA), np.array(tA))
	tF = scipy.interpolate.InterpolatedUnivariateSpline (lAd, tAd)
	band.transfer = tF
	# TODO - check this this won't upset anything ...
	band.minL = max (band.minL, np.min(lAd))
	band.maxL = min (band.maxL, np.max(lAd))

# band carries around the interpolating functions ...
# have to take into account slit tilts when constructing these ...
def bandFromRaytrace (inst, bandName, optData):
	# Definitely not handling the interpolation etc. in the ideal way ftm ...
	band = Band ()
	band.name = bandName
	condBand = optData['band'] == band.name
	bandL = np.extract (condBand, optData['lambda'])
	band.minL = np.min (bandL)
	band.maxL = np.max (bandL)
	#
	yA = np.unique (optData['yin'])
	xPb, yPb = getPolyBasis (5, 5)
	# yAA[i] holds y values for slit with x0 = xA[i]
	#yAA = (xA - xA[:,np.newaxis])/inst.tanSlitTilt
	# yAA[:,i] holds y values for slit with x0 = xA[i]
	#yAA = (xA[:,np.newaxis] - xA)/inst.tanSlitTilt
	def getFits (y):
		cond = np.logical_and (condBand, optData['yin'] == y)
		# convert into suitable form for interpolation ...
		# remember the details of this conversion!
		px = inst.px2pxn (np.extract(cond, optData['px']))
		py = inst.py2pyn (np.extract(cond, optData['py']))
		x = np.extract (cond, optData['xin'])
		l = np.extract (cond, optData['lambda'])
		bfit1 = glmSep2D (xPb, yPb)
		bfit1.initFromPts (px, x, py)
		bfit2 = glmSep2D (xPb, yPb)
		bfit2.initFromPts (px, x, l)
		return bfit1, bfit2
	pyF, lF = unzip ([getFits(y) for y in yA])
	# So - for each x in xA, go across the yAA[i] row
	# and calculate the relevent values ...
	# TODO - think about this ...
	npx = 100
	pxA = np.linspace(-1.0, 1.0, npx)
	# Also, what range of xA? Go with original ftm ...
	x0A = np.unique (optData['xin'])
	ny = len(yA)
	nx = len(x0A)
	#pxAtile = np.tile (pxA, nx)
	# xAA[:,i] holds x values for slit with x0 = xA[i]
	# Diagram all this out? Probably useful ...
	xAA = x0A + inst.tanSlitTilt * yA[:,np.newaxis]
	def getFitsTilt (fA):
		# TODO - figure out what the numpy alternative
		# to this idiom is
		#vAA = np.array([fA[i](pxA, xAA[i]) for i in xrange(ny)])
		# TODO - document broadcasting sematics for fA etc.
		vAA = np.array([fA[i](pxA[:,np.newaxis], xAA[np.newaxis,i]) for i in xrange(ny)]).reshape((ny, npx, nx))
		#pdb.set_trace()
		fits = [glmSep2D (xPb, yPb) for i in xrange(nx)]
		for i in xrange(nx):
			# initialise for a given x0
			#fits[i].initFromPts (pxA[:,np.newaxis], yA[np.newaxis,:], vAA[:,:,i].ravel())
			coords = np.broadcast_arrays (pxA[:,np.newaxis], yA[np.newaxis,:])
			fits[i].initFromPts (coords[0].ravel(), coords[1].ravel(),
					vAA[:,:,i].transpose().ravel())
			#pdb.set_trace()
		return [fit.bh for fit in fits]
	x0, x1 = -1.0, 1.0
	band.pyF3 = TrivariateInterp (xPb, yPb)
	band.lF3 = TrivariateInterp (xPb, yPb)
	band.pyF3.interpBivariateFits (x0A, getFitsTilt(pyF), x0, x1)
	band.lF3.interpBivariateFits (x0A, getFitsTilt(lF), x0, x1)
	return band


