# Simulation of MOSFIRE output, spectroscopic mode 
# Robert Lasenby 2009

from simulator import *
import os

mosfire = Instrument()
# dimensions of detector, in pixels
mosfire.nPx = 2048
mosfire.nPy = 2048
mosfire.pixSize = 18 # width (and height) of a pixel in microns
mosfire.focalSurfaceScale = 0.7238 # mm per arcsecond
mosfire.pixScale = 0.18 # arcsec per (18um) pixel
mosfire.slitTilt = 4*pi/180.0 # angle of slits relative to detector y axis
mosfire.tanSlitTilt = tan(mosfire.slitTilt) # angle of slits relative to detector y axis
mosfire.slitWidthNominal = 0.7 # arcsec
mosfire.nBars = 46
mosfire.barPitch = 8.0 # arcsec : spatial y distance between centres of bars
mosfire.barGap = 0.7 / mosfire.focalSurfaceScale # arcsec : width of blocked-out gap between bars
mosfire.barAperture = mosfire.barPitch - mosfire.barGap
mosfire.fieldAngle = 60 * 3.072 # arcsec per unit field angle
mosfire.yFWHM = 0.5 # arcsec : FWHM of point spread function in spatial direction,
				  # spread due to seeing, basically
# FIXME - total hack atm - not based on any physical reasoning or data at all
mosfire.slitFalloffScale = 0.5 * mosfire.pixScale # scale (in arcseconds) of slit edge falloff
mosfire.barOffset = 0.5 * mosfire.barPitch * mosfire.tanSlitTilt / mosfire.fieldAngle

# TODO - should have this as default stuff, I spose ...
mosfire.slitWidth = [mosfire.slitWidthNominal/mosfire.fieldAngle] * mosfire.nBars
# and should also allow to just ask for a long slit, etc.
mosfire.slitX = [slitYPos(mosfire, slit)[2]*mosfire.tanSlitTilt for slit in xrange(mosfire.nBars)]

# Default nonsense here at first
mosfire.targetId = ["No target"]*mosfire.nBars
mosfire.targetPriority = [100]*mosfire.nBars
mosfire.targetYOffset = [0.0]*mosfire.nBars

# figures for Keck
mosfire.mirrorArea = 78.5 # m^2
mosfire.mirrorReflectivity = 0.97 * 0.97
mosfire.intensityConv = mosfire.mirrorArea * mosfire.mirrorReflectivity

# detector properties
mosfire.ePerCount = 2.0 # electrons per detector count

# TODO - check how good an approximation these are ...
# values from communication by ccs
mosfire.anamorphicFactors = {"H" : 1.357, "K" : 1.357, "Y" : 1.335, "J" : 1.335}

if os.environ.has_key('MOSFIRE_DATA'):
	path = os.environ['MOSFIRE_DATA']
else:
	path = 'data'

skyBgFile = os.path.join(path, 'nearIR_skybg_16_15_stripped.dat')
transferFiles = {'K' : os.path.join(path,"K_tp_tot.dat"),
		'H' : os.path.join(path, "H_tp_tot.dat"),
		'J' : os.path.join(path, "J_tp_tot.dat"),
		'Y' : os.path.join(path, "Y_tp_tot.dat")}
raytraceFile = os.path.join(path,"raytrace-1.0.txt")
detectorQEFile = os.path.join(path, "MOSFIRE-5_2000nm_GLS4.fits")

# TODO - slitX and slitY seem to be exchanged here ...
def readMascgenOutput (fname):
	f = open (fname, 'r')
	def readLn(l) : return [int(l[0]), int(l[1]), int(l[2]),
			float(l[3]), int(l[4]), int(l[5]), float(l[6]),
			float(l[7]), float(l[8]), l[9], int(float(l[10])), float(l[11]),
			int(l[12]), int(l[13]), float(l[14]), int(l[15]), int(l[16]), float(l[17]),
			float(l[18]), float(l[19]), float(l[20]), float(l[21])]
	d = [readLn(l.split()) for l in f]
	dt = np.dtype([('slit', 'i4'), ('raH', 'i4'), ('raM', 'i4'),
		('raS', 'f8'), ('decD', 'i4'), ('decM', 'i4'),
		('decS', 'f8'), ('width', 'f8'), ('length', 'f8'),
		('name', 'a8'), ('priority', 'i4'), ('targetY', 'f8'),
		('objRaH', 'i4'), ('objRaM', 'i4'), ('objRaS', 'f8'),
		('objDecD', 'i4'), ('objDecM', 'i4'), ('objDecS', 'f8'),
		('slitX', 'f8'), ('slitY', 'f8'), ('objX', 'f8'),
		('objY', 'f8')])
	da = np.rec.array (d, dtype=dt)
	f.close()
	return da

# returns things in arcsecs 
# TODO - change sin stuff to tan (after asking about the physical measurements etc.)
def slitParamsFromMascgen (inst, da):
	slit = 0
	slitWidth = []
	slitX = []
	# need to find contiguous stuff ...
	for d in da:
		n = int (round ((d['length'] + 2.0*inst.barGap) / inst.barPitch))
		slitWidth += [d['width']/inst.fieldAngle]*n
		# TODO - get this from RaDec stuff
		# Where's the telescope pointing? And how is it aligned wrt Ra 
		# and Dec? Ask about all this ...
		x0 = d['slitY']
		slitX += [(x0 + inst.tanSlitTilt * (i - (n-1)/2.0) * inst.barPitch) / inst.fieldAngle
				for i in xrange(n)]
	#return slitWidth, slitX
	inst.slitWidth = slitWidth
	inst.slitX = slitX

def updateInstFromMascgen (inst, da):
	slitParamsFromMascgen (inst, da)
	inst.targetId = da['name']
	inst.targetPriority = da['priority']
	inst.targetYOffset = da['objY'] - da['slitY']

def cleanOpticalData (optData):
	# hack to clean up the bad values from the raytrace
	# this is specific to raytrace-1.0.txt
	cond1 = optData['px'] > -30
	cond2 = np.logical_and (optData['slit'] <= 10, optData['py'] < 0)
	cond3 = np.logical_and (optData['slit'] >= 36, optData['py'] > 0)
	cond = np.logical_and (cond1, np.logical_not (np.logical_or (cond2, cond3)))
	return np.extract (cond, optData)

# example of black body spectrum generation
def blackBodySpectrumExample (band, T):
	# spacing same as for the Gemini sky background file,
	# just for consistency
	ffL = np.arange (band.minL, band.maxL, 0.00002)
	ffI = blackBodySpectrum (T, ffL)
	# normalise so that max value is 1
	maxI = np.max (ffI)
	return zip (ffL, ffI / maxI)

def calcCountImage (bandName):
	skyBg = loadSkyBg (skyBgFile)
	transfer = loadTransferFn (transferFiles[bandName])
	optData0 = loadRaytraceData (raytraceFile)
	optData1 = cleanOpticalData (optData0)
	processOpticalData (mosfire, optData1)
	band = getBand (mosfire, bandName)
	applyTransferFunction (band, transfer, skyBg)
	computeSlitOverlaps (mosfire)
	computeSlitGroups (mosfire)
	return drawSlits (mosfire, band)

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
	# TODO - this probably wants refactoring ...
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
	c10 = pyfits.Column (name='target id', format='10A', array=inst.targetId )
	c11 = pyfits.Column (name='target priority', format='J', array=inst.targetPriority)
	c12 = pyfits.Column (name='target location', format='D', array=inst.targetYOffset)
	slitHDU = pyfits.new_table([c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12])
	hdulist = pyfits.HDUList ([imHDU, slitHDU])
	hdulist.writeto (fname)
