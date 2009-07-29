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
mosfire.sinSlitTilt = sin(mosfire.slitTilt) # angle of slits relative to detector y axis
mosfire.slitWidthNominal = 0.7 # arcsec
mosfire.nBars = 46
mosfire.barPitch = 8.0 # arcsec : spatial y distance between centres of bars
mosfire.barGap = 0.5 / mosfire.focalSurfaceScale # arcsec : width of blocked-out gap between bars
mosfire.barAperture = mosfire.barPitch - mosfire.barGap
mosfire.fieldAngle = 60 * 3.072 # arcsec per unit field angle
#mosfire.yFWHM = 0.5 # arcsec : FWHM of point spread function in spatial direction,
				  # spread due to seeing, basically
# FIXME - total hack atm - not based on any physical reasoning or data at all
mosfire.slitFalloffScale = 0.5 * mosfire.pixScale # scale (in arcseconds) of slit edge falloff
mosfire.barOffset = 0.5 * mosfire.barPitch * mosfire.sinSlitTilt / mosfire.fieldAngle

# TODO - should have this as default stuff, I spose ...
mosfire.slitWidth = [mosfire.slitWidthNominal/mosfire.fieldAngle] * mosfire.nBars
# and should also allow to just ask for a long slit, etc.
mosfire.slitX = [slitYPos(mosfire, slit)[2]*mosfire.sinSlitTilt for slit in xrange(mosfire.nBars)]

# figures for Keck
mosfire.mirrorArea = 78.5 # m^2
mosfire.mirrorReflectivity = 0.97 * 0.97
mosfire.intensityConv = mosfire.mirrorArea * mosfire.mirrorReflectivity

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

def cleanOpticalData (optData):
	# hack to clean up the bad values from the raytrace
	# this is specific to raytrace-1.0.txt
	cond1 = optData['px'] > -30
	cond2 = np.logical_and (optData['slit'] <= 10, optData['py'] < 0)
	cond3 = np.logical_and (optData['slit'] >= 36, optData['py'] > 0)
	cond = np.logical_and (cond1, np.logical_not (np.logical_or (cond2, cond3)))
	return np.extract (cond, optData)

# example of black body spectrum generation
def blackBodySpectrumExample (band):
	# spacing same as for the Gemini sky background file,
	# just for consistency
	ffL = np.arange (band.minL, band.maxL, 0.00002)
	# temperature of 1500K (arbitrary value
	# - ask about correct number)
	T = 1500
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
