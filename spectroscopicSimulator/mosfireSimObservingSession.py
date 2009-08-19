# Given some mask design inputs, simulates the files that would
# be produced by an observing session
# Robert Lasenby 2009

from mosfireSim import *
import sys
import getopt

import pdb

# Right - to start with, let's just simulate the spectrum
# change with time ... something fairly simple ...

# percentage change in spectral intensity at wavelength l
# at time t, relative to t=0, for sky background
# Current implementation needs to be swapped out for one
# based on some actual data ...
# t in seconds, l in um
def iChange (t, l):
	return 1 + 0.01 * sin(t/60.0) * (l - 1.0)

# TODO - refactor this stuff
def getCounts (imPhotIntegratedFlux, qe, eGain, readNoise):
	# TODO - check whether we want to use distribution approximations
	imPhotCount = np.random.poisson (imPhotIntegratedFlux)
	# for pixel with QE = p, the number of electrons has a
	# a binomial distribution Bin(n_phot, p)
	# Since np.random.binomial doesn't like n=0,
	# we need to be a little tricky here ...
	#imECount = np.random.binomial (imPhotCount, qe)
	qe1 = np.where (imPhotCount == 0, 0.0, qe)
	photCount1 = np.maximum (imPhotCount, 1)
	imECount = np.random.binomial (photCount1, qe1)
	# We need to make some assumption about the read noise here ...
	# Inverse gain appears to operate basically like a truncation,
	# but then we have the noise stuff ... don't have wonderful
	# data on this, but I spose we go with something like ...
	# Can things go negative? I spose they might be able to ... not that great
	# if they do, of course, but anyway ... just go with a normal
	# distribution, then ...
	# Go with same noise for every pixel atm ...
	imCount = np.rint( np.random.normal (imECount*eGain, readNoise))
	return imCount

# use Fowler-N sampling, producing an inferred variance image
def getCountsFowler (imPhotIntegratedFlux, qe, eGain, readNoise, fowlerN):
	imPhotCount = np.random.poisson (imPhotIntegratedFlux)
	qe1 = np.where (imPhotCount == 0, 0.0, qe)
	photCount1 = np.maximum (imPhotCount, 1)
	imECount = np.random.binomial (photCount1, qe1)
	# Should really add together fowlerN things (for rounding properties),
	# but this'll do ftm
	imCount = np.rint( np.random.normal (fowlerN*imECount*eGain, sqrt(fowlerN)*readNoise)) / float(fowlerN)
	# inferred inverse-variance image
	# gives no-prior value
	#imWeight = (fowlerN-1)/(readNoise * readNoise * np.random.chisquare (fowlerN-1, size=imCount.shape))
	# totally prior:
	imWeight = np.tile (1.0 / (readNoise * readNoise), imCount.shape)
	return imCount, imWeight

class Exposure ():
	def __init__ (self, name, time, startTime=0.0, nodX=0.0, nodY=0.0):
		self.name = name
		self.time = time
		self.nodX = nodX
		self.nodY = nodY

def main ():
	# so - how to deal with the options here?
	# need to specify the sky spectrum, etc. etc.
	opts, args = getopt.getopt (sys.argv[1:], "d:m:bf:s:")
	outputdir = "session1"
	skyBgIm = None
	flatFieldIm = None
	outputBaseImages = False
	for o, a in opts:
		if o == "-m":
			# get input from mask configuration file
			da = readMascgenOutput (a)
			updateInstFromMascgen (mosfire, da)
		if o == "-d":
			# set output directory
			outputdir = a
		if o == "-b":
			outputBaseImages = True
		if o == "-s":
			skyBgIm = pyfits.getdata (a)
		if o == "-f":
			flatFieldIm = pyfits.getdata (a)
	if not os.path.exists (outputdir):
		os.mkdir (outputdir)
	#
	bandName = "K"
	skyBg = loadSkyBg (skyBgFile)
	transfer = loadTransferFn (transferFiles[bandName])
	optData0 = loadRaytraceData (raytraceFile)
	optData1 = cleanOpticalData (optData0)
	processOpticalData (mosfire, optData1)
	band = getBand (mosfire, bandName)
	applyTransferFunction (band, transfer, skyBg)
	skySI = band.spectralIntensity
	computeSlitOverlaps (mosfire)
	computeSlitGroups (mosfire)
	# Firstly, get reference sky background image
	if skyBgIm == None: skyBgIm = drawSlits (mosfire, band)
	# Also want the p_x, p_y -> lambda map for this slit config
	# ... won't have values in between the slits, but that doesn't
	# really matter, hopefully ...
	#lAA = blah
	# Detector qe
	qe0 = pyfits.getdata (detectorQEFile)
	qe = np.clip (qe0, 0.0, 1.0)
	#Firstly, then, the flat field
	flatLampSpec = blackBodySpectrumExample (band, 2000)
	applyTransferFunction (band, transfer, flatLampSpec)
	ffSI = band.spectralIntensity
	if flatFieldIm == None : flatFieldIm = drawSlits (mosfire, band)
	if outputBaseImages:
		saveAsFits (skyBgIm, "skyBg1.fits")
		saveAsFits (flatFieldIm, "ff1.fits")
	#flatFieldCounts = getCounts (flatFieldIm*1000.0, qe, 1.0/mosfire.ePerCount, 1.0)
	# NB - readNoise parameter is the standard deviation
	# What's a sensible value for this? Should probably put
	# into mosfireSim.py as part of the instrument params
	flatFieldCounts, flatFieldWeights = getCountsFowler (flatFieldIm*1000.0, qe, 1.0/mosfire.ePerCount, 10.0, 8)
	# need to see what plan is wrt flat field extensions / headers etc.
	saveAsFitsWithExtensions (mosfire, flatFieldCounts,
			os.path.join(outputdir, "flatField.fits"), [flatFieldWeights])
	# How's the flat field stuff going to be done -
	# are we going to just get one low-noise flat field,
	# or are we expected to combine multiple noisy
	# flat fields into one good one? Ask ...
	# well, it probably depends on the whole approach,
	# I spose ... maybe we have a flat field for every exposure?
	# Is that the plan?
	#
	# Wrt providing a line list (and for arc lamp exposures ...)
	# - how should that work?
	#
	# Then, for each exposure, we want to apply
	# the necessary corrections to this, and
	# add the necessary object spectra on top ...
	# Basically, then, we need a list of exposures
	# and their properties ...
	# At the moment, we seem to be assuming that we're
	# just getting a sequence of fairly normal-looking
	# exposures ... correct? Let's see ...
	# If we're dealing with only a few nod positions, probably
	# want to pre-render the frames? 
	# TODO - deal with this properly
	nExposures = 1
	exposures = [Exposure(name=str(i),time=60.0) for i in xrange(nExposures)]
	for i in xrange(nExposures):
		exposures[i].nodY = 0.25*mosfire.barAperture*(1 - (i%2))
		exposures[i].nodX = mosfire.tanSlitTilt * exposures[i].nodY
	for e in exposures:
		# form the object images ...
		imPhot = skyBgIm
		# Okay ... something like ... well, we need to find
		# the correct y coordinate for the object (assume
		# that we keep it suitably centred on slit)
		# In terms of which slit the object goes in
		# for the case of multiple slits, we want to use
		# the centre of the combined slit
		# for the pointObjYSamp function, and can use
		# any slit for pxYImage stuff. Okay ... let's do this ...
		# add in the object images
		nTargets = len(mosfire.targetYOffset)
		for i in xrange (nTargets):
			print i
			#if i >= 35 : pdb.set_trace()
			slitN = mosfire.slitGroups[i][0]
			# TODO - refactor 0.7 and 20 out
			# Also, should really obtain x from our target list ...
			#xC = np.mean ([mosfire.slitX[g] for g in mosfire.slitGroups[i]])
			xC = mosfire.targetX[i]/mosfire.fieldAngle
			yC = mosfire.targetY[i]/mosfire.fieldAngle
			yA, eA, xA, eAx = pointObjYSamp (mosfire, xC, yC, 0.7, 20)
			iAA = pxYImage (mosfire, band, slitN, xA, yA)*eA[:,np.newaxis]
			pyA, iAAd = distortY (mosfire, band, xA, yA, iAA)
			#pdb.set_trace()
			imPhot[pyA] += iAAd
		#saveAsFits (imPhot, e.name + ".fits")
		#imCount = getCounts (imPhot*e.time, qe, 1.0/mosfire.ePerCount, 10.0)
		imCount, imWeights = getCountsFowler (imPhot*e.time, qe, 1.0/mosfire.ePerCount, 10.0, 8)
		saveAsFitsWithExtensions (mosfire, imCount, 
				os.path.join(outputdir, e.name + ".fits"), [imWeights])
	# Now, outputting the plan file ... using configparser stuff ...
	# well, to start with, just write plain text stuff ...
	planfileName = os.path.join (outputdir, "session.plan")
	planfile = open (planfileName, 'w')
	planfile.write ("# Files generated by mosfireSimObservingSession.py, revision alpha\n")
	planfile.write ("[files]\n")
	planfile.write ("FLATNAME : flatField.fits\n")
	planfile.write ("SCIENCENAME : ")
	for e in exposures[:-1] : planfile.write (e.name + ".fits , ")
	planfile.write (exposures[-1].name + ".fits\n")
	planfile.close()
	# Okay ... how about cosmic ray hits as well? Well ... leave that
	# alone for the moment?

if __name__ == "__main__":
	main()
