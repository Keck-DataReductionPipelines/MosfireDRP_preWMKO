# Simple "driver" for calibration routines
# Robert Lasenby 2009

from flatFieldEdgeTracing import *
from mosfireRaytraceModel import *
import pyfits

bandName = "K"
inst = mosfire
optData0 = loadRaytraceData (raytraceFile)
optData1 = cleanOpticalData (optData0)
processOpticalData (inst, optData1)
band = bandFromRaytrace (inst, bandName, optData1)

ffFile = '/home/robert/Projects/mosfire/drp/mosfire/spectroscopicSimulator/session5/flatField.fits'
sciFile = '/home/robert/Projects/mosfire/drp/mosfire/spectroscopicSimulator/session5/0.fits'

ffHduList = pyfits.open (ffFile)
sciHduList = pyfits.open (sciFile)
slitConfig = slitConfigFromFits (ffHduList)
qe0 = pyfits.getdata (detectorQEFile)

pyFA = calibFlat (inst, band, ffHduList, qe0, shape1)

# utility function to write ds9 reg file
def writeRegLines (f, xA, yA, colour):
	for i in xrange(len(xA)-1):
		# Off-by-one correction needed for ds9 pixel coords
		px0 = xA[i]*1024 + 1024 + 1
		px1 = xA[i+1]*1024 + 1024 + 1
		py0 = yA[i]*1024 + 1024 + 1
		py1 = yA[i+1]*1024 + 1024 + 1
		f.write ("line " + str(px0) + " " + str(py0) + " " + str(px1) + " " + str(py1) + " # color=" + colour + "\n")

def drawSlitEdges (f, slit, pyF):
	# TODO - refactor this pattern?
	pbX = inst.px2pxn (np.arange (inst.nPx))
	y0, y1 = inst.slitEdges (slit)
	pbY = pyF (pbX, np.ones_like(pbX) * y0)
	ptY = pyF (pbX, np.ones_like(pbX) * y1)
	writeRegLines (f, pbX, pbY, "green")
	writeRegLines (f, pbX, ptY, "red")

freg = open ("all_slits_calib_test.reg", "w")

for slit, pyF in zip (slitConfig, pyFA):
	drawSlitEdges (freg, slit, pyF)

freg.close()
