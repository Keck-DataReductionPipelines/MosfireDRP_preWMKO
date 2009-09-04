# Simulation of MOSFIRE output, spectroscopic mode 
# Robert Lasenby 2009

from mosfireSim import *
import sys
import getopt

helpstring = """Simple spectroscopic simulation.
Options:
-b [name] : set band to [name]
-o [file] : set output file to [file]
-r        : randomise slit placement
-t [val]  : set exposure time to [val] seconds
-p        : apply counting noise
-q [file] : apply detector qe from fits file [file]
-m [file] : read in mask parameters from mascgen output file [file]
"""

def main ():
	opts, args = getopt.getopt (sys.argv[1:],"hb:o:rt:pqm:")
	if opts == [] or opts[0][0] == "-h":
		print helpstring
		return
	outName = "output.fits"
	bandName = "K"
	exposureTime = 1
	poissonNoise = False
	useDetectorQE = False
	for o,a in opts:
		if o == "-b": 
			bandName = a
		elif o == "-o": 
			outName = a
		elif o == "-r":
			mosfire.slitX = 0.6 * np.random.rand (mosfire.nBars) - 0.3
		elif o == "-t":
			exposureTime = float(a)
		elif o == "-p":
			poissonNoise = True
		elif o == "-q":
			useDetectorQE = True
			if a != '':
				detectorQEFile = a
		elif o == "-m":
			da = readMascgenOutput (a)
			updateInstFromMascgen (mosfire, da)
	im = exposureTime * calcCountImage(bandName)
	if useDetectorQE:
		qe0 = pyfits.getdata (detectorQEFile)
		qe = np.clip (qe0, 0.0, 1.0)
		im *= qe
	if poissonNoise:
		im = genPoissonVals (im)
	#saveAsFits (im, outName)
	saveAsFitsWithExtensions (mosfire, im, outName)

if __name__ == "__main__":
	main()
