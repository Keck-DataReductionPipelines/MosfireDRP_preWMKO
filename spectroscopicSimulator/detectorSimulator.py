# Simulation of MOSFIRE output, spectroscopic mode
# Robert Lasenby 2009

from math import *
import numpy as np

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
