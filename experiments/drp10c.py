
'''

Written around 26 May 2011 by npk.

code is used to test B-splines outlier rejection using univariate splines.


Created drp10a.py and documented on pg 53 of drp notebook 1

'''

import numpy as np
import time

import pylab as pl
import scipy as sp
from scipy import interpolate as II
from scipy import signal
import spline

import pyfits as pf

import MOSFIRE.Options as options
from MOSFIRE import IO, Fit, Wavelength

path = "/scr2/mosfire/c9_npk/npk_calib4_q1700_pa_0/"

mdat = IO.readmosfits(path + "m110323_2737.fits")
fdat = IO.readfits(path + "pixelflat_2d_H.fits")
ldat = IO.readfits(path + "lambda_solution_m110323_2737.fits")

gain = 2.15

dat = mdat[1]/fdat[1] * gain
lam = ldat[1]

yroi=slice(707,736)
slit = dat[yroi,:]
lslit = lam[yroi,:]


print "Constructing Sampled Spectrum"
ss = []
ls = []
for j in xrange(slit.shape[0]):
    ss.extend(slit[j,:])
    ls.extend(lslit[j,:])

ss = np.array(ss)
ls = np.array(ls)

srt = np.argsort(ls)
ls = ls[srt]
ss = ss[srt]

okl = np.append(np.diff(ls) > 0, False)
OK = (np.isfinite(ls)) & (np.isfinite(ss)) & (okl)

ls = ls[OK]
ss = ss[OK]

ltick = lslit[10,:]
stick = slit[10,:]

print "Rejection"
med_ss = signal.medfilt(ss, 13)
dev = (ss - med_ss) / np.sqrt(np.abs(med_ss))
OK = np.abs(dev) < 5
BAD = np.abs(dev) >= 5

print "Interpolation"
samp = range(5, len(ls[OK])-5, 5)
bs = II.InterpolatedUnivariateSpline(ls[OK], ss[OK], k=3)
bs = II.LSQUnivariateSpline(ls[OK], ss[OK], ls[OK][samp], k=3)

pl.figure(3)
pl.clf()
pl.plot(ls, med_ss)
pl.scatter(ls[OK], ss[OK])
pl.scatter(ls[BAD], ss[BAD], color='red')
pl.plot(ltick, bs(ltick))

pl.show()


pl.clf()
pl.scatter(ltick, stick, s=30, color='yellow')
pl.plot(ls, ss, 'x')
pl.xlim([15403, 15442]),
pl.ylim([-50, 700])
pl.plot(ltick, bs(ltick), 'red')

ltick = lslit[3,:]
stick = slit[3,:]

pl.scatter(ltick, stick, s=30, color='orange')
pl.plot(ltick, bs(ltick), 'red')


pl.show()

pl.figure(2)
pl.plot(ls, ss)
pl.plot(ls, bs(ls))
pl.ylim([-200, 800])
pl.show()


print "Constructing background subtracted image"
output = slit.copy()
for y in xrange(slit.shape[0]):
    output[y,:] -= bs(lslit[y,:])

hdu = pf.PrimaryHDU(np.array([output/np.std(np.abs(slit))]))
hdu.writeto("test.fits")
