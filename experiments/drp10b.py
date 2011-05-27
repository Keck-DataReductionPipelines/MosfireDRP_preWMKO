
'''

Written around 26 May 2011 by npk.

code is used to test B-splines outlier rejection


Created drp10a.py and documented on pg 53 of drp notebook 1

'''

import numpy as np
import time

import pylab as pl
import scipy as sp
from scipy import interpolate as II
import spline

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
print "Spline"

ltick = lslit[10,:]
stick = slit[10,:]


tick = time.time()
smooth = 13
samp = range(2, len(ls)-2,smooth)
bs = II.LSQUnivariateSpline(ls, ss, ls[samp], k=3)
print ".. took {0} s".format(time.time()-tick)
lnew = lslit[10,:]
lnew = np.arange(np.min(ls), np.max(ls), .1)
ynew = bs(lnew)

pl.clf()
pl.scatter(ltick, stick, s=30, color='yellow')
pl.plot(ls, ss, 'x')
pl.scatter(ls[samp], ss[samp], color='red')
pl.xlim([15403, 15422]),
pl.ylim([-50, 700])
pl.text(15418, 268, "S: {0}".format(smooth))


devs = (ss-bs(ls)) / np.sqrt(np.abs(bs(ls)))

OK = np.abs(devs) < 5

samp = range(2, len(ls[OK])-2,smooth)
bs = II.LSQUnivariateSpline(ls[OK], ss[OK], ls[OK][samp], k=3)
pl.plot(ls, bs(ls), 'purple')


pl.show()

pl.figure(2)
pl.plot(ls, ss)
pl.plot(ls, bs(ls))
pl.ylim([-200, 800])
pl.show()
