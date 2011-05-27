
'''

Written around 20 May 2011 by npk.

code is used to test B-splines and knot positions of c9 test data.

Created drp10a.py and documented on pg 49-52 of DRP notebook 1

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

dat = mdat[1]/fdat[1]
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
#bs = II.LSQUnivariateSpline(ls, ss, np.arange(np.min(ls)+5, np.max(ls)-5, .13))

ltick = lslit[10,:]
stick = slit[10,:]

pl.figure(1)
pl.clf()
for i in range(1, 7):
    tick = time.time()
    smooth = (i-1) * 3 + 1
    samp = range(2, len(ls)-2,smooth)
    bs = II.LSQUnivariateSpline(ls, ss, ls[samp], k=3)
    print ".. took {0} s".format(time.time()-tick)
    lnew = lslit[10,:]
    lnew = np.arange(np.min(ls), np.max(ls), .1)
    ynew = bs(lnew)

    pl.subplot(3,2,i)
    pl.scatter(ltick, stick, s=30, color='yellow')
    pl.plot(ls, ss, 'x')
    pl.plot(lnew, ynew, color='red')
    pl.scatter(ls[samp], ss[samp], color='red')
    pl.xlim([15403, 15422]),
    pl.ylim([-50, 400])
    pl.text(15418, 268, "S: {0}".format(smooth))


pl.show()
