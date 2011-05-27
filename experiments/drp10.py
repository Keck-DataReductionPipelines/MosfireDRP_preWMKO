

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

DRAW = False
if DRAW:
    pl.figure(1)
    pl.clf()
    pl.subplot(2,2,1)
    lamroi = slice(500,600)
    pl.imshow(mdat[1][yroi,lamroi])
    pl.subplot(2,2,2)
    pl.imshow(fdat[1][yroi,lamroi])
    pl.colorbar()
    pl.subplot(2,2,3)
    pl.imshow(slit[:,lamroi])
    pl.colorbar()
    pl.subplot(2,2,4)
    pl.imshow(lslit[:,lamroi])
    pl.colorbar()

    pl.show()


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

if DRAW:
    pl.figure(2)
    pl.clf()
    pl.scatter(ls,ss,s=10)
    pl.xlim([lslit[:,lamroi].min(), lslit[:,lamroi].max()])
    pl.xlim([17151, 17194])
    pl.ylim([-50,2050])

print "Spline"
tick = time.time()
#bs = II.LSQUnivariateSpline(ls, ss, np.arange(np.min(ls)+5, np.max(ls)-5, .13))
samp = range(2, len(ls)-2,55)
bs = II.LSQUnivariateSpline(ls, ss, ls[samp])
print ".. took {0} s".format(time.time()-tick)
lnew = lslit[10,:]
lnew = np.arange(np.min(ls), np.max(ls), .1)
ynew = bs(lnew)


lr = lslit[10,:]
yr = bs(lr)






if DRAW:
    pl.plot(lnew,ynew,color='r')
    pl.show()
