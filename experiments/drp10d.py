
'''

Written around 26 May 2011 by npk.

code is used to test B-splines outlier rejection with BIVARIATE Splines


Created drp10a.py and documented on pg 53 of drp notebook 1

'''

import numpy as np
import time
import os

import pylab as pl
import scipy as sp
from scipy import interpolate as II
from scipy import signal
import scipy.ndimage
import spline

import pyfits as pf

import MOSFIRE.Options as options
from MOSFIRE import IO, Fit, Wavelength

path = "/users/npk/desktop/c9_reduce/npk_calib4_q1700_pa_0/"
path = "/scr2/mosfire/c9_npk/npk_calib4_q1700_pa_0/"

mdat = IO.readmosfits(path + "m110323_2737.fits")
fdat = IO.readfits(path + "pixelflat_2d_H.fits")
ldat = IO.readfits(path + "lambda_solution_m110323_2737.fits")

gain = 2.15

dat = mdat[1]/fdat[1] * gain
lam = ldat[1]

yroi=slice(707,736)
xroi=slice(0,2048)
slit = dat[yroi,xroi]
lslit = lam[yroi,xroi]

med = sp.ndimage.median_filter(slit, size=(3,3))
deviation = (slit-med)/np.sqrt(np.abs(slit))
OK = (np.abs(deviation) < 5).flatten()
print slit.shape


ltick = lslit[10,:]
stick = slit[10,:]

print "Interpolation"
xx = np.arange(slit.shape[1]*1.)
yy = np.arange(slit.shape[0]*1.)

X,Y = np.meshgrid(xx,yy)
yy = (Y.flatten())[OK]
ls = (lslit.flatten())[OK]
ss = (slit.flatten())[OK]

print "meshgrid created"
bs = II.bisplrep(yy, ls, ss, kx=1, ky=3, task=-1,
        tx=np.array([1.,2.,max(yy)-2,max(yy)-1]),
        ty=np.arange(np.min(ls), np.max(ls)),
        nxest=10,
        nyest=4000)

pl.figure(1)
pl.clf()
pl.plot(ltick, stick)
pl.plot(ltick, II.bisplev(10, ltick, bs))
pl.show()

if True:
    print "Constructing background subtracted image"
    output = slit.copy()
    for i in xrange(output.shape[0]):
        output[i,:] -= II.bisplev(i,lslit[i,:], bs)

    hdu = pf.PrimaryHDU(np.array([output/np.std(np.abs(slit))]))
    try: os.remove("std.fits")
    except: pass
    hdu.writeto("std.fits")

    hdu = pf.PrimaryHDU(np.array([output]))
    try: os.remove("resid.fits")
    except: pass
    hdu.writeto("resid.fits")

    hdu = pf.PrimaryHDU(np.array([slit]))
    try: os.remove("raw.fits")
    except: pass
    hdu.writeto("raw.fits")
