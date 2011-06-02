
'''

Written around 26 May 2011 by npk.

code is used to test B-splines outlier rejection with UNIVARIATE Splines


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

mdat = IO.readmosfits(path + "m110323_2738.fits")
fdat = IO.readfits(path + "pixelflat_2d_H.fits")
ldat = IO.readfits(path + "lambda_solution_m110323_2737.fits")

gain = 2.15

dat = mdat[1]/fdat[1] * gain
lam = ldat[1]

yroi=slice(748, 822)
yroi=slice(173, 203)
yroi=slice(707,736)
xroi=slice(0,2048)
slit = dat[yroi,xroi]
lslit = lam[yroi,xroi]

median_slit = sp.ndimage.median_filter(slit, size=(11,3))


ltick = lslit[10,:]
stick = slit[10,:]

print "Creating transmission function"
transmission = np.sum(slit,1)
transmission /= np.median(transmission)
slitpos = np.arange(len(transmission))
pp = np.poly1d(np.polyfit(slitpos, transmission, 3))
print pp

pl.figure(2)
pl.scatter(slitpos, transmission)
pl.plot(slitpos, pp(slitpos))
pl.show()

print "Interpolation"
xx = np.arange(slit.shape[1]*1.)
yy = np.arange(slit.shape[0]*1.)

X,Y = np.meshgrid(xx,yy)

ls = lslit.flatten()
ss = (slit*pp(Y)).flatten()

sort = np.argsort(ls)
ls = ls[sort]
ss = ss[sort]

OK = np.diff(ls) > 0.001
OK = np.append(OK,False)


for i in range(5):
    print "meshgrid created"

    knots = np.arange(np.min(ls)+20, np.max(ls)-20, .08)
    bs = II.splrep(ls[OK], ss[OK], k=3, task=-1, t=knots)


    print "Constructing background subtracted image"
    output = slit.copy()
    model = slit.copy()
    for i in xrange(output.shape[0]):
        output[i,:] -= II.splev(lslit[i,:], bs) * pp(i)
        model[i,:] = II.splev(lslit[i,:], bs) * pp(i)

    std = output/np.sqrt(np.abs(model))
    OK = OK & (std < 10).flatten()
    print len(np.where(OK)[0])



pl.figure(1)
pl.clf()
pl.plot(ltick, stick,'x')
pl.plot(ltick, II.splev(ltick, bs))
pl.show()

if True:
    hdu = pf.PrimaryHDU(np.array([output/np.std(np.abs(slit))]))
    try: os.remove("f_std.fits")
    except: pass
    hdu.writeto("f_std.fits")

    hdu = pf.PrimaryHDU(np.array([output]))
    try: os.remove("f_resid.fits")
    except: pass
    hdu.writeto("f_resid.fits")

    hdu = pf.PrimaryHDU(np.array([slit]))
    try: os.remove("f_raw.fits")
    except: pass
    hdu.writeto("f_raw.fits")

    hdu = pf.PrimaryHDU(np.array([slit-output]))
    try: os.remove("f_model.fits")
    except: pass
    hdu.writeto("f_model.fits")

    hdu = pf.PrimaryHDU(np.array([lslit]))
    try: os.remove("f_lambda.fits")
    except: pass
    hdu.writeto("f_lambda.fits")

    mask = np.reshape(OK, slit.shape)
    mask = mask.astype(np.int16)
    hdu = pf.PrimaryHDU(np.array([mask]))
    try: os.remove("f_mask.fits")
    except: pass
    hdu.writeto("f_mask.fits")


