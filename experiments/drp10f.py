
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
import matplotlib

import pyfits as pf

import MOSFIRE.Options as options
from MOSFIRE import IO, Fit, Wavelength

pl.ion()

path = "/users/npk/desktop/c9_reduce/npk_calib3_q1700_pa_0/"
path = "/scr2/mosfire/c9_npk/npk_calib3_q1700_pa_0/"

mdat = IO.readmosfits(path + "m110323_2718.fits")
fdat = IO.readfits(path + "pixelflat_2d_H.fits")
ldat = IO.readfits(path + "lambda_solution_m110323_2718.fits")

gain = 2.15

dat = mdat[1]/fdat[1] * gain
dat = mdat[1] * gain
lam = ldat[1]

dat[np.logical_not(np.isfinite(dat))] = 0.


yroi=slice(86, 160)
yroi=slice(173, 203)
yroi=slice(1015, 1090)
xroi=slice(0,2048)
slit = dat[yroi,xroi]
ivar = 1/(slit)
ivar[np.logical_not(np.isfinite(ivar))] = 0
lslit = lam[yroi,xroi]
roi = np.abs(lslit - 16944.9) < 20
ivar[roi] = 0

fslit = fdat[1][yroi, xroi]

median_slit = sp.ndimage.median_filter(slit, size=(5,5))

deviations = np.abs(slit-median_slit)/np.sqrt(np.abs(median_slit)+5)
deviations[np.logical_not(np.isfinite(deviations))] = 0


ltick = lslit[10,:]
stick = slit[10,:]


print "Interpolation"
xx = np.arange(slit.shape[1]*1.)
yy = np.arange(slit.shape[0]*1.)

X,Y = np.meshgrid(xx,yy)

ls = lslit.flatten()
ss = slit.flatten()
ys = Y.flatten()

sort = np.argsort(ls)
ls = ls[sort]
ss = ss[sort]
ys = ys[sort]

diff = np.append(np.diff(ls), False)
OK = (diff > 0.001) & (deviations.flatten() < 35)


print "Creating transmission function"

transmission = np.sum(slit*ivar,1) / np.sum(ivar,1)
transmission /= np.median(transmission)
slitpos = np.arange(len(transmission))

pl.figure(2)
pl.clf()


troi = np.abs(lslit.flatten() - 15050) > 5
tran_ss = slit.flatten()
tran_ss[troi] = 0
t2 = np.sum(tran_ss.reshape(slit.shape), 1)
t2 /= np.median(t2)
pp = np.poly1d(np.polyfit(slitpos, t2, 2))

pl.plot(slitpos, t2, 'green')

pl.scatter(slitpos, transmission)
pl.plot(slitpos, pp(slitpos))
print pp

ss = (slit / pp(Y)).flatten()
ss = ss[sort]

knots = np.arange(np.min(ls)+20, np.max(ls)-20, .3)

for i in range(2):
    print "meshgrid created"
    print len(np.where(OK)[0])

    bs = II.splrep(ls[OK], ss[OK], k=3, task=-1, t=knots)


    print "Constructing background subtracted image"
    model = II.splev(lslit.flatten(), bs)
    model = model.reshape(slit.shape)
    model *= pp(Y)

    output = slit - model


    std = np.abs(output)/(np.sqrt(np.abs(slit)))
    tOK = (std < 20).flatten() & np.isfinite(std).flatten()
    OK = OK & tOK[sort]



pl.figure(1)
pl.clf()
pl.plot(ls, ss,'x')
pl.plot(ls[OK], ss[OK], 'o')
knots2 = np.arange(np.min(ls)+20, np.max(ls)-20, .1)
bse = II.splev(knots2, bs)
pl.plot(knots2, bse)
ns = np.sqrt(np.abs(bse) + 5**2) * 3
pl.plot(knots2, bse + ns)
pl.plot(knots2, bse - ns)
pl.ylim([-1000,15000])

vmn = -50
vmx = 100
roi = (ltick > 15000) & (ltick < 15377)
pl.figure(3)
pl.clf()
pl.subplot(6,1,1)
pl.imshow(slit[:, roi], vmin=vmn, vmax=100*vmx)
pl.colorbar()
pl.subplot(6,1,2)
pl.imshow(output[:, roi], vmin=vmn, vmax=vmx, cmap=matplotlib.cm.gray)
pl.colorbar()
pl.subplot(6,1,3)
pl.imshow(output[:, roi]/np.sqrt(np.abs(slit[:, roi])), vmin=-1, vmax=1)
#pl.imshow(deviation[:, roi]/np.sqrt(np.abs(med[:, roi])), vmin=-3, vmax=10)
pl.colorbar()
pl.subplot(6,1,4)
pl.imshow(np.resize(OK, slit.shape)[:,roi], vmin=0, vmax=1)
pl.colorbar()


pl.figure(2)
w = 1/np.sqrt(np.abs(slit.flatten()))
w[sort][np.logical_not(OK)] = 0
w = np.reshape(w, slit.shape)
w[np.logical_not(np.isfinite(w))] = 0
t2 = np.sum(slit*w, 1)/np.sum(w,1)

pl.plot(slitpos, t2/np.median(t2), 'yellow')

pl.figure(3)
pl.subplot(6,1,5)
pl.imshow(w[:,roi], vmin=0, vmax=1)
pl.colorbar()


pl.subplot(6,1,6)
pl.imshow(fslit[:,roi], vmin=.97, vmax=1/.97)
pl.colorbar()

if False:
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


