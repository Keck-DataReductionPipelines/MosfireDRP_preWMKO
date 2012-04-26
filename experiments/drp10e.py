
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

fdat[1][np.logical_not(np.isfinite(fdat[1]))] = 1.0
dat = mdat[1]/fdat[1] * gain
lam = ldat[1]

yroi=slice(748, 822)
yroi=slice(1016, 1091)
yroi=slice(707,736)
yroi=slice(173, 203)
yroi=slice(10, 72) # Double bottom slit
yroi=slice(1014, 1093) # middle slit
yroi=slice(1813, 1887) # double high slit
yroi=slice(1634, 1666) # single slit high
yroi=slice(1903, 2020) # Triple top slit
yroi=slice(86, 160) # double low slit
xroi=slice(0,2048)
slit = dat[yroi,xroi]
lslit = lam[yroi,xroi]

for i in xrange(lslit.shape[0]):
	lslit[i,:] -= 0.006 * i  * 0

med = sp.ndimage.median_filter(slit, size=(5,5))
deviation = np.abs(slit-med)/np.sqrt(np.abs(slit))
deviation[np.logical_not(np.isfinite(deviation))] = 0
OK = (deviation >= 0).flatten()
print slit.shape


ltick = lslit[10,:]
stick = slit[10,:]

print "Interpolation"
xx = np.arange(slit.shape[1]*1.)
yy = np.arange(slit.shape[0]*1.)

X,Y = np.meshgrid(xx,yy)
yy = Y.flatten()
ls = lslit.flatten()
ss = slit.flatten()

lknots = np.arange(np.min(lslit)+5, np.max(lslit)-5, 0.7)

sknots = np.array([12., 13., 14., max(yy)-15., max(yy)-14.,max(yy)-13.])

w = 1/np.sqrt(np.abs(ss))
w[np.logical_not(np.isfinite(w))] = 0.
w += 1e-1

for i in range(3):
    print "Creating Spline"
    bs = II.bisplrep(yy[OK], ls[OK], ss[OK], w=w[OK], kx=1, ky=3, task=-1,
            tx=sknots,
            ty=lknots,
	    quiet=1,
            nxest=len(sknots)+2,
            nyest=len(lknots)+10)

    print "Constructing background subtracted image"
    output = slit.copy()
    model = slit.copy()
    for i in xrange(output.shape[0]):
        model[i,:] = II.bisplev(i,lslit[i,:], bs)

    output -= model
    std = np.abs(output)/np.sqrt(np.abs(model))
    OK = OK & (std < 20).flatten()
    print len(np.where(OK)[0])


pl.figure(1)
pl.clf()
pl.plot(ls, ss, 'x')
pl.plot(ls[OK], ss[OK], '.')
lt = np.arange(lknots.min(), lknots.max(), .7)
pl.plot(lt, II.bisplev(7, lt, bs), 'r')
pl.plot(lt, II.bisplev(13, lt, bs), 'b')
pl.ylim([-5000,25000])

pl.figure(2)
pl.clf()
roi = np.abs(ls-15050.5) < 1
pl.scatter(yy[roi], ss[roi], c=ls[roi])
ypos = np.arange(0,np.max(yy),.1)
pl.plot(ypos, II.bisplev(ypos, 15050., bs))
pl.plot(ypos, II.bisplev(ypos, 15049., bs))
pl.plot(ypos, II.bisplev(ypos, 15051., bs))

vmn = -50
vmx = 100
roi = (ltick > 14600) & (ltick < 15160)
pl.figure(3)
pl.clf()
pl.subplot(4,1,1)
pl.imshow(slit[:, roi], vmin=vmn, vmax=vmx)
pl.colorbar()
pl.subplot(4,1,2)
pl.imshow(output[:, roi], vmin=vmn, vmax=vmx, cmap=matplotlib.cm.gray)
pl.colorbar()
pl.subplot(4,1,3)
pl.imshow(output[:, roi]/np.sqrt(np.abs(med[:, roi])), vmin=-1, vmax=1)
#pl.imshow(deviation[:, roi]/np.sqrt(np.abs(med[:, roi])), vmin=-3, vmax=10)
pl.colorbar()
pl.subplot(4,1,4)
pl.imshow(np.resize(OK, slit.shape)[:,roi], vmin=0, vmax=1)
pl.colorbar()

pl.figure(4)
pl.clf()
for i in xrange(1,50,10):
	roi = (np.abs(ls - 15050.5) < 4) & (yy == i)
	pl.plot(ls[roi], ss[roi])

if True:
    hdu = pf.PrimaryHDU(np.array([output/np.std(np.abs(slit))]))
    try: os.remove("e_std.fits")
    except: pass
    hdu.writeto("e_std.fits")

    hdu = pf.PrimaryHDU(np.array([output]))
    try: os.remove("e_resid.fits")
    except: pass
    hdu.writeto("e_resid.fits")

    hdu = pf.PrimaryHDU(np.array([slit]))
    try: os.remove("e_raw.fits")
    except: pass
    hdu.writeto("e_raw.fits")

    hdu = pf.PrimaryHDU(np.array([slit-output]))
    try: os.remove("e_model.fits")
    except: pass
    hdu.writeto("e_model.fits")

    hdu = pf.PrimaryHDU(np.array([lslit]))
    try: os.remove("e_lambda.fits")
    except: pass
    hdu.writeto("e_lambda.fits")

    mask = np.reshape(OK, slit.shape)
    mask = mask.astype(np.int16)
    hdu = pf.PrimaryHDU(np.array([mask]))
    try: os.remove("e_mask.fits")
    except: pass
    hdu.writeto("e_mask.fits")


