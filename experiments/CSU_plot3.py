
import numpy as np
import pylab as pl
import scipy as sp
import scipy.io
import os
from MOSFIRE import CSU
from matplotlib.backends.backend_pdf import PdfPages

from IPython.Shell import IPShellEmbed

start_shell = IPShellEmbed()

def fit_line_with_sigclip(xs, data, i = 0):

        ps = np.polyfit(xs, data, 1)
        pf = np.poly1d(ps)


        residual = np.abs(pf(xs) - data)
        sd = np.std(residual)
        
        ok = np.where(residual < 3.0*sd)[0]

        ps = np.polyfit(xs[ok], data[ok], 1)
        pf = np.poly1d(ps)
        return [pf, ok]

path = "/users/npk/dropbox/mosfire/cooldown 9/csu_meas/" 
proto = "m11031%1.1i_%4.4i.fits.sav.mat"

pl.ion()
rs = []
exclude = [266, 391,392,393,394,395,426]
for i in range(212,554):
        if i in exclude: continue
        if os.path.exists(path + proto % (2, i)):
                ld = scipy.io.loadmat(path + proto % (2, i))
        elif os.path.exists(path + proto % (3, i)):
                ld = scipy.io.loadmat(path + proto % (3, i))
        else:
                print i
                raise Exception("No luck")

        ld["img_num"] = i
        rs.append(ld)



bar_posns = [[] for i in range(92)]

for r in rs:
        for i in range(92):
                assert(r["bars"][i] == i+1)
                p = r["poss_mm"][i]
                d = r["deltas_mm"][i]

                if not np.isfinite(p): continue
                if not np.isfinite(d): continue
                bar_posns[i].append([p, d, r["bars"][i]])



fits = []
ff = []
cnt = 1


rmss = []

c = '''struct barInfo movingBars[] = 
{
'''

s = " Coefficients fit March 15th 2011 by npk\n"
s = " Slopes and RMSs are reported between wavelength\n"
s = "  pixel values of 500 < pix < 1700\n\n"
s = " All units in microns (um)\n"
s += "b#  slope      offset RMS\n"
for bar_vals in bar_posns:
        a = np.array(bar_vals)
        assert(cnt == a[0,2])
        x = a[:,0].ravel()
        y = a[:,1].ravel()

        roi = (x>67) & (x<210)
        [f, ok] = fit_line_with_sigclip(x[roi], y[roi])

        # This is a hack for bar #4
        #if cnt == 4:
                #f[0] = -0.480

        ff.append(f)
        fits.append([f[0],f[1]])


        sd = np.std(y[roi]-f(x[roi]))
        rmss.append(sd)
        s += "%2.0i % 8.7f % 4.3f %3.0f\n" % (cnt, 1-f[1], -f[0], sd*1000)
        c += "        { %2.0i,      0,  % 9.9f,  % 4.3f},\n" % (cnt, 1-f[1], -f[0])
        cnt += 1

s += "       Mean RMS at Keck: %3.0f [um] \n" % (np.mean(rmss) * 1000) 
s += "   Mean RMS at Detector: %1.2f [pix] \n" % (np.mean(rmss)/7.25/0.018) 

c += "};"
try:
        f = open("csu_coefficients.c","w")
        f.writelines(c)
        f.close()
except:
        print "Couldn't write csu_coefficients.c"
try:
        f = open("CSU_coefficients_march_15th_2011.txt", "w")
        f.writelines(s)
        f.close()
except:
        print "Can't write file"


pl.figure(1, figsize=(9,7))
ds = []
means = np.zeros(92)
for bar in range(92):
        bv = np.array(bar_posns[bar])
        pos = bv[:,0]
        delts = bv[:,1]
        r = (pos < 12) | (pos > 250)
        if r.any():
                delts = delts[r]
                if means[bar] == 0:
                        means[bar] = delts.mean()
                delts -= means[bar]
                ds.extend(delts.tolist())

ds = np.array(ds)
pl.clf()
pl.hist(ds*1000,40,color='w')
pl.xlim([-30,30])
pl.text(-20, 2200, r"SD: %3.1f $\mu$" % np.std(ds*1000))
pl.title("Residual measurements for stationary bars placed at 10 mm and 255 mm")

pl.hist(ds*1000,40,color='w')
pl.xlim([-30,30])
pl.xlabel(r"Residual position [$\mu$]")
pl.ylabel(r"$N$")
pdf = PdfPages("CSU_fits_mar_15_2011.pdf")
pl.figure(1).savefig(pdf, format="pdf")


# Measured by hand
groups = [25, 50, 75, 100, 125, 150, 172, 196, 220, 244]

pl.figure(1, figsize=(9,7))
pl.clf()
pl.ylim([49,-3])
pl.xlim([0, 270])
pl.axvline(67,color='black',ls='-.',lw=.5)
pl.axvline(210,color='black',ls='-.',lw=.5)
scale = 200
vlines = []
fits = np.array(fits)
for bar in range(92):
        bv = np.array(bar_posns[bar])
        pos = bv[:,0]
        delts = bv[:,1]
        yshift = (bar % 2) / 1.5
        for group in groups:
                p = np.where((np.isfinite(pos)) & (np.abs(pos - group) < 5))[0]
                position = pos[p].mean()
                delta = delts[p].mean() - ff[bar](position)

                if delta > 0: pl.arrow(position, bar/2.-yshift, delta*scale, 0, color='red')
                else: pl.arrow(position, bar/2.-yshift, delta*scale, 0, color='blue')

                if np.abs(delta*scale) > 60:
                        print "b%2.0i p%6.1i g%4.0i % 5.2f" % (bar, position, group, delta)


pl.text(10,0,"0.1 arcsecond")
pl.arrow(10, 1, .0725*scale, 0, head_width=0.5)
pl.xlabel("CSU Position [mm]")
pl.ylabel("Slit #")
pl.figure(1).savefig(pdf, format="pdf")


fits = np.array(fits)
prevfit = fits[3,0] 
fits[3,0]=0.48
#### FITS Offsets
pl.figure(1, figsize=(9,7))
pl.clf()
odd = np.arange(0, 92, 2)
even = np.arange(1, 92, 2)
pl.plot(fits[odd,0],'^-k')
pl.plot(fits[even,0],'o-k')
pl.title("Slit constant coefficient")
pl.xlabel("Slit #")
pl.ylabel("Bar offset [mm]")
pl.legend(['odd bar', 'even bar'])

pl.figure(1).savefig(pdf, format="pdf")

#FITS Coefficients
pl.figure(1, figsize=(9,7))
pl.clf()
odd = np.arange(0, 92, 2)
even = np.arange(1, 92, 2)

pl.plot(fits[odd,1],'^-k')
pl.plot(fits[even,1],'o-k')
pl.legend(['odd bar', 'even bar'])
pl.title("Slit slope")
pl.xlabel("Slit #")
pl.ylabel("Linear coefficient - 1")

pl.figure(1).savefig(pdf, format="pdf")

fits[3,0] = prevfit
if True:

        i = 0
        for bar_vals in bar_posns:
                i += 1
                a = np.array(bar_vals)
                x = a[:,0].ravel()
                y = a[:,1].ravel()

                roi = (x>67) & (x<210)
                [f, ok] = fit_line_with_sigclip(x[roi], y[roi])
                sd = np.std(y[roi]-f(x[roi]))

                if i % 2 == 1:
                        pl.figure(1, figsize=(9,7))
                        pl.clf()
                        pl.subplot(2,1,1)
                else:
                        pl.subplot(2,1,2)

                ok = (x > 12) & (x < 250)
                pl.plot(x[ok], y[ok] - f(x[ok]), 'x')
                if i % 2 == 0: pl.xlabel("x [mm]")
                pl.ylabel("Delta Position [mm]")
                pl.axvline(67,color='black',ls='-.',lw=.5)
                pl.axvline(210,color='black',ls='-.',lw=.5)
                pl.axhline(0,color='black',lw=.5)
                pl.title("b%2.0i  RMS: %2.0f [micron]  %3.2f [pix] (around best fit)" % (i, sd * 1000, sd/7.25/.018))
                pl.ylim([-.1,.1])
                pl.xlim([0,270])
                pl.text(65, .08, "Delta = % 7.2e x + %3.3f" % (f[1], f[0]))
                if i % 2 == 0:
                        pl.figure(1).savefig(pdf, format="pdf")


pdf.close()

