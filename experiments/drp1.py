'''

DRP Experiment code.

npk March 23rd 2011

'''
import MOSFIRE
from MOSFIRE import Fit, IO
import numpy as np, pylab as pl

pl.ion()
(header, data, targs, ssl, msl, asl) = IO.readfits_all("/users/npk/desktop/c9/m110323_2690.fits")

reload(Fit)
reload(IO)

(x,y) = (936, 85) # Central pixels picked by eye
(x,y) = (936, 1807) # Central pixels picked by eye
(x,y) = (936, 1006) # Central pixels picked by eye
(x,y) = (936, 85+43.532*10) # Central pixels picked by eye

pl.figure(1)
pl.clf()

roi_width = 15 

yposs = []
widths = []
xposs = []

for i in range(-49, 50):
        delt = 1900/100.
        xp = x+delt*i
        v = data[y-roi_width:y+roi_width, xp-2:xp+2]
        v = np.median(v, axis=1)
        pl.plot(v, '*')
        pl.ylabel("slice")

        xposs.append(x+delt*i)


        ff = Fit.do_fit(v, residual_fun=Fit.residual_disjoint_pair)
        if (0 < ff[4] < 4):
                xs = np.arange(len(v))
                pl.plot(xs, Fit.fit_disjoint_pair(ff[0], xs))

                yposs.append(ff[0][1] - roi_width)
                widths.append(ff[0][5])
                #print "%4.1f" % ff[0][1]
        else:
                print "Skipping: %i" % (x+delt*i)

(xposs, yposs, widths) = map(np.array, (xposs, yposs, widths))
i = 6
fun = np.poly1d(Fit.polyfit_clip(xposs, yposs, i))
wfun = np.poly1d(Fit.polyfit_clip(xposs, widths, i))

pl.figure(2)
pl.clf()
mn = np.min(fun(xposs))
pl.plot(xposs, yposs-mn, 'b.')
pl.plot(xposs, fun(xposs)-mn, 'r')
pl.ylabel("ypos")

res = fun(xposs) - yposs
sd = np.std(res)
ok = np.abs(res) < 2*sd


print "Fit Residuals Sigma: %5.3f P2V: %5.3f" % (np.std(res), res.max()-res.min())
print "Clipped Resid Sigma: %5.3f P2V: %5.3f" % (np.std(res[ok]), res[ok].max()-res[ok].min())


regout = ""

delt = 2048./50.
for i in range(50):
        x = delt * i + 1
        sx = x + 1
        ex = x + delt + 1
        sy = y + fun(sx) + 1
        ey = y + fun(ex) + 1
        #regout += "line(startx, starty, endx, endy)"
        regout  += "line(%f, %f, %f, %f)\n" % (sx, sy, ex, ey)

        sy += wfun(sx)
        ey += wfun(ex)
        #regout += "line(startx, starty, endx, endy)"
        regout  += "line(%f, %f, %f, %f) # color=blue\n" % (sx, sy, ex, ey)


f=open("/users/npk/desktop/ds9.reg","w")
f.write(regout)
f.close()
