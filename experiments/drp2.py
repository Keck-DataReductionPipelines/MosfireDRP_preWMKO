'''

DRP Experiment code. Measure slit edges on a set of bars based on header

npk March 23rd 2011

Conclusions:
        5-order polynomial is sufficient

'''
import MOSFIRE
import time
from MOSFIRE import Fit, IO
import numpy as np, pylab as pl

(header, data1, targs, ssl, msl, asl) = IO.readfits_all("/users/npk/desktop/c9/m110326_3242.fits")
data = data1

ssl = ssl[ssl["Slit_Number"] != ' ']
numslits = np.round(np.array(ssl["Slit_length"], dtype=np.float) / 7.02)

for i in range(len(ssl)):
        print ssl[i]["Target_Name"], numslits[i]


reload(Fit)
reload(IO)

pl.figure(1)
pl.clf()

roi_width = 15 


regout = ""
y = 2028
toc = 0
slits = []
top = [0., 2048.]
for target in range(len(ssl)-1):
        y -= 44.25 * numslits[target]
        print "-------------==========================------------------"
        print "Finding Slit Edges for %s starting at %4.0i" % (ssl[target]["Target_Name"], y)
        tock = time.clock()

        yposs = []
        widths = []
        xposs = []

        x = 936
        for i in range(-49, 50):
                delt = 1900/100.
                xp = x+delt*i
                v = data[y-roi_width:y+roi_width, xp-2:xp+2]
                v = np.median(v, axis=1)

                ff = Fit.do_fit(v, residual_fun=Fit.residual_disjoint_pair)
                if (0 < ff[4] < 4):
                        xposs.append(x+delt*i)
                        xs = np.arange(len(v))

                        yposs.append(ff[0][1] - roi_width)
                        widths.append(ff[0][5])
                else:
                        print "Skipping: %i" % (x+delt*i)

        (xposs, yposs, widths) = map(np.array, (xposs, yposs, widths))

        fun = np.poly1d(Fit.polyfit_clip(xposs, yposs, 3))
        wfun = np.poly1d(Fit.polyfit_clip(xposs, widths, 3))
        res = fun(xposs) - yposs
        sd = np.std(res)
        ok = np.abs(res) < 2*sd

        print fun
        print "Clipped Resid Sigma: %5.3f P2V: %5.3f" % (np.std(res[ok]), res[ok].max()-res[ok].min())


        delt = 2048./50.
        bottom = fun.c.copy() 
        bottom[-1] += y + 1
        slits.append({"targ": ssl[target]["Target_Name"], "top": np.poly1d(top), "bottom": np.poly1d(bottom)})

        top = wfun.c.copy() + fun.c.copy()
        top[-1] += y + 1

        for i in range(50):
                x = delt * i + 1
                sx = x + 1
                ex = x + delt + 1
                sy = y + fun(sx) + 1
                ey = y + fun(ex) + 1

                #regout += "line(startx, starty, endx, endy)"
                if i == 25:
                        regout  += "line(%f, %f, %f, %f) # text={S%2.0i (%s)}\n" % (sx, sy, ex, ey, target+1, ssl[target]["Target_Name"])
                else:
                        regout  += "line(%f, %f, %f, %f)\n" % (sx, sy, ex, ey)

                sy += wfun(sx)
                ey += wfun(ex)
                #regout += "line(startx, starty, endx, endy)"
                regout  += "line(%f, %f, %f, %f) # color=blue\n" % (sx, sy, ex, ey)


        tic = time.clock()

        print " .... %5.1f s elapsed." % (tic - tock)
        print

# The file slit-edges goes to drp3.py
np.save("slit-edges", slits)
f=open("/users/npk/desktop/ds9.reg","w")
f.write(regout)
f.close()
