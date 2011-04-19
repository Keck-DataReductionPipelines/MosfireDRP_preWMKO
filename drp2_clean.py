'''

DRP Experiment code. Measure slit edges on a set of bars based on header

npk March 23rd 2011
npk April 14th 2011 -- Convert to a more readable code

Conclusions:
        5-order polynomial is sufficient

'''
import MOSFIRE
import time
from MOSFIRE import Fit, IO
import numpy as np, pylab as pl
pl.ion()
pl.figure(1)
pl.clf()

(header, data1, targs, ssl, msl, asl) = IO.readfits_all("/users/npk/desktop/c9/m110326_3242.fits")
data = data1

ssl = ssl[ssl["Slit_Number"] != ' ']
numslits = np.round(np.array(ssl["Slit_length"], dtype=np.float) / 7.02)

for i in range(len(ssl)):
        print ssl[i]["Target_Name"], numslits[i]


reload(Fit)
reload(IO)

roi_width = 15 


def fit_edges(data, y)
        '''
        fit_edges finds slit edges near a y position

        data is a MOSFIRE data frame 2k x 2k, units are irrelevant
        y is the rough starting position near the center of the frame

        returns three vectors
        xposs -- The x location of the fit
        yposs -- The y location of the fit
        widths -- The distance between two slit edges
        '''
        global roi_width
        yposs = []
        widths = []
        xposs = []
        x = 936
        for i in range(-49, 50):
                delt = 1900/100.
                xp = x+delt*i
                v = data[y-roi_width:y+roi_width, xp-2:xp+2]
                v = np.median(v, axis=1)

                if np.median(v) < 300:
                        continue

                ff = Fit.do_fit(v, residual_fun=Fit.residual_disjoint_pair)
                if (0 < ff[4] < 4):
                        xposs.append(x+delt*i)
                        xs = np.arange(len(v))

                        yposs.append(ff[0][1] - roi_width)
                        widths.append(ff[0][5])
                else:
                        print "Skipping: %i" % (x+delt*i)

        return map(np.array, (xposs, yposs, widths))



regout = ""
y = 2028
toc = 0
slits = []
top = [0., 2048.]
results = []
for target in range(len(ssl)-1):
        y -= 44.25 * numslits[target]
        print "-------------==========================------------------"
        print "Finding Slit Edges for %s starting at %4.0i" % (ssl[target]["Target_Name"], y)
        tock = time.clock()


        (xposs, yposs, widths) = fit_edges(data, y)

        result = {"y": y}
        result["xposs"] = xposs
        result["yposs"] = yposs
        result["widths"] = widths

        fun = np.poly1d(Fit.polyfit_clip(xposs, yposs, 5))
        wfun = np.poly1d(Fit.polyfit_clip(xposs, widths, 5))
        res = fun(xposs) - yposs
        sd = np.std(res)
        ok = np.abs(res) < 2*sd

        result["fun"] = fun
        result["wfun"] = wfun
        result["widths"] = widths
        result["sd"] = sd
        result["ok"] = ok

        print fun
        print "Clipped Resid Sigma: %5.3f P2V: %5.3f" % (np.std(res[ok]), res[ok].max()-res[ok].min())


        delt = 2048./50.
        bottom = fun.c.copy() 
        bottom[-1] += y + 1

        result["targ"] = ssl[target]["Target_Name"]

        slits.append({"targ": ssl[target]["Target_Name"], "top": np.poly1d(top), "bottom": np.poly1d(bottom)})

        top = wfun.c.copy() + fun.c.copy()
        top[-1] += y + 1

        result["top"] = np.poly1d(top)
        result["bottom"] = np.poly1d(bottom)

        results.append(result)


        tic = time.clock()

        print " .... %5.1f s elapsed." % (tic - tock)
        print

pl.ion()


# The file slit-edges goes to drp3.py
np.save("slit-edges", slits)
f=open("/users/npk/desktop/ds9.reg","w")
f.write(regout)
f.close()
