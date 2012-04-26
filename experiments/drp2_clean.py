'''

===================
MOSFIRE Flat Fields 
===================




npk April 14th 2011

'''
import MOSFIRE
import time
from MOSFIRE import Fit, IO
import numpy as np, pylab as pl

reload(Fit)
reload(IO)

if __name__ == "__main__":
        (header, data1, targs, ssl, msl, asl) = IO.readfits_all("/users/npk/desktop/c9/m110326_3242.fits")
        data = data1

        ssl = ssl[ssl["Slit_Number"] != ' ']
        numslits = np.round(np.array(ssl["Slit_length"], dtype=np.float) / 7.02)

        for i in range(len(ssl)):
                print ssl[i]["Target_Name"], numslits[i]





        Outputs:
        xposs []: Array of x positions along the slit edge [pix]
        yposs []: The fitted y positions of the "top" edge of the slit [pix]
        widths []: The fitted delta from the top edge of the bottom [pix]
        '''

        x = 936
        for i in range(-49, 50):
                delt = 1900/100.
                xp = x+delt*i
                v = data[y-roi_width:y+roi_width, xp-2:xp+2]
                v = np.median(v, axis=1)

                if np.median(v) < 200:
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

def fit_edge_poly(xposs, yposs, widths, order):
        '''
        fit_edge_poly fits a polynomial to the measured slit edges.
        This polynomial is used to extract spectra.

        input-
        xposs, yposs [N]: The x and y positions of the slit edge [pix]
        widths [N]: the offset from end of one slit to beginning of another [pix]
        order: the polynomial order
        '''

        fun = np.poly1d(Fit.polyfit_clip(xposs, yposs, order))
        wfun = np.poly1d(Fit.polyfit_clip(xposs, widths, order))
        res = fun(xposs) - yposs
        sd = np.std(res)
        ok = np.abs(res) < 2*sd

        return (fun, wfun, res, sd, ok)


def data_quality_plots(results):



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

        yposs = []
        widths = []
        xposs = []

        # TODO: Make these options
        roi_width = 15
        order = 5

        (xposs, yposs, widths) = find_edges(data, y, roi_width)
        (fun, wfun, res, sd, ok) = fit_edge_poly(xposs, yposs, widths, order)

        bottom = fun.c.copy() 
        top = wfun.c.copy() + fun.c.copy()

        bottom[-1] += y + 1
        top[-1] += y + 1

        result = {"y": y}
        result["targ"] = ssl[target]["Target_Name"]
        result["xposs"] = xposs
        result["yposs"] = yposs
        result["widths"] = widths
        result["fun"] = fun
        result["wfun"] = wfun
        result["sd"] = sd
        result["ok"] = ok
        result["top"] = np.poly1d(top)
        result["bottom"] = np.poly1d(bottom)
        result["poly_order"] = order

        results.append(result)

        print fun
        print "Clipped Resid Sigma: %5.3f P2V: %5.3f" % (np.std(res[ok]), res[ok].max()-res[ok].min())

        tic = time.clock()

        print " .... %4.2f s elapsed." % (tic - tock)
        print

np.save("slit-results", results)
