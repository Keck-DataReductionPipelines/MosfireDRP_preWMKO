'''

Written March 3rd 2011 by npk

This file is part of an experiment to measure CSU positions from test data.
It is used in conjunction with CSU_plot3.py and CSU_plot_confirm.py
'''

import sys, datetime, getpass, os
import numpy as np, pyfits as pf
import pylab as pl
import scipy.ndimage.filters, scipy.io
from pyraf import iraf
from MOSFIRE import CSU, Detector, IO, Fit

from IPython.Shell import IPShellEmbed

start_shell = IPShellEmbed()

NaN = np.nan


def save_region(bs, fname):
        s = bs.to_ds9_region()
        try:
                f = open(fname, "w")
                f.writelines(s)
                f.close()
        except:
                return Exception("Could not write barset region file to: %s" % fname)

def fit_line_with_sigclip(xs, data, i=0):


        ps = Fit.do_fit_edge(xs, data)
        pf = lambda x: Fit.fit_bar_edge(ps, x)


        residual = np.abs(pf(xs) - data)
        sd = np.std(residual)

        ok = np.where(residual < 2.5 * sd)[0]
        if len(ok) == 0:
                return [lambda x: NaN, []]

        ps = Fit.do_fit_edge(xs[ok], data[ok])
        pf = lambda x: Fit.fit_bar_edge(ps, x)

        return [pf, ok]



if False:
        def fit_line_with_sigclip(xs, data, i = 0):

                ps = np.polyfit(xs, data, 1)
                pf = np.poly1d(ps)


                residual = np.abs(pf(xs) - data)
                sd = np.std(residual)
                
                ok = np.where(residual < 2.5*sd)[0]

                ps = np.polyfit(xs[ok], data[ok], 1)
                pf = np.poly1d(ps)
                return [pf, ok]

                if len(ok) == len(residual):
                        return [pf, ok]
                elif i > 2:
                        return [pf, ok]
                else:
                        return fit_line_with_sigclip(xs[ok], data[ok], i+1)

def median_tails(v):
        a = np.median(v[0:2])
        b = np.median(v[-3:-1])

        t = v - np.float(a+b)/2.

        return t

def make_slice(pos, offset, w, h):
        '''Returns [ [xslice, ylsice], [x0, x1, y0, x1] ] where
        xslice is used as Array[x0:x1]
        yslice is used as Array[y0:y1]'''
        x0 = round(pos[0]- w + offset)
        if x0 < 0: x0 = 0
        x1 = round(pos[0] + w + offset)
        if x1 > 2047: x1 = 2047
        if x0 > x1: x0 = x1
        xs = slice(x0, x1)

        y0 = round(pos[1]-h)
        if y0 < 0: y0 = 0
        y1 = round(pos[1]+h)
        if y1 > 2047: y1 = 2047
        if y0 > y1: y0 = y1
        ys = slice(y0,y1)


        return [[xs,ys],[x0,x1,y0,y1],round(pos[1])]

#k = np.array([[1,1,1],[1,1,1],[1,1,1]], dtype=np.float32)
#k /= k.sum()
#data = scipy.ndimage.filters.convolve(data,k)
#data = scipy.ndimage.filters.median_filter(data,size=(11,11))
#data = scipy.ndimage.filters.gaussian_filter(data,3)

deg = np.pi/180.
np.set_printoptions(precision=2)
np.set_printoptions(suppress=True)

reload(CSU)
reload(IO)
reload(Fit)
reload(Detector)

pl.ion()


def sigclip(data, low=4, high=4):
        c = data.ravel()
        delta = 1
        while delta:
                s = c.std()
                m = c.mean()
                size = c.size

                c = c[(c > (m - s*low)) & (c < (m + s*high))]

                delta = size-c.size

        return c.mean()


def plot_stamps(stamp, ps):
        global plot_arr
        pl.figure(1)
        pl.subplot(*plot_arr)
        pl.imshow(stamp)

def plot_linefits(y, ff):
        global plot_arr
        pl.figure(3)
        pl.subplot(*plot_arr)
        ps = Fit.do_fit(y, Fit.residual_single)[0]

        xx = np.arange(len(y))
        xxx = np.arange(0,len(y)-1,.1)
        fun = Fit.fit_single(ps,xxx)
        pl.plot(xx,y, '*')
        pl.plot(xxx, fun)
        x0 = ff(0)
        x1 = ff(0) + ps[4]
        pl.plot([x0, x0], [0, Fit.fit_single(ps, x0)])
        pl.plot([x1, x1], [0, Fit.fit_single(ps, x1)])

        pl.figure(4)
        pl.subplot(*plot_arr)
        fun = Fit.fit_single(ps,xx)
        r = (y - fun)/np.sqrt(np.abs(fun))
        pl.plot(xx, r, '*')
        pl.ylim([-5,5])
        
def plot_ridgeline(fits, ok, xs, ff, bar):
        global plot_arr

        pl.figure(5)
        if plot_arr[2] == 10: pl.title("Bar Edge Ridgeline: Delta pixels v Pixel on slit")
        pl.subplot(*plot_arr)

        pl.plot(xs, fits[:,1] - ff(xs))
        pl.plot(xs, fits[:,1] - ff(xs), '*-')
        pl.plot(xs[ok], fits[ok,1] - ff(xs[ok]), 'or-')
        pl.ylim([-.1,.1])
        pl.title("%2i" % bar)

        pl.figure(6)
        pl.subplot(*plot_arr)

        pl.plot(xs, fits[:,1])
        pl.plot(xs, fits[:,1], '*-')
        pl.plot(xs[ok], fits[ok,1], 'or-')
        pl.title("%2i" % bar)



def plot_widths(fits, ok, xs):
        global plot_arr

        # fits[4] is slit width
        pl.figure(7)
        pl.subplot(*plot_arr)

        pl.plot(xs, fits[:,4], '*-')
        pl.plot(xs[ok], fits[ok, 4], 'or-')

def is_in_bounds(extent):
        if extent[1] > 2047: return False
        if extent[0] < 0: return False
        if extent[2] > 2047: return False
        if extent[0] < 0: return False
        if extent[0] == extent[1]: 
                return False
        if extent[2] == extent[3]: 
                return False

        for i in [0,1]:
                for j in [2,3]:
                        if not CSU.in_field(extent[i],extent[j]): return False
        return True

def go(fname):
        global plot_arr

        print fname
        
        (header, data) = IO.readfits(fname)
        bs = CSU.Barset()
        bs.set_header(header)

        cntfit = 1
        for i in range(1,8):
                continue
                pl.figure(i)
                pl.clf()

        first_time = True

        bars = []
        means = []
        sds = []
        deltas = []
        deltas_mm = []
        poss = []
        poss_mm = []
        poss_y = []
        request = []
        qs = []
        bars = range(1, 93)
        for bar in bars:
                pos = bs.get_bar_pix(bar)
                if bar % 8 == 0:
                        print "%2.0i: (%7.2f, %7.2f)" % (bar, pos[0], pos[1])

                width = 19 
                [[xslice, yslice],extent,ystart] = make_slice(pos,0,width,30)


                if not is_in_bounds(extent):
                        fits = [0,0,0,0,0]
                        [ff,ok] = [np.poly1d(0,0), []]
                        means.append(fits)
                        sds.append(fits)
                        drop_this = True
                else:
                        drop_this = False
                        fits = []
                        ys = np.arange(-10,10, dtype=np.float32)
                        for i in ys:
                                tofit = data[ystart-i, xslice]
                                y = median_tails(tofit)

                                ps = Fit.do_fit(y, Fit.residual_single)
                                fits.append(ps[0])
                        
                        fits = np.array(fits)
                        fits[:,1] += 1

                        # fit to the ridgeline
                        [ff, ok] = fit_line_with_sigclip(ys, fits[:,1])
                        m = [np.mean(fits[:,i]) for i in range(5)]
                        s = [np.std(fits[:,i]) for i in range(5)]
                        means.append(m)
                        sds.append(s)

                        #if bar == 44: start_shell()
                        #if bar == 43: start_shell()

                        # End of measurement logic, PLOTS

                        plot_arr = [12, 8, cntfit]
                        cntfit += 1
                        if False:
                                plot_stamps(data[yslice, xslice], ps[0])
                                y0 = median_tails(data[ystart, xslice])
                                plot_linefits(y0, ff)
                                plot_ridgeline(fits, ok, ys, ff, bar)
                                plot_widths(fits, ok, ys)

                        # End of plots, BOOKEEPING


                slit_center_offset = pos[1] - ystart
                fc = ff(slit_center_offset)
                slit_center_pos = np.float(extent[0] + fc )

                if drop_this: 
                        poss.append(NaN)
                        poss_y.append(NaN)
                        poss_mm.append(NaN)
                else: 
                        poss.append(slit_center_pos)
                        poss_y.append(ystart)
                        poss_mm.append(CSU.csu_pix_to_mm_poly(slit_center_pos, ystart)[0])

                delta = np.float(slit_center_pos - pos[0])
                if drop_this: 
                        deltas.append(NaN)
                        deltas_mm.append(NaN)
                else: 
                        deltas.append(delta)
                        b = CSU.csu_pix_to_mm_poly(slit_center_pos + delta, ystart)[0]
                        deltas_mm.append(b - poss_mm[-1])

                q = np.float(np.degrees(np.tan(ff(1)-ff(0))))
                if drop_this: qs.append(NaN)
                qs.append(q)

                if False: #is_in_bounds(extent):
                        pl.figure(2)
                        pl.text(pos[0], pos[1], 'b%2.0i: w=%3.2f p=%5.2f q=%3.2f d=%1.3f' % (bar, np.mean(fits[:,4]), extent[0]+ff(0), q, delta), fontsize=11, family='monospace', horizontalalignment='center')
                        


        means = np.array(means)
        f = lambda x: np.array(x).ravel()
        sds = f(sds)
        deltas = f(deltas)
        poss = f(poss)
        poss_y = f(poss_y)
        poss_mm = f(poss_mm)
        deltas_mm = f(deltas_mm)
        qs  = f(qs)
        bars = f(bars)
        pl.draw()


        if False:
                [pf, ok] = fit_line_with_sigclip(bars, deltas)
                print "** Residual: %4.3f [pix]" % np.std(deltas[ok]-pf(bars[ok]))
                evens = range(0,92,2)
                odds = range(1,92,2)

                pl.plot(bars[evens], deltas[evens],'r*-')
                pl.plot(bars[odds], deltas[odds],'b*-')
                pl.title("Red: even, blue: odd")
                pl.xlabel("Bar #")
                pl.ylabel("Delta pixel (Measured - Predicted)")
                pl.plot(bars, pf(bars))


        fout = "/users/npk/dropbox/mosfire/cooldown 9/csu_meas/" + fname.split("/")[-1] + ".sav"
        print "saving"
        tosav = {"bars": bars, "request": bs.pos, "deltas_mm": deltas_mm, "poss": poss, "poss_mm": poss_mm, "deltas": deltas, "means": means, "qs": qs}
        scipy.io.savemat(fout, tosav)
        save_region(bs, "/users/npk/dropbox/mosfire/cooldown 9/csu_meas/" + fname.split("/")[-1] + ".reg")
        print "saved"

        regout = "/users/npk/dropbox/mosfire/cooldown 9/csu_meas/" + fname.split("/")[-1] + ".meas.reg"
        pairs = np.array([poss,poss_y]).transpose()
        s = CSU.to_ds9_region(pairs, dash=0, color="blue", label=False)
        try:
                f = open(regout, "w")
                f.writelines(s)
                f.close()
        except:
                print "Couldn't write: %s" % regout


        return [tosav, bs]


                
fname = "/users/npk/desktop/c8/m101029_0233.ref.fits"
fname = "/users/npk/desktop/c9/m110311_0130.fits"
fname = "/users/npk/desktop/c9/m110312_0156.fits"
fname = "/users/npk/desktop/c9/m110312_%4.4i.fits"
fname = "/users/npk/desktop/c9/m110312_%4.4i.fits"

path  = "/users/npk/desktop/c9/"
def generate_fname(num):
        global path
        files = os.listdir(path)
        for fn in files:
                if fn.find("12_%4.4i" % num) > 0:
                        return fn
                if fn.find("13_%4.4i" % num) > 0:
                        return fn
                if fn.find("16_%4.4i" % num) > 0:
                        return fn
                if fn.find("17_%4.4i" % num) > 0:
                        return fn
                if fn.find("18_%4.4i" % num) > 0:
                        return fn


# 426 is the result of some bug where headers are not populated.
exclude = [426,392,393,394,395] 
#for i in range(156, 192):

if True:
        for i in range(212, 555):
                if i in exclude: continue
                fname = generate_fname(i)

                if fname is None:
                        print "Ignoring file #%i" % i
                        continue

                print fname
                go(path + fname)
else:

        #odd = generate_fname(300)
        #even = generate_fname(403)

        #for i in range(1550,1568):
        for i in range(1781, 1793):
                print i
                [t,b] = go(path + generate_fname(i))
                break

        #[to, bo] = go(path+odd)
        #[te, be] = go(path+even)




