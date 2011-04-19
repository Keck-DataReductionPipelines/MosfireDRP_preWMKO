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
from scipy.signal.signaltools import correlate
from scipy.ndimage.filters import median_filter

from IPython.Shell import IPShellEmbed

reload(Fit)
reload(IO)

start_shell = IPShellEmbed()

global data, header, targs, ssl, msl, asl
try:
        data
except:
        print "reading data"
        (header, data1, targs, ssl, msl, asl) = IO.readfits_all("/users/npk/desktop/c9/m110323_2707.fits")
        data = data1

ssl = ssl[ssl["Slit_Number"] != ' ']
numslits = np.round(np.array(ssl["Slit_length"], dtype=np.float) / 7.02)

for i in range(len(ssl)):
        print ssl[i]["Target_Name"], numslits[i]



regout = ""
for target in range(len(ssl)-1):
        y -= 44.25 * numslits[target]
        print "-------------==========================------------------"
        print "Extracting slit for %s starting at %4.0i" % (ssl[target]["Target_Name"], y)
        tock = time.clock()

        yposs = []
        widths = []
        xposs = []


        assert(slits[target]["targ"] == ssl[target]["Target_Name"])

        roi = slits[target]["targ"] == msl["Target_in_Slit"]
        loc = np.array(msl[roi]["Position_of_Slit"], dtype=np.float)
        print "Located : ", np.median(loc)

        x = np.arange(2048)

        top = slits[target]["top"](x)
        bottom = slits[target]["bottom"](x)
        mn = np.floor(bottom.min())
        mx = np.ceil(top.max())

        bbox = [slice(mn,mx), slice(0,2048)]
        slitlet = data[bbox]

        for i in range(2048):
                slitlet[0:np.ceil(bottom[i])-mn,i]=0
                slitlet[np.floor(top[i])-mn:slitlet.shape[0],i]=0
