'''

DRP Experiment code. Take experiment drp3 and make a allslits file

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


reload(Fit)
reload(IO)



regout = ""
y = 2028
toc = 0
slits = np.load("slit-edges.npy")

width = 5000 
write_pos = 2500
out = np.zeros([write_pos,width])

fiducial_spec = np.array([])

pl.ion()
pl.figure(1)
pl.clf()


(minl, maxl) = (90000,0)

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

        dy = slitlet.shape[0]
        mid = dy/2.
        spec = slitlet[mid-3:mid+3, :]
        spec = np.median(spec, axis=0)

        continuum = median_filter(spec, 250)
        if fiducial_spec.shape == (0,):
                fiducial_spec = spec
                fiducial_spec /= continuum
                fiducial_spec /= fiducial_spec.max()

        soc = spec/continuum
        soc /= soc.max()
        as_per_pix = .18/1.4
        shift = -np.median(loc)/as_per_pix

        lags = np.arange(shift-300,shift+300,dtype=np.int)
        xc = Fit.xcor(soc[800:1300], fiducial_spec[800:1300], lags)
        #start_shell()
        shift = lags[np.argmax(xc)]



        center = width/2 - shift 
        print center, write_pos

        xs = np.arange(len(spec))
        pl.plot(xs - shift, spec/continuum)

        try:
                out[(write_pos - dy):write_pos, center-1024:center+1024] = slitlet
                if minl > center-1024: minl = center-1024
                if maxl < center+1024: maxl = center+1024
        except:
                print "Skipping %i %i %i" % (center, shift, write_pos)
                pass
        #regout += "line(startx, starty, endx, endy)"
        # text{%s}\n" 
        regout += "line(0, %i, 2200, %i)\n" % (write_pos-dy,write_pos-dy)
        cd = np.float(ssl[target]["Target_to_center_of_slit_distance"])
        cd /= .18
        regout += "line(0, %i, 2200, %i) # color=blue text={%s} \n" % (write_pos-dy/2 + cd, write_pos-dy/2 + cd, ssl[target]["Target_Name"])

        write_pos -= dy



try:
        f = open("demo.reg",'w')
        f.write(regout)
        f.close()
except:
        print "Could not write reg file"

out = out[:,minl:maxl]
