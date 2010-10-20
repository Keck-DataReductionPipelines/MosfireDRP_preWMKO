# -*- coding: utf-8 -*-
from numpy import *
from pylab import *
from scipy import *
import pyfits
from pyraf import iraf
import pydrizzle
from math import * 
import nmpfit_sy
import numpy as np
import pylab as py
import random

#This code loads a fits image (identified on line 21) of the slit mask and computes the slit locations.
#The fit is made using a gaussian. Header information of expected bar locations be is loaded for comparison.
#Output includes a tab seperated text file and a png figure of the slit locations.
#Version 2.0
#
# 2010-10-19  GNM

f = pyfits.open('m100917_0316.ref.fits')
prihdr = f[0].header
scidata = f[0].data

def fitGaussianMP(p0=None,data=None,quiet=0):
    """Fits Gaussian using mpfit.

    Inputs should be given as p0=[A, mu, sigma] and
    data=[x, occ, err], where  x is the deviation in
    units of sigma and occ is the occupation of the
    histogram bin.    Returns object of class mpfit.
    From S. Yelda, 10.07.2010.
    """

#    print 'Initial Guess:'
#    print '   A     = %6.2f' % p0[0]
#    print '   mu    = %5.3f' % p0[1]
#    print '   sigma = %5.3f' % p0[2]

#    # Remove data with zero error (residuals go as 1/err)
#    x   = (data[0])[np.nonzero(data[2])]
#    occ = (data[1])[np.nonzero(data[2])]
#    err = (data[2])[np.nonzero(data[2])]

    # Set data to be passed to fit
    functargs = {'x':L,'occ':M,'err':np.sqrt(M)}

    # Set initial values and limits (no limits on parameters)
    pinfo = [{'value':0,'fixed':0,'limited':[0,0],
	      'limits':[0,0]}]*len(p0)
    for ii in range(len(p0)):
        pinfo[ii]['value'] = p0[ii]

    # Use mpfit fitting algorithm to fit parameters
    m = nmpfit_sy.mpfit(fitfuncGauss, p0, functkw=functargs, parinfo=pinfo,
		    quiet=quiet)
    if (m.status <= 0):
        print 'Error message = ', m.errmsg
    
    return m

def fitfuncGauss(p, fjac=None, x=None, occ=None, err=None):
    """Find residuals of Gauss fit.

    For Gaussian, p should be list of form [A, mu, sigma],
    while data should be list of form [x, occ, err], where 
    x is the deviation in units of sigma and occ is the
    occupation of the histogram bin.
    From S. Yelda, 10.07.2010.
    """

    num = len(L)
    model = np.zeros(num, dtype=float)
    devs = np.zeros(num, dtype=float)

    # Set parameters
    # Include normalization constant A as the data
    # is assumed not to be normalized
    A = p[0]
    mu = p[1]
    sigma = p[2]

    model = modelGaussian(L, A, mu, sigma)
    residuals = (occ - model)/err
    status = 0

    return [status, residuals]


def modelGaussian(L, A, mu, sigma): #From S. Yelda, 10.07.2010.
    model = rowmedian + (A/(sigma*np.sqrt(2*np.pi)))*np.exp(-.5*((L-mu)/sigma)**2)
    return model

def getBars():
    bars = zeros(92)
    for i in range(1,93):
        bars[i-1] = prihdr["B%2.2iPOS" % i]
    return bars
    
def slitBars():
    bars = getBars()
    slitbar = zeros(46)
    for i in range(0,46):
        slitbar[i] = bars[i*2] + ((bars[i*2+1]-bars[i*2])/2)
    return slitbar

def write_to_file():
    out = open('slitcheck.txt','w')
    hdr = '%4s  %4s  %4s  %10s  %4s  %4s  %11s  %4s  %4s  %6s\n'
    fmt = '%2.2i  %4.2i  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f  %7.2f\n'
    out.write(hdr % ('Spec','Row','Bar','Bar Center','Bar', 'Edge', 'Slit Center', 'Edge', 'FWHM','Offset'))
    
    for n in range(0,46):
        out.write(fmt % (n+1, y[n], bars_pix[n*2+1], BSC[n], bars_pix[n*2], Sedge2[n], SC[n], Sedge1[n], FWHM[n], SE[n]))

    out.close()

def make_image_plot(): # make figure and save to file
    prop = matplotlib.font_manager.FontProperties(size=10) 
    py.figure(5,figsize=(8,8))
    py.clf()
    py.plot(SC,y,'r,',lw=0.5)
    py.plot(Sedge1,y,'r|',lw=1.5)
    py.plot(Sedge2,y,'r|',lw=1.5)
    py.plot(BSC,y,'b,',lw=0.5)
    py.plot(bars_pix,ybars,'b|',lw=1.5)
    py.plot(badflag_x,badflag_y,'kx',lw=3)
    py.errorbar(SC,y,ones(46),SEP,',')
    py.axis([0,2048,0,2048])  
    bw = imread('Hx.png')
    imshow(bw,aspect='equal',origin='upper', extent=(0,2048,-190,2048))
    py.title('Imaged Bar Locations and Header Bar Positions')
    py.annotate('Header', xy=(15, 15), xytext=(1500, -350), color='b')
    py.annotate('Image', xy=(15, 15), xytext=(1800, -356), color='r')
    py.annotate('Slit 1', xy=(15, 15), xytext=(2060, 2022))
    py.annotate('Bar 1', xy=(15, 15), xytext=(2060, 1978))
    py.annotate('Slit 46', xy=(15, 15), xytext=(2060, 22))
    py.annotate('Bar 91', xy=(15, 15), xytext=(2060, -22))
    py.xlabel('Pixels')
    py.ylabel('Pixels')
    py.savefig('slitcheck_image.png')
    py.close(5)
    
offset = 44 # The height of the bottom spectrum from the bottom of the image.
height = 44 # 2048 - edges / 46
X_center = 1033.0 # CSU x center on detector from Chuck's analysis
Y_center = 1033.0 # CSU y center on detector from Chuck's analysis
magnif = (0.018/0.18025)/0.72515
pixmm = 1/0.018053 # Should be 0.018, Removes linear term in offset, as distance from center increases due to poor row definition.
rotation = 0.49
slitcenters_y = [] # Not used. Will be used when flat images are used to define the slit locations
loc_median = [ ] 
slitnum = 1 
slitnumber = [ ]
SE = [ ]#Offset between expected and actual slit locations.
BSC = [ ] #Bar slit center
Sedge1 = [ ] #Slit edge from image
Sedge2 = [ ] #Slit edge form image
SC = [ ] #Slit center in image
y = [ ]#rows where slits are defined
ybars = [ ]# same as y with double entries since 2 bars per row.
badflag_x = [ ]# saves info for slits that cause errors in the code
badflag_y = [ ]# saves info for slits that cause errors in the code
im_median = median(scidata)
FWHM = [ ]
bars_pix = zeros(92)
SEP = [ ]
eplot = 200

print ' ------------------------- '

for spec in range(len(scidata[0])-10, offset,-height):
    bars = getBars()
    if spec > 1024:
        bars_pix [slitnum*2-2] = X_center + ((137.4-bars[slitnum*2-2])*pixmm*magnif) - (spec-1024)*tan(radians(rotation))
        bars_pix [slitnum*2-1] = X_center + ((137.4-bars[slitnum*2-1])*pixmm*magnif) - (spec-1024)*tan(radians(rotation))
    else:
        bars_pix [slitnum*2-2] = X_center + ((137.4-bars[slitnum*2-2])*pixmm*magnif) - (spec-1024)*tan(radians(rotation))
        bars_pix [slitnum*2-1] = X_center + ((137.4-bars[slitnum*2-1])*pixmm*magnif) - (spec-1024)*tan(radians(rotation))
    slitbar_pix = bars_pix[slitnum*2-1] + ((bars_pix[slitnum*2-2] - bars_pix[slitnum*2-1])/2)
    rowmedian = median(scidata[spec,:])
    edge = 2 * rowmedian  #Should be replaced with formal definition using FWHM or something.
    
    L = [ ]
    i = 0
    for val in scidata[spec,:]: # This loop makes a first pass at finding slits, but also identifies hot pixels.
        i = i+1
        if val >= 8*rowmedian:
            L.append(i)
    loc = median(L)
    width = len(L)

    if rowmedian > 8*im_median: # This loop removes mostly open slits from the analysis due to difficulties in gaussian fitting.
        print 'Slit', slitnum, 'is too wide and has been flagged as bad. Row:', spec
        SE.append(1024)
        SEP.append(0)
        BSC.append(1024)
        Sedge1.append(2038) #Slit edge from image
        Sedge2.append(10) #Slit edge form image
        SC.append(1024) #Slit center from image
        y.append(spec)
        ybars.append(spec)
        ybars.append(spec)
        slitnumber.append(slitnum)
        badflag_x.append(1024)
        badflag_y.append(spec)
        FWHM.append(500)
        slitnum = slitnum + 1
        continue
    
    else: # This loop is the primary code.
        L = [ ]
        M = [ ]
        i = 0
        for val in scidata[spec,:]: # This loop tkae the first pass from above and clips it down to just the actual slit.
            i = i+1
            if val >= edge:
                if i in range(loc-width,loc+width):
                    L.append(i)
                    M.append(val)
        loc_mean = mean(L) 
        loc_width = std(L)
        
        if L == [ ]: # This loop identifies closed slits and removes them from analysis.
            print 'Slit', slitnum, 'was not found and has been flagged as bad. Row:', spec
            SE.append(1024)
            SEP.append(0)
            BSC.append(1024)
            Sedge1.append(2038) #Slit edge from image
            Sedge2.append(10) #Slit edge form image
            SC.append(1024) #Slit center from image
            y.append(spec)
            ybars.append(spec)
            ybars.append(spec)
            slitnumber.append(slitnum)
            badflag_x.append(1024)
            badflag_y.append(spec)
            FWHM.append(0.01)
            slitnum = slitnum + 1
            continue


        print "For spec", slitnum , "at row", spec, ":"
        slitbar = slitBars()
        p0 = [max(M), mean(L), std(L)]
        rfit = fitGaussianMP(p0,[L, M, np.sqrt(M)],1) # Gaussian fit to slit array.
        m = rfit
        p = m.params          # Best-fit parameters
        perr = m.perror       # Error in parameter fits from covariance matrix
        m.dof = len(L)-len(p) # Number of degrees of freedom
        Rchi2 = m.fnorm/m.dof # Reduced Chi^2 statistic
        amp = rowmedian + (p[0]/(p[2]*np.sqrt(2*np.pi)))*np.exp(-.5*((p[1]-p[1])/p[2])**2)

        slit_error = p[1] - slitbar_pix # determines the slit offsets for V_2.0

        print 'Row Median (Counts)= %5.2f' % rowmedian
        print 'Slit Center(pixels)= %5.3f +/- %5.3f' % (p[1],perr[1])
        print 'FWHM               = %5.3f +/- %5.3f' % (p[2],perr[2])
        print 'Gaussian Amplitude = %5.3f +/- %5.3f' % (amp,sqrt(amp))
#       print 'chi^2  = %5.2f' % m.fnorm
#       print 'dof    = %2d' % m.dof
        print 'Rchi^2             = %5.2f' % Rchi2
        print 'Slit Offset        = %5.3f' % slit_error
#        print 'Bar Postion        = %5.3f' % bars_pix[slitnum*2-1]
#        print 'Bar Slit Center    = %5.3f' % slitbar_pix
#        print 'Bar Position       = %5.2f' % bars_pix[slitnum*2-2]
        print ' ------------------------- '
    
        width = len(L)
        fwhm = p[2]
        edge1 = p[1] + fwhm
        edge2 = p[1] - fwhm
        
        # Saves the computed values into arrays for figure and text files.
        SE.append(slit_error) #slit offsets
        SEP.append(slit_error*eplot)
        BSC.append(slitbar_pix)#Slit center from bars
        Sedge1.append(edge1) #Slit edge from image
        Sedge2.append(edge2) #Slit edge form image
        SC.append(p[1]) #Slit center from image
        y.append(spec)#Row where computation made
        ybars.append(spec)
        ybars.append(spec)
        slitnumber.append(slitnum)
        FWHM.append(fwhm)
        slitnum = slitnum + 1
    

make_image_plot()
write_to_file()

print 'Median Slit Offset = %5.3f' % median(SE)
