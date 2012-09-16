
import getpass
import os
import sys


try: import pyfits
except: print "Could not import pyfits, please reinstall stscipython"

try: import stsci
except: print "STScI python not properly installed"


print "Installing mospy to /usr/local/bin"

try:
    f = open("apps/mospy")
    lines = f.readlines()
    f.close()
except:
    print """Could not open or read the apps/mospy file. Make sure
that you run this file from the MOSFIRE directory"""

f = open("/usr/local/bin/mospy", "w")

for line in lines:
    line.replace("AAAA", os.getcwd())
    f.write(line)

f.close()
os.system('chmod a+x /usr/local/bin/mospy')


print """Default directories:

    Place raw data files here:
 /scr2/mosfire

    Reduced data will be placed
 /scr2/%s/mosfire_redux

    Place the bad pixel mask in
 /scr2/mosfire/badpixels/

""" % (getpass.getuser())

yorn = raw_input("Would you like to download the bad pixel mask [y/n]?")
if yorn == 'y' or yorn == 'Y':
    os.system("curl -O http://mosfire.googlecode.com/files/badpix_10sep2012.fits")

    os.system("mkdir -p /scr2/mosfire/badpixels")
    os.system("mv badpix_*.fits /scr2/mosfire/badpixels/")



