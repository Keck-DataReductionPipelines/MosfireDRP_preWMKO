
import getpass
import os
import sys


try: import pyfits
except: print "Could not import pyfits, please reinstall stscipython"

try: import stsci
except: print "STScI python not properly installed"


print "Installing mospy to /usr/local/bin"

try:
    f = open("apps/mospy_mac")
    lines = f.readlines()
    f.close()
except:
    print """Could not open or read the apps/mospy file. Make sure
that you run this file from the MOSFIRE directory"""

os.system("cp apps/mospy_mac /usr/local/bin/mospy")
os.system('chmod a+x /usr/local/bin/mospy')


print """Default directories:
    The DRP will be moved to:
 ~/mosdrp/DRP

    Place raw data files here:
 ~/mosdrp/data

    Reduced data will be placed
 ~/mosdrp/output

    The bad pixel mask in
 ~/mosdrp/badpixels/

    If mospy issues a command not found, place
      /usr/local/bin in your PATH
""" 

os.system("mkdir -p $HOME/mosdrp/DRP")
os.system("mkdir -p $HOME/mosdrp/data")
os.system("mkdir -p $HOME/mosdrp/output")
os.system("mkdir -p $HOME/mosdrp/badpixels")
os.system("cp -r ~/mosfire/* $HOME/mosdrp/DRP/")
os.system("cp -r ~/mosfire/.hg $HOME/mosdrp/DRP/")

yorn = raw_input("Would you like to download the bad pixel mask [y/n]?")
bpm = 'badpix_10sep2012.fits'
if yorn == 'y' or yorn == 'Y':

    os.system("curl -O http://mosfire.googlecode.com/files/%s" % bpm)
    os.system("mv badpix_*.fits $HOME/mosdrp/badpixels/")


try:
    f = open(os.path.expanduser("~/mosdrp/DRP/MOSFIRE/Options.py"))
    lines = f.readlines()
    f.close()
except:
    print "Could not open and read Options.py in DRP directory"
    sys.exit()

f = open(os.path.expanduser("~/mosdrp/DRP/MOSFIRE/Options.py"), "w")

path = os.path.expanduser('~')
for line in lines:
    sp = line.split('=')
    if sp[0].rstrip() == 'indir': 
        sp[1] = " '%s/mosdrp/data'" % (path)
    if sp[0].rstrip() == 'outdir': 
        sp[1] = " '%s/mosdrp/output'" % (path)
    if sp[0].rstrip() == 'path_bpm': 
        sp[1] = " '%s/mosdrp/badpixels/%s'" % (path, bpm)

    outstr = '='.join(sp) + "\n"
    f.write(outstr)

f.close()


print("I recommend that you remove the ~/mosfire directory now.")
