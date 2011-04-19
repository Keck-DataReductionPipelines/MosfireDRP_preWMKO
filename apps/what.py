'''
MOSFIRE 'what' command:
        Spits out an informative summary of files
        in the current directory. Or files selected
        via a glob.

        i.e. what *0311.fits

npk - March 23 2011
'''

import MOSFIRE
import glob
import sys


files = []
if len(sys.argv) == 1:
        files = glob.iglob('*')
else:
        for i in range(1, len(sys.argv)):
                files.extend(glob.iglob(sys.argv[i]))


#print "filename              object  exptime                maskname lamp  filt   Turret"
for fname in files:

        try:
                (header, data) = MOSFIRE.IO.readfits(fname)
        except IOError:
                print "Skipping %s" % fname
                continue

        lamps = ""
        try:
                if header["pwstata7"] == 1:
                        lamps += header["pwloca7"][0:2]
                if header["pwstata8"] == 1:
                        lamps += header["pwloca8"][0:2]
        except KeyError:
                lamps = "???"
                
        header.update("lamps", lamps)
        try:
                print "%(datafile)12s %(object)25s %(truitime)6.1f %(maskname)25s %(lamps)3s %(filter)6s %(mgtname)7s" % (header)
        except:
                try:
                        print "%(datafile)12s %(object)25s %(truitime)6.1f %(lamps)3s %(filter)6s %(mgtname)7s" % (header)
                except:
                        print "%s Skipped" % fname