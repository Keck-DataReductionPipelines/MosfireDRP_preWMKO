
'''
    MOSFIRE DRP Test

    23 Apr 2012 npk: Repackaged for experimentation
     4 Apr 2012 npk: Flat field test wrapper
'''



import MOSFIRE
import time
import os
from MOSFIRE import Flats, Options, IO
import numpy as np, pylab as pl
import pyfits

reload(Flats)
reload (IO)


if __name__ == "__main__":


    if True: # NGC 5053 - J
        name = "NGC5053"
        band = "J"
        nums = [36, 37, 38]
        roll = -5 # See ngc5053.py
    elif False: # test_marc Y
        name = "test_marc"
        band = "Y"
        nums = [31, 32, 33, 34, 35]
        roll = 0
    elif False: # test_marc J
        name = "test_marc"
        band = "J"
        nums = [26, 27, 28, 29]
        roll = 0
    elif False: # test_marc H
        name = "test_marc"
        band = "H"
        nums = [19, 20, 21, 22, 23, 24]
        roll = 0

    flatlist = []

    if roll == 0:
        for num in nums:
            flatlist.append("/users/npk/desktop/5apr/m120406_%4.4i.fits" % num)
    else:
        for num in nums:
            header, ff, bs, targs, ssl, msl, asl = IO.readmosfits(
                    "/users/npk/desktop/5apr/m120406_%4.4i.fits" % num)
            hdu = pyfits.PrimaryHDU(np.roll(ff, -5, 0), header)
            hdulist = pyfits.HDUList([hdu])

            for tbl in [targs, ssl, msl, asl]:
                hdu = pyfits.new_table(tbl)
                hdulist.append(hdu)

            fn = "/users/npk/desktop/5apr/roll_m120406_%4.4i.fits" % num
            os.remove(fn)
            hdulist.writeto(fn)

            flatlist.append(fn)

    
    print flatlist

    Flats.handle_flats(flatlist, name, band,  Options.flat)




