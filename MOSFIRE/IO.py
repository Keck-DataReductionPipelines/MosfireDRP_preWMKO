'''
MOSFIRE Input/Output Utility Code
Written March 2, 2011 by npk

Provides tools to read fits files and parse their headers.
'''

import numpy as np
import pyfits as pf
import CSU

import unittest



def readfits(path):
        '''Read a fits file from path and return a tuple of (header, data).'''
        hdulist = pf.open(path)
        header = hdulist[0].header
        data = hdulist[0].data

        return (header, data)

def parse_header_for_bars(header):
        '''Parse {header} and convert to an array of CSU bar positions in mm. If the positon is negative it means the barstat is not OK'''

        poss = []
        posfmt = "B%2.2iPOS"
        statfmt = "B%2.2iSTAT"
        for i in range(1,CSU.numbars+1):
                p = posfmt % i
                s = statfmt % i
                pos = np.float32(header[p])
                if header[s] != 'OK':
                        pos *= -1
                poss.append(pos)

        if len(poss) != CSU.numbars:
                raise CSU.MismatchError("Found %i bars instead of %i" % (lens(poss), CSU.numbars))
                

        return np.array(poss)


class TestIOFunctions(unittest.TestCase):

        def setUp(self):
                pass

        def test_readfits(self):
                self.assertTrue(True)

if __name__ == '__main__':
        unittest.main()
