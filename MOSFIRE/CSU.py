

''' 
MOSFIRE CSU Utlity Code
        Includes physical parameters of CSU.

Created March 2, 2011 by npk

numslits, numbars give information about the number of slits and bars in the CSU.

Note, for consistency with the hardware bars and slits are indexed from 1.

tempscale is the thermal scaling factor to shrink room temperature linear dimensions to 120 K. The number 0.99646 is from R. Weber's spreadsheet "MOSFIRE Thermal Dimension Scaling Factors.xls", referenced from "Thermophysical properties of matter, Vol 12-13, Thermal Expansion".

demagnification (7.24254), center (1042.99 pix, 1035.88 pix), and rotation (0.240 deg) is measured by ccs in January and Feb 2011 using pinhole mask data taken during the eighth cooldown. These are described in "README.focal_plane_mapping.txt".

bar pitch is 5.8 mm which is related to pixels using the demagnification and temperature scale.

'''
import Detector, IO
import numpy as np
import unittest

class MismatchError(Exception):
        '''The code expected a CSU with 46 slits, but found something else.'''

        def __init__(self, value):
                self.parameter = value

numslits = 46
numbars = numslits * 2

tempscale = 0.99646 

def mm_to_pix(mm):
        return mm/demagnification/Detector.pixelsize

mm = 1
rotation = np.radians(0.2442)
demagnification = 7.24254
center = (1042.99, 1035.88)
barpitch_mm = (5.8 * mm * tempscale)/demagnification
barpitch_pix = mm_to_pix(5.8 * mm * tempscale)

def csu_mm_to_pix(x_mm, slitno):
        '''Convert a slit's position into a pixel value. This is a linear approximation to a sixth order polynomial fit by ccs.'''

        # _kfp is keck focal plane
        # Not sure where to apply tempscale
        # The x_kfp has a "fudge factor"
        x_kfp = x_mm*tempscale  - 4.4*1.2 *mm
        y_kfp = 5.8*mm * tempscale * (numslits - slitno + 0.94)

        # _mfp is Mosfire focal plane
        # Convert the mms into pixels
        x_mfp = 2047 - mm_to_pix(x_kfp)
        y_mfp = mm_to_pix(y_kfp)


        # Rotate around the center
        x_mfp -= center[0]
        y_mfp -= center[1]

        x_mfp = np.cos(rotation)*x_mfp - np.sin(rotation)*y_mfp
        y_mfp = np.sin(rotation)*x_mfp + np.cos(rotation)*y_mfp

        x_mfp += center[0]
        y_mfp += center[1]

        return np.array([x_mfp, y_mfp])
        
def bar_to_slit(x):
        '''Convert a bar #(1-92) to a slit(1-46) number'''
        if (x < 1) or (x > numbars):
                raise MismatchError("Not indexing CSU properly")
        return int(x+1)/2


class Barset:
        '''Barset provides convenience functions around a CSU slitmask'''

        pos = [] 

        def __init__(self):
                pass

        def set_header(self, header):
                '''Passed "header" a FITS header dictionary and converts to a Barset'''
                self.pos = IO.parse_header_for_bars(header)

        def set_mms(self, positions):
                '''Reads a list of [mm] * numbars'''

                if len(positions) != numbars:
                        raise MismatchError("Found %i bars instead of %i" % (lens(poss), numbars))

                self.pos = positions

        def get_bar_pix(self, bar):
                '''Return the pixel position of bar(1-92)'''
                slit = bar_to_slit(bar)
                return csu_mm_to_pix(self.pos[bar-1], slit)


class TestCSUFunctions(unittest.TestCase):

        def setUp(self):
                pass

        def test_bar_to_slit(self):
                sa = self.assertTrue
                sa(bar_to_slit(1) == 1)
                sa(bar_to_slit(2) == 1)
                sa(bar_to_slit(91)==46)
                sa(bar_to_slit(92)==46)
                sa(bar_to_slit(92.)==46)
                sa(bar_to_slit(1.)==1)
                sa(bar_to_slit(1.5)==1)
                sa(bar_to_slit(2.)==1)
                self.assertRaises(MismatchError, bar_to_slit, (-1, 0, 93, 94))

        def test_Barset(self):
                b = Barset()
                pos = np.arange(92)
                b.set_mms(pos)
                self.assertTrue((b.pos == pos).all())

                p1 = b.get_bar_pix(1)
                self.assertTrue((p1==csu_mm_to_pix(pos[0], 1)).all())

                # FIXME More tests here

        
if __name__ == '__main__':
        unittest.main()