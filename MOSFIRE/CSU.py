

''' 
MOSFIRE CSU Utlity Code
        Includes physical parameters of CSU.

Created March 2, 2011 by npk

numslits, numbars give information about the number of slits and bars in the CSU.

Note, for consistency with the hardware bars and slits are indexed from 1.

tempscale is the thermal scaling factor to shrink room temperature linear dimensions to 120 K. The number 0.99646 is from R. Weber's spreadsheet "MOSFIRE Thermal Dimension Scaling Factors.xls", referenced from "Thermophysical properties of matter, Vol 12-13, Thermal Expansion".

demagnification (7.24254), and center_pix (1042.99 pix, 1035.88 pix) is measured by ccs in January and Feb 2011 using pinhole mask data taken during the eighth cooldown. These are described in "README.focal_plane_mapping.txt".

bar pitch is 5.8 mm which is related to pixels using the demagnification and temperature scale.

'''
import Detector, IO
import numpy as np
import unittest
from pyraf import iraf


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
demagnification = 7.24254
center_pix = (1042.986, 1035.879)
barpitch_mm = (5.8 * mm * tempscale)/demagnification
barpitch_pix = mm_to_pix(5.8 * mm * tempscale)

def in_field(px, py):
        '''Determines if the pixel coordinate (x,y) is within the circular 1150 pix FOV'''

        x = px-center_pix[0]
        y = py-center_pix[1]
        dist = np.sqrt(x*x + y*y)
        if dist < 1150.: return True
        return False

def csu_pix_to_mm_poly(x_pix, y_pix):
        (x_kfp, y_kfp) = mosfire_geoxytran(x_pix, y_pix, direction="backward")
        centerx = 137.400

        x_mm = centerx - x_kfp
        y_mm = y_kfp

        return (x_mm, y_mm)

def csu_mm_to_pix_poly(x_mm, slitno):
        '''Uses ccs fits in ../platescale directory'''
        # _kfp is keck focal plane
        centerx = 137.400
        x_kfp = (centerx - x_mm) 
        y_kfp = 5.8 * (numslits/2. - slitno + 0.35)  * tempscale


        return mosfire_geoxytran(x_kfp, y_kfp)


def csu_mm_to_pix(x_mm, slitno, Please_Use=False):
        '''Convert a slit's position into a pixel value. This is a linear approximation to a sixth order polynomial fit by ccs.
        Positions are index from 1: 1 .. 2048
        '''

        if Please_Use==False:
                raise Exception("Use csu_mm_to_pix_poly (a polynomial fit) rather than csu_mm_to_pix (a linear fit)")
                return

        # _kfp is keck focal plane
        centerx = 137.400
        x_kfp = (centerx - x_mm) * tempscale
        y_kfp = 5.8*mm * (numslits/2. - slitno + 0.35) * tempscale

        #
        return  mosfire_geoxytran(x_kfp, y_kfp, database="../platescale/linear_pix2mm_120k.db", transform="linear_pix2mm_120k")


def mosfire_geoxytran(x_kfp, y_kfp, transform="final.pix2mm.4.972.120k", database="../platescale/10March2011.4.972.db", direction="forward"):
        '''Conveninece wrapper around IRAF geoxytran'''
        iraf.images()

        pars = iraf.geoxytran.getParList()
        iraf.geoxytran.unlearn()
        t = iraf.geoxytran("STDIN", "STDOUT", Stdin=["%f %f" % (x_kfp, y_kfp)], Stdout=1,
                        database=database,
                        transform=transform,
                        direction=direction)

        iraf.geoxytran.setParList(pars) 
        (x,y) = np.array(t[0].split(), dtype=np.float64)

        return (x,y)

def mosfire_geoxytrans(x_kfp, y_kfp, transform="final.pix2mm.4.972.120k", database="../platescale/10March2011.4.972.db", direction="forward"):
        '''Conveninece wrapper around IRAF geoxytran'''
        iraf.images()

        pars = iraf.geoxytran.getParList()
        iraf.geoxytran.unlearn()
        ins = []
        for i in range(len(x_kfp)):
                ins.append("%f %f" % (x_kfp[i], y_kfp[i]))
        results = iraf.geoxytran("STDIN", "STDOUT", Stdin=ins, Stdout=1,
                        database=database,
                        transform=transform,
                        direction=direction)

        iraf.geoxytran.setParList(pars) 

        poss = []
        for result in results:
                poss.append(np.array(result.split(), dtype=np.float64))

        return np.array(poss)

        
def bar_to_slit(x):
        '''Convert a bar #(1-92) to a slit(1-46) number'''
        if (x < 1) or (x > numbars):
                raise MismatchError("Not indexing CSU properly")
        return int(x+1)/2

def to_ds9_region(poss, dash=1, color="green", label=True):
        s = []
        d = np.radians(4.2)
        dx = barpitch_pix/2. * np.sin(d)
        dy = barpitch_pix/2. * np.cos(d)

        for i in range(1,numbars+1):
                pos = poss[i-1]
                if not np.isfinite(pos[0]): continue
                if not np.isfinite(pos[1]): continue

                ln = [pos[0]+dx,pos[1]-dy,pos[0]-dx,pos[1]+dy]
                if label:
                        s.append("line(%6.3f, %6.3f, %6.3f, %6.3f) # line=0 0 color=%s text={b%2.0i} dash=%1i fixed=1 edit=0 move=0 rotate=0 \n" % (ln[0], ln[1], ln[2], ln[3], color, i, dash))
                else:
                        s.append("line(%6.3f, %6.3f, %6.3f, %6.3f) # line=0 0 color=%s dash=%1i fixed=1 edit=0 move=0 rotate=0 \n" % (ln[0], ln[1], ln[2], ln[3], color, dash))

        return s


class Barset:
        '''Barset provides convenience functions around a CSU slitmask'''

        pos = [] 
        pos_pix = []

        def __init__(self):
                pass

        def set_header(self, header):
                '''Passed "header" a FITS header dictionary and converts to a Barset'''
                self.pos = np.array(IO.parse_header_for_bars(header))
                self.set_pos_pix()

        def set_mms(self, positions):
                '''Reads a list of [mm] * numbars'''

                if len(positions) != numbars:
                        raise MismatchError("Found %i bars instead of %i" % (lens(poss), numbars))

                self.pos = np.array(positions)
                self.set_pos_pix()

        def set_pos_pix(self):
                # _kfp is keck focal plane
                centerx = 137.400
                x_kfp = (centerx - self.pos) 
                slitno = np.ceil(np.arange(1, numbars+1)/2.)
                y_kfp = 5.8 * (numslits/2. - slitno + 0.35)  * tempscale

                self.pos_pix = mosfire_geoxytrans(x_kfp, y_kfp)


        def to_ds9_region(self):
                poss = []

                for i in range(1,numbars+1):
                        poss.append(self.get_bar_pix(i))

                return to_ds9_region(poss)

        def get_bar_pix(self, bar):
                '''Return the pixel position of bar(1-92)'''
                return self.pos_pix[bar-1]


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

        def test_bar_mm(self):
                sa = self.assertTrue


                # Values are taken from ccs
                p0 = mosfire_geoxytran(0,0)
                sa(np.abs(p0[0] - center_pix[0]) < 1e-6)
                sa(np.abs(p0[1] - center_pix[1]) < 1e-6)

        def test_Barset(self):
                b = Barset()
                pos = np.arange(92)
                b.set_mms(pos)
                self.assertTrue((b.pos == pos).all())

                p1 = b.get_bar_pix(1)
                #self.assertTrue((p1==csu_mm_to_pix(pos[0], 1)).all())

                # FIXME More tests here

        
if __name__ == '__main__':
        unittest.main()
