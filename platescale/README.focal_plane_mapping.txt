 Pinhole Mask Mapping of Focal Plane 

    C. Steidel, 12 January 2011 
     updates 15 January 2011, with final fidiucal positions; using m110115_0135.fits at
   fiducial FCS corrected imaging position.

  A sequence of ~100 FCS-dithered frames were taken on 7 January (152-253) with the
telescope at zenith and FCS centered around the nominal center "thetax=0.0 thetay=0.0"
which is represented by frames m*0152 and m*0153.  The relative shift of all frames
relative to 152 was then calculated (shifts_vs_152.dat), and the shifts were applied
to the images and a stack registered to the fidicual (frame 0152) was created, called
fcs_shift.avg.0152.fits.  The pinholes were then measured in this stacked image, with
the resulting positions called all_spots_xsort_0152.xy (and all_spots_xsort_0152.reg
is a corresponding ds9 format reg file). Another version of this list was created
which was inversely sorted in detector ypix position to correspond to the order of
measurements on the Excel spreadsheet. The file "pinholes_xy_XYmm.dat" is a text-formatted
file with the corresponding xpix, ypix, Xmm, Ymm where Xmm and Ymm are the actual measured
pinhole positions, where inches were converted to mm (mm = inches*25.40). 
[Note that this was initially done with room-temp CMM measurements, and then repeated
with appropriate values for the pinhole mask at 120k. All summary numbers below 
us the appropriate 120k values] 
 
  ** 15 Jan: used pinholes_xy_XrotYrot_120k.dat as master coords file.** 
***Use fcs_finalshift_pix2mm.db and fcs_finalshift_mm2pix.db, which are summarized
in pix2mm.finalshift.sum and mm2pix.finalshift.sum. 

  In order to convert the focal plane coords into the coordinate system used by the
CSU (i.e., not rotated as the pinhole mask pattern), the Xmm,Ymm coords were rotated
by -5 degrees, the nominal rotation designed into the pinhole mask.  This was done
via the SM script "do_rot" in the macro set plots.sm (included here). 

do_rot    	data pinholes_xy_XYmm.dat
		# rotate Xmm Ymm pinhole positions by -5 degrees to align with CSU coords
		read {xmm 3 ymm 4}
		set theta=-5.0*3.14159265/180.
		set xmm_prime=xmm*cos(theta)-ymm*sin(theta)
		set ymm_prime=xmm*sin(theta)+ymm*cos(theta)
		print XmmYmm_XrotYrot.dat '%10f %10f %10f %10f\n' {xmm ymm xmm_prime ymm_prime}
		print Xrot_Yrot.dat '%10f %10f \n' {xmm_prime ymm_prime}
		

  A new reference file was created with the following format:

   xpix ypix  Xrot Yrot  

where (xpix,ypix) is the measured position of a pinhole image, and (Xrot,Yrot)
is the position of a pinhole in mm using the CSU coordinate frame (i.e., a rotated
version of the CMM frame used for the measurements), centered at (0,0).

  I then used the IRAF task "geomap" to map the two coordinate systems to one another. 
I made one mapping mm to x,y (in the database file fcs_finalshift_mm2pix.db, the other fcs_finalshift_pix2mm.db, 
though one can use the task "geoxtyran" to go either forward or backward using the single solution. 

After some experimentation, I ended up using legendre polynomials in both coords, 
with order 6 and crossterms, and the residuals to the solution are about 6 microns 
rms in the CSU plane, or 0.08 pix going from the CSU plane to the detector. 
(see the summary files, pix2mm.finalshift.sum and mm2pix.finalshift.sum for resume of the
fit).  As part of solving the transformation, one
gets an axis rotation, a translation, and a scale factor which are summarized here:
Note that the numbers quoted are the best linear values, but the actual solution
includes higher order terms (see below for how to use the solutions). 


I repeated exactly the process above with the measurements from
the CSU after scaling the measurements from 293k to 120k, involving applying
a scaling factor of 0.99646 (i.e., measures shrink).  

Coordinate list: fcs_shift.finalshift.XrotYrotxy.120k.dat  Transform: mm2pix.finalshift
    Results file: mm2pix.finalshift.sum
Coordinate mapping status
    Xin and Yin fit rms: 0.04616478  0.05277011
Coordinate mapping parameters
    Mean Xref and Yref: 0.4451433  -3.111243
    Mean Xin and Yin: 1046.501  1012.034
    X and Y shift: 1042.987  1035.884  (xin  yin)
    X and Y scale: 7.670743  7.67076  (xin / xref  yin / yref)
    X and Y axis rotation: 359.75525  359.75653  (degrees  degrees)

Coordinate list: fcs_shift.finalshift.xyXrotYrot.120k.dat  Transform: pix2mm.finalshift
    Results file: pix2mm.finalshift.sum
Coordinate mapping status
    Xin and Yin fit rms: 0.006013579  0.00687888
Coordinate mapping parameters
    Mean Xref and Yref: 1046.501  1012.034
    Mean Xin and Yin: 0.4451433  -3.111243
    X and Y shift: -136.5421  -134.4611  (xin  yin)
    X and Y scale: 0.1303655  0.1303651  (xin / xref  yin / yref)
    X and Y axis rotation: 0.24474  0.24348  (degrees  degrees)

**** New imaging center with FCS correction is (1042.99, 1035.88)  ****

This suggests that 1 pixel maps to 0.1303655 mm in the CSU focal plane, for a net
scale change of 7.24254 (mm at detector to mm in focal plane).   The
accuracy of the solution is +/-6 microns at the CSU, or 0.24 mil, which
is hitting up against the accuracy of the CMM measurements, as expected. 

The axis rotation, 0.245 degrees, is 
just the overall rotation of the image on the detector relative to the CSU coord
axes (as we have previously measured). Note that in spec mode, the spectra are
within about 0.08 degrees of being lined up with the detector rows (we chose to
optimize the detector alignment for spectroscopy). 

The solutions look something like this:

# Sun 14:34:16 16-Jan-2011
begin	mm2pix.finalshift
	xrefmean	0.4451433
	yrefmean	-3.111243
	xmean		1046.501
	ymean		1012.034
	geometry	general
	function	legendre
	xshift		1042.987
	yshift		1035.884
	xmag		7.670743
	ymag		7.67076
	xrotation	359.7552
	yrotation	359.7565
	xrms		0.04616478
	yrms		0.05277011
	surface1	11
			2.	2.
			2.	2.
			2.	2.
			0.	0.
			-137.	-137.
			137.	137.
			-137.	-137.
			137.	137.
			1042.987	1035.884
			1050.882	4.488986
			-4.465763	1050.885
	surface2	29
			2.	2.
			6.	6.
			6.	6.
			2.	2.
			-137.	-137.
			137.	137.
			-137.	-137.
			137.	137.
			-0.005353444	0.04255729
			0.6608887	-0.01530954
			-0.05799868	0.1157489
			0.9264673	-0.01261903
			0.01994258	0.0797478
			0.2227303	0.01168982
			-0.01945649	0.6234929
			0.03669446	0.2193253
			-0.04854986	1.709188
			0.03470791	0.1229889
			-0.04925432	0.3735662
			-0.01216916	0.1628895
			1.833832	-0.02406128
			0.1073147	0.2629326
			1.076982	-0.09890223
			-0.009983247	0.848157
			-0.02049647	0.05725569
			5.596653E-4	0.9405158
			0.06984585	0.0326806
			0.4765024	-0.04174766
			-0.01029157	0.225348

The name of the solution in this case is "mm2pix.finalshift". 

***********************************************
   Closest Linear Solution
***********************************************

 Note that a linear solution allowing for no distortion and a single scale
change and rotation has rms of about 24 microns at the CSU focal plane 
 (i.e., just a rotation and scale change is allowed).  
The mm to pixels transformation has an rms of about 0.18 pixels over
the whole field. The linear solutions have been included in this
directory.  

Coordinate list: fcs_shift.finalshift.xyXrotYrot.120k.dat  Transform: pix2mm.finalshift.linear
    Results file: pix2mm.finalshift.linear.sum
Coordinate mapping status
    Xin and Yin fit rms: 0.02339263  0.02286088
Coordinate mapping parameters
    Mean Xref and Yref: 1046.501  1012.034
    Mean Xin and Yin: 0.4451433  -3.111243
    X and Y shift: -136.5435  -134.4629  (xin  yin)
    X and Y scale: 0.1303655  0.1303655  (xin / xref  yin / yref)
    X and Y axis rotation: 0.24420  0.24420  (degrees  degrees)


Coordinate list: fcs_shift.finalshift.XrotYrotxy.120k.dat  Transform: mm2pix.finalshift.linear
    Results file: mm2pix.finalshift.linear.sum
Coordinate mapping status
    Xin and Yin fit rms: 0.1794415  0.1753713
Coordinate mapping parameters
    Mean Xref and Yref: 0.4451433  -3.111243
    Mean Xin and Yin: 1046.501  1012.034
    X and Y shift: 1042.985  1035.885  (xin  yin)
    X and Y scale: 7.670743  7.670743  (xin / xref  yin / yref)
    X and Y axis rotation: 359.75580  359.75580  (degrees  degrees)


*********************************************************************8

Once we are happy with the solution, it is easy to use it to translate det x,y to
focal plane mm. An example parameter file for geoxytran is attached here. 

        input = "pinholes_xpix_ypix.dat" Input coordinate files to be transformed
       output = "STDOUT"        Output transformed coordinate files
     database = "fcs_finalshift.pix2mm.db" The GEOMAP database file
   transforms = "pix2mm.finalshift"      Names of the coordinate transforms in the database
    (geometry = "geometric")    Transformation type (linear,geometric)
   (direction = "forward")      Transformation direction (forward|backward)
        (xref = INDEF)          X input origin in reference units
        (yref = INDEF)          Y input origin in reference units
        (xmag = INDEF)          X scale in output units per reference unit
        (ymag = INDEF)          Y scale in output units per reference unit
   (xrotation = INDEF)          X axis rotation in degrees
   (yrotation = INDEF)          Y axis rotation in degrees
        (xout = INDEF)          X output origin in output units
        (yout = INDEF)          Y output origin in output units
      (xshift = INDEF)          X origin shift in output units
      (yshift = INDEF)          Y origin shift in output units
     (xcolumn = 1)              Input column containing the x coordinate
     (ycolumn = 2)              Input column containing the y coordinate
    (calctype = "real")         Data type for evaluation coordinates
     (xformat = "")             Output format of the x coordinate
     (yformat = "")             Output format of the y coordinate
(min_sigdigit = 7)              Minimum precision of output x and y coordinates
        (mode = "ql")           

 This file will take a list of (x y) pairs, one per line, and print out to STDOUT
the predicted focal plane position in Xmm,Ymm of the proper CSU plane. 
