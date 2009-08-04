/**
 * 
 * MASCGEN: MOSFIRE Automatic Slit Configuration GENerator. 
 * See help text (execute with "help" argument) for usage, input/output, and
 * options. Mascgen computes the slit configuration that maximizes total 
 * priority. It creates a slit list, slit region file, and bar position list.
 * Author: Christopher Klein
 * Summer 2007 
 * 
 * **/

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.*;


public class Mascgen_Delta {
	// First, define some dimensions of the CSU.
	// We use dimensions as given in mm, but then multiply by the conversion
	// factor.
	final static double MM_TO_AS = 1.37896;
	final static double SINGLE_SLIT_HEIGHT = 5.1 * MM_TO_AS;
	final static double OVERLAP = 0.7 * MM_TO_AS;
	// The CSU has 46 rows for slits, spaced by 45 "full" overlap regions, 
	// with one extra "half" overlap region on the very top and another 
	// "half" overlap region on the very bottom. Thus, we compute the total
	// CSU height: 46 * SINGLE_SLIT_HEIGHT + 45 * OVERLAP + OVERLAP / 2 + 
	// OVERLAP / 2 = 46 * (SINGLE_SLIT_HEIGHT + OVERLAP).
	final static double CSU_HEIGHT = 46 * (SINGLE_SLIT_HEIGHT + OVERLAP);
	// We assume the CSU to be square, with width = height.
	final static double CSU_WIDTH = CSU_HEIGHT;
	// Set the default slit width, for use if invalid slit width is entered.
	final static double DEF_SLIT_WIDTH = 0.7;
	// Set the CSU Focal Plane radius.
	final static double CSU_FP_RADIUS = 204;
	// Set the bar tilt angle in degrees.
	final static double BAR_TILT = 4;
	
    static Region orArray[] = new Region[47];
    static Region rrArray[] = new Region[46];
    static RaDec finalFieldCenter = new RaDec();
    static double finalPA;
    static AstroObj[] astroObjArrayOriginal;
    static AstroObj[] optimumAstroObjArray;
	
	
	
	// Blank constructor for mascgen_beta.
	public Mascgen_Delta() {
	}
	
	
	// The optimize method takes in an array of AstroObjs and the user-defined
	// CSU Parameters, center position (in right ascension and declination),
	// and position angle (in degrees). It returns an array of AstroObjs (of 
	// length less than or equal to 46) which compose the slit mask 
	// configuration with the higest total priority. If this array is passed to
	// the method slitConfigurationGenerator, it will return the corresponding
	// optimal slit configuration.
	public AstroObj[] optimize(AstroObj[] astroObjArrayOriginal,
			CSUParameters verifiedCSUData, 
			RaDec centerPosition, 
			double pa) {
		
        double minLegalX = verifiedCSUData.getMinLegalX();
        double maxLegalX = verifiedCSUData.getMaxLegalX();
        double xCenter = verifiedCSUData.getXCenter();
        
        // "Hard copy" the input AstroObj array so that subsequent optimize
        // calls are disrupted by any changes to the original 
        // astroObjArrayOriginal as read in.
        AstroObj[] astroObjArray = new AstroObj[astroObjArrayOriginal.length];
        for  (int i = 0; i < astroObjArray.length; i++) {
        	astroObjArray[i] = new AstroObj();
        	astroObjArray[i].setObjName(astroObjArrayOriginal[i].getObjName());
        	astroObjArray[i].setObjPriority(
        			astroObjArrayOriginal[i].getObjPriority());
        	astroObjArray[i].setObjMag(astroObjArrayOriginal[i].getObjMag());
        	astroObjArray[i].setRaHour(astroObjArrayOriginal[i].getRaHour());
        	astroObjArray[i].setRaMin(astroObjArrayOriginal[i].getRaMin());
        	astroObjArray[i].setRaSec(astroObjArrayOriginal[i].getRaSec());
        	astroObjArray[i].setDecDeg(astroObjArrayOriginal[i].getDecDeg());
        	astroObjArray[i].setDecMin(astroObjArrayOriginal[i].getDecMin());
        	astroObjArray[i].setDecSec(astroObjArrayOriginal[i].getDecSec());
        	astroObjArray[i].setEpoch(astroObjArrayOriginal[i].getEpoch());
        	astroObjArray[i].setEquinox(astroObjArrayOriginal[i].getEquinox());
        	astroObjArray[i].setJunk1(astroObjArrayOriginal[i].getJunk1());
        	astroObjArray[i].setJunk2(astroObjArrayOriginal[i].getJunk2());
        	astroObjArray[i].setWcsX(astroObjArrayOriginal[i].getWcsX());
        	astroObjArray[i].setWcsY(astroObjArrayOriginal[i].getWcsY());   
        }
        
        // Transform the entire astroObjArray into the CSU plane by subtracting
        // the center coordinate from each AstroObj's xCoordinate and 
        // yCoordinate and putting these into the ObjX and ObjY.
        for (int i = 0; i < astroObjArray.length; i++) {
//        	 Rotate the objects in the CSU plane by the Position Angle.
            /** Objects were read in with coordinate system origin at center of 
             *  CSU field. The optimize method runs with the coordinate system 
             *  origin in the lower left. So, simply add CSU_WIDTH / 2 to the x  
             *  position and CSU_HEIGHT / 2 to the y position of each object. **/
        	double xOld = astroObjArray[i].getWcsX() 
			- centerPosition.xCoordinate;
        	double yOld = astroObjArray[i].getWcsY() 
			- centerPosition.yCoordinate;
        	double theta = pa * Math.PI / 180;
        	astroObjArray[i].setObjX(xOld * Math.cos(theta) 
        			- yOld * Math.sin(theta) + CSU_WIDTH / 2);
        	astroObjArray[i].setObjY(xOld * Math.sin(theta) 
        			+ yOld * Math.cos(theta) + CSU_HEIGHT / 2);
        }
        double circleOriginX = CSU_WIDTH / 2;
        double circleOriginY = CSU_HEIGHT / 2;
        
//	    for (int i = 0; i < astroObjArray.length; i++) {
//    	astroObjArray[i].printObj();
//	    }
//	    
//	    System.out.println();

	    ArrayList<AstroObj> astroObjArrayList = new ArrayList<AstroObj>();
        /** Crop out all AstroObjs in the astroObjArray that lie outside the 
         * focal plane circle, defined by CSU_FP_RADIUS centered at the origin
         * (CSU_WDITH / 2, CSU_HEIGHT / 2). **/
        /** Crop out all AstroObjs in the astroObjArray that have x coordinate
         * positions outside of the legal range. **/
        /** Crop out all AstroObjs in the astroObjArray that have x coordinate
         * positions outside of the CSU Plane. **/
	    for (int i = 0; i < astroObjArray.length; i++) {
	    	if (isInCircle(astroObjArray[i], circleOriginX, circleOriginY) && 
	    			isLegalInX(astroObjArray[i], minLegalX, maxLegalX) && 
	    			(astroObjArray[i].getObjX() > 0) && 
        			(astroObjArray[i].getObjX() < CSU_WIDTH)) {
	    		astroObjArrayList.add(astroObjArray[i]);
	    	}
	    }
	    /** Assign each AstroObj to its correct RowRegion or OverlapRegion.
         * Then crop out the objects that fall into OverlapRegions 0 or 46.
         * These are the regions at the very top and bottom, and objects in 
         * these regions cannot be assigned slits. **/
	    for (int i = 0; i < astroObjArrayList.size(); i++) {
	    	for (int j = 0; j < 47; j++){
        		if ((orArray[j].getMinY() < astroObjArrayList.get(i).getObjY()) &&
        				(astroObjArrayList.get(i).getObjY() < orArray[j].getMaxY()))
        				astroObjArrayList.get(i).setObjOR(j);
            	if (j < 46 && 
            			(rrArray[j].getMinY() < astroObjArrayList.get(i).getObjY()) &&
            			(astroObjArrayList.get(i).getObjY() < rrArray[j].getMaxY()))
            			astroObjArrayList.get(i).setObjRR(j);
            	}
        	if ((astroObjArrayList.get(i).getObjOR() == 0) || 
        			(astroObjArrayList.get(i).getObjOR() == 46)) {
        		astroObjArrayList.remove(i);
        	}
	    }
	    AstroObj[] astroObjArray4 = 
	        astroObjArrayList.toArray(new AstroObj[astroObjArrayList.size()]);
        
        // Thus far, MASCGEN has cropped out the following objects from the 
        // input ObjectList: 
        // 1) Objects outside the CSU Focal Plane circle
        // 2) Objects outside the legal x-coordinate range (the exact range was
        //    specified by the user)
        // 3) Objects in OverlapRegions 0 or 46
        
        /** Scan through each RowRegion and OverlapRegion and find the object(s)
         * with the highest priority in each each region. This is essentially
         * the first pass algorithm, and the result is highPriorityArray1. **/
        AstroObj[] highPriorityArray1 = new AstroObj[91];
        for (int j = 0; j < 46; j++) {
        	highPriorityArray1[2 * j] = 
        		getHighestPriorityInRow(astroObjArray4, j, xCenter);
        	if (j != 45) {
        		AstroObj blankTemp = new AstroObj();
        		blankTemp.setObjOR(j + 1);
        		highPriorityArray1[2 * j + 1] = blankTemp;
        	}
        }
        // Perform a similar operation on the overlap region objects.
        AstroObj[] overlapRegionHighPriorityArray1 = new AstroObj[45];
        for (int j = 1; j < 46; j++) {
        	overlapRegionHighPriorityArray1[j - 1] = 
        		getHighestPriorityInOverlap(astroObjArray4, j, xCenter);
        }
        
		/** Run the second pass algorithm. **/
        // Create the "second pass" array by comparing the sum of the 
        // priorities of objects in two neighboring rows with the priority of 
        // the object in the overlap region in the middle. If the object in the
        // overlap region's priority is larger, "link" the two row regions to 
        // make a double-length slit. This sacrifices the two objects in the 
        // rows for the one in the overlap, but it increases the total priority
        // of the configuration.
        // First, sort overlapRegionHighPriorityArray1 by priority and closeness 
		// to x-center to create array oRHPA1s.
        List<AstroObj> overlapRegionHPList1 = 
        	Arrays.asList(overlapRegionHighPriorityArray1);
        Collections.sort(overlapRegionHPList1, 
        		new PriorityAndXOrderCompare(xCenter));
        AstroObj oRHPA1s[] = 
        	overlapRegionHPList1.toArray(new AstroObj[]{});
        // Then, do the comparison and insertion into the new hPArray2.
        // The employed method here is the simple nearest-neighbors-only linkage
        // scheme.
        AstroObj blank = new AstroObj();
        AstroObj hPArray2[] = new AstroObj[highPriorityArray1.length];
        System.arraycopy(highPriorityArray1, 0, hPArray2, 0,
        		highPriorityArray1.length);
        AstroObj[] rowRegionhighPriorityArray1 = new AstroObj[46];
        int rowRegionhighPriorityArray1Index = 0;
        for (int i = 0; i < highPriorityArray1.length; i++) {
        	if ((highPriorityArray1[i].isNotBlank()) || (highPriorityArray1[i].getObjOR() == -1)) {
        		System.arraycopy(highPriorityArray1, i, rowRegionhighPriorityArray1, 
        				rowRegionhighPriorityArray1Index,1);
        		rowRegionhighPriorityArray1Index++;
        	}
        }
        // Here starts the big for loop with if statements to compare and insert
        // when doing so raises total priority.
        for (int i = 0; i < oRHPA1s.length; i++) {
        	int oRNumber = oRHPA1s[i].getObjOR();
        	if ((oRNumber < 0) || (oRNumber > 45)) {
        		// There's no good object in the OR (it's filled with "blank"),
        		// so skip it.
        		continue;
        	}
        	if ((oRHPA1s[i].getObjPriority() > 
    			(hPArray2[(oRNumber - 1) * 2].getObjPriority() + 
    			hPArray2[oRNumber * 2].getObjPriority())) && 
    			oRNumber == 1) {
        		if (hPArray2[oRNumber * 2 + 1].isBlank()) {
    					hPArray2[oRNumber * 2 - 1] = oRHPA1s[i];
    					hPArray2[oRNumber * 2 - 2] = blank;
    					hPArray2[oRNumber * 2] = blank;
    					hPArray2[oRNumber * 2 + 1] = blank;
        		}
        	}
        	if ((oRHPA1s[i].getObjPriority() > 
				(hPArray2[(oRNumber - 1) * 2].getObjPriority() + 
				hPArray2[oRNumber * 2].getObjPriority())) && 
				oRNumber == 45) {
        		if (hPArray2[oRNumber * 2 - 3].isBlank()) {
					hPArray2[oRNumber * 2 - 1] = oRHPA1s[i];
					hPArray2[oRNumber * 2 - 2] = blank;
					hPArray2[oRNumber * 2] = blank;
					hPArray2[oRNumber * 2 - 3] = blank;
    			}
        	}
        	if (oRNumber > 1 && oRNumber < 45 && (oRHPA1s[i].getObjPriority() > 
        		(hPArray2[(oRNumber - 1) * 2].getObjPriority() + 
        		hPArray2[oRNumber * 2].getObjPriority())) && 
        		(hPArray2[oRNumber * 2 + 1].isBlank() && 
        		hPArray2[oRNumber * 2 - 3].isBlank())) {
        			hPArray2[oRNumber * 2 - 1] = oRHPA1s[i];
        			hPArray2[oRNumber * 2 - 2] = blank;
        			hPArray2[oRNumber * 2] = blank;
        			hPArray2[oRNumber * 2 - 3] = blank;
        			hPArray2[oRNumber * 2 + 1] = blank;
        	}		
        }
        // Then, copy over the non-blank objects into the high-priority 
        // hPArray2Short array. This is what is eventually returned by optimize.
        int finalNumTargets = 0;
        for (int i = 0; i < hPArray2.length; i++) {
        	if (hPArray2[i].isNotBlank()) {
        		finalNumTargets++;
        	}
        }
        AstroObj[] hPArray2Short = new AstroObj[finalNumTargets];
        int hPArray2ShortIndex = 0;
        for (int i = 0; i < hPArray2.length; i++) {
        	if (hPArray2[i].isNotBlank()) {
        		System.arraycopy(hPArray2, i, hPArray2Short, hPArray2ShortIndex, 
        				1);
        		hPArray2ShortIndex++;
        	}
        }
        // hPArray2 is the object list of all the objects that are included into
        // the final slit configuration. The total priority of hPArray2 is the
        // total priority of the final slit configuration. hPArray2Short just 
        // has the blanks cut out. Both hPArray2 and hPArray2Short contain the
        // same data, but the later is easier to read, and it is the input used
        // by the remainder of the program to expand slits where possible and 
        // produce the final slit configuration.
        // But, before returning hPArray2Short, optimize has to transform back
        // from it's lower-left origin coordinate system to the center of CSU
        // origin coordinate system.
        for (int i = 0; i < hPArray2Short.length; i++) {
        	hPArray2Short[i].
        	setObjX(hPArray2Short[i].getObjX() - CSU_WIDTH / 2);
        	hPArray2Short[i].
        	setObjY(hPArray2Short[i].getObjY() - CSU_HEIGHT / 2);
        }
        return hPArray2Short;
	}

	
	// Make the bar list with bar x-positions for each row. There is a total of
	// 92 bars in the MOSFIRE CSU, odd-numbered bars extend from the left and 
	// even-numbered bars extend from the right.
	private double[] barPositionArray = new double[92];                              
	public double[] getBarPositionArray() {
		return barPositionArray;
	}
	
	
	/** Additonal Functions **/
	// Print the help text.
	private static void PrintHelpText() {
		java.text.DecimalFormat thirdPlace = 
			new java.text.DecimalFormat("0.000");
		System.out.println("MASCGEN: MOSFIRE Automatic Slit " +
				"Configuration GENerator" +
				"\n\nUsage: java mascgen_beta [object " +
				"list] [x range] [x center] [slit width] \n\t[dither space] " +
				"[center RaHour] [center RaMin] [center RaSec] " +
				"\n\t[center DecDeg] [center DecMin] [center DecSec]" +
				"\n\t[num x steps] [x step size] [num y steps] [y step size]" +
				"\n\t[center position angle] [num pa steps] [pa step size] " +
				"\n\t[slit list] [slit region file] [bar position list]" +
				"\n\t\tOr, omit the center Ra and Dec coordinates and Mascgen " 
				+ "\n\t\twill compute and use the weighted-priority center " +
				"position. " +
				"\n\nobject list: input text file of target objects. The " +
				"columns (separated by\n\twhitespace) are 1) object name, " 
				+ "2) object priority, 3) object magnitude \n\t4) object " +
				"RaHour, 5) object RaMin, 6) object RaSec, 7) object DecDeg, " +
				"\n\t8) object " +
				"DecMin, 9) object DecSec, 10) epoch, 11) equinox, 12) junk1, "
				+ "\n\tand 13) junk2. " +
				"Do not use priority values < or = 0.\n" +
				"x range: set the width of the legal x coordinate range in " +
				"arc minutes (double)." +
				"\nx center: set the center of the legal x coordinate range " +
				"in arc minutes \n\t(double). \nslit width: set the global slit " +
				"width in arc minutes (this is applied to all\n\tgenerated " +
				"slits) (double). \ndither space: set the minimum distance in " 
				+ "arc seconds from the top/bottom of the\n\tslit an object may"
				+ " legally be placed (this will be the \"buffer zone\" " +
				"\n\tthat permits dithering or nodding without the object" +
				" moving behind a \n\tbar) (double). " +
				"\ncenter RaHour: desired center position Right Ascension " +
				"hours (integer)." +
				"\ncenter RaMin: desired center position Right Ascension " +
				"minutes (integer)." +
				"\ncenter RaSec: desired center position Right Ascension " +
				"seconds (double)." +
				"\ncenter DecDeg: desired center position Declination " +
				"degrees (integer)." +
				"\ncenter DecMin: desired center position Declination arc " +
				"minutes (integer)." +
				"\ncenter DecSec: desired center position Declination arc " +
				"seconds (double)." +
				"\nnum x steps: number of x iterations to perform above and " +
				"below the center " +
				"\n\tposition (integer). The actual number of calculated " +
				"iterations is \n\t2 * (num x steps) + 1." +
				"\nx step size: size of each x step in arc minutes (double)." + 
				"\nnum y steps: number of y iterations to perform above and " +
				"below the center " +
				"\n\tposition (integer). The actual number of calculated " +
				"iterations is \n\t2 * (num y steps) + 1." + 
				"\ny step size: size of each y step in arc minutes (double)." +
				"\ncenter position angle: desired center position angle in " +
				"degrees (double)." +
				"\nnumber pa steps: number of position angle iterations to " +
				"perform above and below " +
				"\n\tthe center position (integer). The actual number of " +
				"calculated \n\titerations is 2 * (num pa steps) + 1." + 
				"\npa step size: size of each pa step in degrees (double)." +
				"\nslit list: name of file to write the ouput slit list to." +
				"\nslit region file: name of file to write the ouput SAOImage "
				+ "Ds9 region file to." +
				"\nbar position list: output list of CSU coordinates (in mm) " +
				"for all 92 bars. " +
				"\n\tOdd-numbered bars extend from the left, even from " +
				"the right.");
		System.out.println("\nPlease note that the CSU total width is " +
				thirdPlace.format(CSU_WIDTH / 60) + " arc minutes (restricts x "
				+ "\n\trange and slit width parameters) and that maximum " +
				"slit size is " +
				thirdPlace.format(SINGLE_SLIT_HEIGHT) + 
				"\n\tarc seconds between adjacent slits " +
				"(dither space can be at most half \n\tthis). " + 
				"Also, x center cannot be < " + 
				thirdPlace.format(-CSU_WIDTH / 120) + " or > " + 
				thirdPlace.format(CSU_WIDTH / 120) + " arc minutes." + 
				"\n\tEntries that do not conform " +
    			"to these restrictions will cause the program\n\tto " +
    			"exit with an error message.");
		System.out.println("\nThe output of the program is the slit list and " +
				"slit region file. The columns of \n\tthe slit list " +
				"(separated by whitespace) are " +
				"1) slit number, 2) slit \n\tcenter RaHour, 3) slit center " +
				"RaMin, 4) slit center RaSec, 5) slit \n\tcenter DecDeg, 6) " +
				"slit center DecMin, 7) slit center DecSec, 8) slit \n\twidth, "
				+ "9) slit length, 10) object name, 11) object priority, 12) " +
				"object \n\ttarget location (vertical offset from slit " +
				"center in arc seconds), " +
				"\n\t13) object RaHour, 14) object RaMin, 15) object RaSec, " +
				"16) object \n\tDecDeg, 17) object DecMin, and 18) object " +
				"DecSec. The slit region file " +
				"\n\tcan be loaded into SAOImage Ds9 on top of " +
				"an image of the field the " +
				"\n\tinput object list was created from.");
		System.exit(1);
	}
	
	//Read in file and return array of AstroObjs.
	private static AstroObj[] readInFile(String fileName) {
		Scanner fillIn = null;
	    String tempObjName;
	    double tempObjPriority;
	    double tempObjMag;
	    double tempObjRaHour;
	    double tempObjRaMin;
	    double tempObjRaSec;
	    double tempObjDecDeg;
	    double tempObjDecMin;
	    double tempObjDecSec;
	    double tempObjEpoch;
	    double tempObjEquinox;
	    double tempObjJunk1;
	    double tempObjJunk2;
	    ArrayList<AstroObj> astroObjArrayList = new ArrayList<AstroObj>();
	    try {
	    	// First, read in the ObjectList.txt file.
	        fillIn = new Scanner(new BufferedReader(new FileReader(fileName)));   
	    while (fillIn.hasNext()) {
	       	tempObjName = fillIn.next();
	       	tempObjPriority = Double.parseDouble(fillIn.next());
	       	tempObjMag = Double.parseDouble(fillIn.next());
	       	tempObjRaHour = Double.parseDouble(fillIn.next());
	       	tempObjRaMin = Double.parseDouble(fillIn.next());
	       	tempObjRaSec = Double.parseDouble(fillIn.next());
	       	tempObjDecDeg = Double.parseDouble(fillIn.next());
	       	tempObjDecMin = Double.parseDouble(fillIn.next());
	       	tempObjDecSec = Double.parseDouble(fillIn.next());
	       	tempObjEpoch = Double.parseDouble(fillIn.next());
	       	tempObjEquinox = Double.parseDouble(fillIn.next());
	       	tempObjJunk1 = Double.parseDouble(fillIn.next());
	       	tempObjJunk2 = Double.parseDouble(fillIn.next());
	       	astroObjArrayList.add(new AstroObj(tempObjName, 
	       				                       tempObjPriority, 
	       				                       tempObjMag, 
	       				                       tempObjRaHour,
	       				                       tempObjRaMin,
	       				                       tempObjRaSec,
	       				                       tempObjDecDeg,
	       				                       tempObjDecMin,
	       				                       tempObjDecSec,
	       				                       tempObjEpoch,
	       				                       tempObjEquinox,
	       				                       tempObjJunk1,
	       				                       tempObjJunk2));
	    }    
	    {
	        if (fillIn != null)
	            fillIn.close(); // Clean up and close.
	    }
	    } catch (FileNotFoundException e) {
			e.printStackTrace();
		}   
	    AstroObj[] astroObjArray = 
	        astroObjArrayList.toArray(new AstroObj[astroObjArrayList.size()]);
	    return astroObjArray;
	}
	
	// Test AstroObjs in the astroObjArray for position inside the CSU Focal
	// Plane circle.
	private boolean isInCircle(AstroObj obj, double circleOriginX, 
			double circleOriginY) {
		double xDist = Math.abs(obj.getObjX() - circleOriginX);
		double yDist = Math.abs(obj.getObjY() - circleOriginY);
		double diagonalDistance = Math.sqrt(Math.pow(xDist, 2) + 
				Math.pow(yDist, 2));
		if (diagonalDistance < CSU_FP_RADIUS)
			return true;
		else
			return false;
	}
	
	// Test AstroObjs in the astroObjArray2 for legality in the x coordinate.
	// (Does the x position of each object fall in the legal range specified by
	// the user?)
	private boolean isLegalInX(AstroObj obj, double min, double max) {
		if ((obj.getObjX() >= min) && (obj.getObjX() <= max))
			return true;
		else
			return false;
	}
	
	// Sum up the total priorty of objects in an astroObjArray.
	private static double prioritySum(AstroObj array[]) {
		double result = 0;
		for (int i = 0; i < array.length; i++)
			result += array[i].getObjPriority();
		return result;
	}
	
	// Sum up the total weighted RA coordinates of objects in an AstroObjArray.
	public static double raWeightedSum(AstroObj array[]) {
		double result = 0;
		for (int i = 0; i < array.length; i++)
			result += (array[i].getRaHour() * 3600 + 
					array[i].getRaMin() * 60 +
					array[i].getRaSec()) * array[i].getObjPriority();
		return result;
	}
	
	// Sum up the total weighted DEC coordinates of objects in an AstroObjArray.
	public static double decWeightedSum(AstroObj array[]) {
		double result = 0;
		for (int i = 0; i < array.length; i++)
			result += (array[i].getDecDeg() * 3600 + 
					array[i].getDecMin() * 60 +
					array[i].getDecSec()) * array[i].getObjPriority();
		return result;
	}
	
	// Test an AstroObject for lying in a given RowRegion.
	private boolean isInRR(AstroObj obj, int i) {
		if (obj.getObjRR() == i)
			return true;
		else
			return false;
	}
	
	// Test an AstroObject for lying in a given OverlapRegion.
	private boolean isInOR(AstroObj obj, int i) {
		if (obj.getObjOR() == i)
			return true;
		else
			return false;
	}
	
	// Clean up all the Slit Multiple Length numbers and Priorities.
	private void cleanUpSlitMulPri(Slit[] array) {
        for (int i = 0; i < array.length; i++) {
        	if (array[i].getSlitObjPriority() > 0) {
        		for (int j = 0; j < array.length; j++) {
        			if (array[j].getSlitObjName() 
        					== array[i].getSlitObjName()) {
        				array[j].setSlitMul(array[i].getSlitMul());
        				array[j].setSlitObjPriority(array[i].
        						getSlitObjPriority());
        			}
            	}
        	}
        }
	}
	
	// Expand slit1 into slit2 and copy everything over.
	private void expandSlit(Slit s1, Slit s2, double slitWidth) {
		if (s1.getSlitObjName().equals("blank")) {
		}
		else {
			s2.setSlitObjName(s1.getSlitObjName());
	    	s2.setSlitWidth(slitWidth);
	    	s2.setSlitObjX(s1.getSlitObjX());
			s2.setSlitObjY(s1.getSlitObjY());
			s1.setSlitMul(s1.getSlitMul() + 1);
		}
	}
	
	// Return the object with highest priority in specified row region. If there
	// are no objects in the row region, return a blank AstroObj.
    private AstroObj getHighestPriorityInRow(AstroObj[] objectList, int row, 
    		double xCenter) {
    	ArrayList<AstroObj> objsInRow = new ArrayList<AstroObj>();
        for (int i = 0; i < objectList.length; i++) {
        	if (isInRR(objectList[i], row)) {
        		objsInRow.add(objectList[i]);
        		}
        	}
        Collections.sort(objsInRow, new PriorityAndXOrderCompare(xCenter));
        if (objsInRow.isEmpty()) {
        	return new AstroObj();
        }
        else
        	return objsInRow.get(0);
    }
    
	// Return the object with highest priority in specified overlap region. If 
    // there are no objects in the overlap region, return a blank AstroObj.
    private AstroObj getHighestPriorityInOverlap(AstroObj[] objectList, 
    		int overlap, double xCenter) {
    	ArrayList<AstroObj> objsInOverlap = new ArrayList<AstroObj>();
        for (int i = 0; i < objectList.length; i++) {
        	if (isInOR(objectList[i], overlap)) {
        		objsInOverlap.add(objectList[i]);
        		}
        	}
        Collections.sort(objsInOverlap, new PriorityAndXOrderCompare(xCenter));
        if (objsInOverlap.isEmpty()) {
        	return new AstroObj();
        }
        else
        	return objsInOverlap.get(0);
    }
    
    // Produce the slit configuration from a given input array of AstroObjs.
    // This method should be called after optimize has been run, and its input
    // should be the output hPArray2Short from optimize.
    private Slit[] slitConfigurationGenerator(AstroObj[] hPArray2Short, 
    		double slitWidth, double deadSpace, double rowRegionHeight, 
    		RaDec centerPosition, double pa, String barPositionList) {
        /** The rest of the program simply takes the data from hPArray2Short and
         * maps it to a slit configuration. The empty bars are filled
         * by expanding original "singles" into "doubles" or "triples" when
         * possible and by expanding original "doubles" (which were created to
         * include a high-priority object in an overlap region) into "triples"
         * when possible. After this first pass expansion is complete, any
         * remaining empty slits are filled by expanding the nearest slit with
         * the higher-priority object. The slit expansion is conducted so as to 
         * give more vertical space to the objects that are more likely to need 
         * it and then to objects with higher priority. **/
        // First, make a blank Slit Array.
        Slit[] slitArray = new Slit[46];
        for (int i = 0; i < slitArray.length; i++) {
        	slitArray[i] = new Slit();
        	slitArray[i].setSlitNumber(i + 1);
        }
        // Copy the Overlap region objects into their proper slit assignment
        // pair. Object in OverlapRegion# is given slits numbered # and # - 1. 
        // Also copy the Row region objects into their corresponding slit 
        // assignments.
        for (int i = 0; i < hPArray2Short.length; i++) {
        	if (hPArray2Short[i].getObjOR() >= 0) {
        		int OR = hPArray2Short[i].getObjOR();
        		slitArray[OR - 1].setSlitObjName(hPArray2Short[i].getObjName());
        		slitArray[OR].setSlitObjName(hPArray2Short[i].getObjName());
        		slitArray[OR - 1].
        			setSlitObjPriority(hPArray2Short[i].getObjPriority());
        		slitArray[OR - 1].setSlitWidth(slitWidth);
        		slitArray[OR].setSlitWidth(slitWidth);
        		slitArray[OR - 1].setSlitObjX(hPArray2Short[i].getObjX());
        		slitArray[OR].setSlitObjX(hPArray2Short[i].getObjX());
        		slitArray[OR - 1].setSlitObjY(hPArray2Short[i].getObjY());
        		slitArray[OR].setSlitObjY(hPArray2Short[i].getObjY());
        		slitArray[OR - 1].setSlitMul(2);

        	}
        	if (hPArray2Short[i].getObjRR() >= 0) {
        		int RR = hPArray2Short[i].getObjRR();
        		slitArray[RR].setSlitObjName(hPArray2Short[i].getObjName());
        		slitArray[RR].
        			setSlitObjPriority(hPArray2Short[i].getObjPriority());
        		slitArray[RR].setSlitWidth(slitWidth);
        		slitArray[RR].setSlitObjX(hPArray2Short[i].getObjX());
        		slitArray[RR].setSlitObjY(hPArray2Short[i].getObjY());
        		slitArray[RR].setSlitMul(1);
        	}
        }
        cleanUpSlitMulPri(slitArray);
        // Extend slits in length to expand singles into double, doubles into
        // triples, etc.
        // First expand slit 1 into slit 2 and slit 44 into slit 43 if the 
        // first slits are originally singles and if 44 and 43 are unoccupied
        // and if slit 1's object has higher priority than slit 3's (and if 
        // slit 44's object has higher priority than slit 42's).
        // Then, lengthen occupied slits 2 and 45 into 1 and 46, respectively, 
        // if slits 1 and 46 are unoccupied. There is no need to compare 
        // priorities since there can be no objects in slits 0 or 47 (there are
        // no such slits).
        if ((slitArray[1].getSlitMul() == 1) &&
        		(slitArray[2].getSlitMul() == -1) &&
        		(slitArray[1].getSlitObjPriority() > 
        		slitArray[3].getSlitObjPriority())) {
        	expandSlit(slitArray[1], slitArray[2], slitWidth);
        }
        if ((slitArray[44].getSlitMul() == 1) &&
        		(slitArray[43].getSlitMul() == -1) &&
        		(slitArray[44].getSlitObjPriority() > 
        		slitArray[42].getSlitObjPriority())) {
        	expandSlit(slitArray[44], slitArray[43], slitWidth);
        }
        if ((slitArray[1].getSlitMul() == 1)
        		&& (slitArray[0].getSlitMul() == -1)) {
        	slitArray[0].setSlitObjName(slitArray[1].getSlitObjName());
        	expandSlit(slitArray[1], slitArray[0], slitWidth);
        }   
        if ((slitArray[44].getSlitMul() == 1) 
        		&& (slitArray[45].getSlitMul() == -1)) {
        	slitArray[45].setSlitObjName(slitArray[44].getSlitObjName());
        	expandSlit(slitArray[44], slitArray[45], slitWidth);
        }
        cleanUpSlitMulPri(slitArray);
        // Then extend all "singles" into "doubles" (or "triples") when possible 
        // and so that, if there is a conflict, the slit which is extended is 
        // the one that contains the higher-priority object. Note that we at 
        // first ignore the first set of doubles that were created in order to 
        // surround objects in overlap regions. This makes sense since those 
        // original doubles are guaranteed to have their target objects near 
        // their vertical slit center (within the middle overlap region).
        for (int i = 0; i < slitArray.length; i++) {
          	if (slitArray[i].getSlitMul() == 1) {
          		if (i < 44) {
          			if (slitArray[i + 1].getSlitMul() == -1) {
          				if ((slitArray[i].getSlitObjPriority() > 
          				slitArray[i + 2].getSlitObjPriority()) || 
          				(slitArray[i + 2].getSlitMul() == 2)) {
          					expandSlit(slitArray[i], slitArray[i + 1], 
          							slitWidth);
          					slitArray[i + 1].
          					setSlitMul(slitArray[i].getSlitMul());
          				}
          			}
          		}
          		if (i > 1) {
          			if (slitArray[i - 1].getSlitMul() == -1) {
          				if ((slitArray[i].getSlitObjPriority() > 
          				slitArray[i - 2].getSlitObjPriority()) || 
          				(slitArray[i - 2].getSlitMul() == 2)) {
          					expandSlit(slitArray[i], slitArray[i - 1], 
          							slitWidth);
          					slitArray[i - 1].
          					setSlitMul(slitArray[i].getSlitMul());
          				}
          			}
          		}
          	}
        }
        cleanUpSlitMulPri(slitArray);
        // Extend all slits to envelope blanks. After this, every row should be 
        // occupied by a slit.
        for (int j = 0; j < 46; j++) {
        	for (int i = 1; i < slitArray.length - 1; i++) {
        		if ((slitArray[i].getSlitMul() == -1) && 
        				(slitArray[i - 1].getSlitMul() > 0) &&
        				(slitArray[i - 1].getSlitObjPriority() > 
        				slitArray[i + 1].getSlitObjPriority())) {
        			expandSlit(slitArray[i - 1], slitArray[i], slitWidth);
        			cleanUpSlitMulPri(slitArray);
        		}
        		if ((slitArray[i].getSlitMul() == -1) && 
        				(slitArray[i + 1].getSlitMul() > 0) &&
        				(slitArray[i + 1].getSlitObjPriority() > 
        				slitArray[i - 1].getSlitObjPriority())) {
        			expandSlit(slitArray[i + 1], slitArray[i], slitWidth);
        			cleanUpSlitMulPri(slitArray);
        		}
        		
        		if ((slitArray[i].getSlitMul() == -1) && 
        				(slitArray[i + 1].getSlitObjPriority() == 
        					slitArray[i - 1].getSlitObjPriority())) {
        			if ((slitArray[i + 1].getSlitObjY() - deadSpace - (i + 1) 
        					* (SINGLE_SLIT_HEIGHT) - 2 * i * deadSpace) < 
        				(deadSpace + (i + 1) * (SINGLE_SLIT_HEIGHT) + 
        						2 * i * deadSpace - 
        						slitArray[i + 1].getSlitObjY())) {
        				expandSlit(slitArray[i + 1], slitArray[i], slitWidth);
        				cleanUpSlitMulPri(slitArray);
        			}
        			else {
        				expandSlit(slitArray[i - 1], slitArray[i], slitWidth);
        				cleanUpSlitMulPri(slitArray);
        			}
        		}
        	}
        	cleanUpSlitMulPri(slitArray);
        }
        if ((slitArray[45].getSlitMul() == -1) && 
        		(slitArray[44].getSlitMul() > 0)) {
        		expandSlit(slitArray[44], slitArray[45], slitWidth);
        		cleanUpSlitMulPri(slitArray);
        }
        if ((slitArray[0].getSlitMul() == -1) && 
        		(slitArray[1].getSlitMul() > 0)) {
        		expandSlit(slitArray[1], slitArray[0], slitWidth);
        		cleanUpSlitMulPri(slitArray);
        }
        cleanUpSlitMulPri(slitArray);
        // Correct the slitMul values for each slit.
        String tmpObjName;
        int multiple = 0;
        for (int i = 0; i < slitArray.length; i++) {
        	tmpObjName = slitArray[i].getSlitObjName();
        	multiple = 0;
        	for (int j = 0; j < slitArray.length; j++) {
        		if (tmpObjName.equals(slitArray[j].getSlitObjName())) {
        			multiple++;
        		}
        	}
        	slitArray[i].setSlitMul(multiple);
        }
        cleanUpSlitMulPri(slitArray);
        // Calculate the length, y-coordinate, and x-coordinate of each slit.
        for (int i = 0; i < slitArray.length; i++) {
        	if (slitArray[i].getSlitMul() > 0){
        		slitArray[i].setSlitLength(
        				SINGLE_SLIT_HEIGHT * slitArray[i].getSlitMul() 
        				+ OVERLAP * (slitArray[i].getSlitMul() - 1));
        		slitArray[i].setSlitY((deadSpace + rowRegionHeight / 2) 
    					+ i * (rowRegionHeight + 2 * deadSpace) 
    					- CSU_HEIGHT / 2);
        		slitArray[i].setSlitX(slitArray[i].getSlitObjWcsX() - 
        				centerPosition.xCoordinate);
        	}
        }
        // Fill in the barPositionArray with the correct values and write it out
        // to the user-specified file.
        for (int i = 0; i < slitArray.length; i++) {
        	barPositionArray[2*i] = -slitArray[45 - i].getSlitObjX() 
        		- Math.tan(BAR_TILT * Math.PI / 180) 
				* (slitArray[i].getSlitObjY() - slitArray[i].getSlitY()) 
				- slitWidth / 2;
        	barPositionArray[2*i + 1] = -slitArray[45 - i].getSlitObjX() 
        		- Math.tan(BAR_TILT * Math.PI / 180) 
				* (slitArray[i].getSlitObjY() - slitArray[i].getSlitY()) 
				+ slitWidth / 2;
        }
        FileOutputStream out;
		PrintStream p;
		try
		{
			out = new FileOutputStream(barPositionList);
			p = new PrintStream(out);
	    	
		       for (int i = 0; i < barPositionArray.length; i++) {
		        	p.println(barPositionArray[i]);
	    	}
			p.close();
		}
		catch (Exception e)
		{
			System.err.println ("Error writing to file");
		}
        // Make slitArray2, which has none of the blank rows in it. Also, 
		// set all slit numbers to -1 and then .
        int slitArray2Length = 0;
        for (int i = 0; i < slitArray.length; i++) {
        	if (slitArray[i].getSlitMul() > 0) {
        		slitArray2Length++;
        	}
        }
        Slit[] slitArray2 = new Slit[slitArray2Length];
        int slitArray2Index = 0;
        for (int i = 0; i < slitArray.length; i++) {
        	if (slitArray[i].getSlitMul() > 0) {
        		System.arraycopy(slitArray, i, slitArray2, slitArray2Index, 1);
        		slitArray2Index++;
        	}
        	slitArray2[i].setSlitNumber(-1);
        }
        int renumberSlits = 0;
        int newSlitNum = 1;
        while (renumberSlits < slitArray2.length) {
        	slitArray2[renumberSlits].setSlitNumber(newSlitNum);
        	renumberSlits = renumberSlits 
        		+ slitArray2[renumberSlits].getSlitMul();
        	newSlitNum++;
        }
        // Make slitArray3, which is simply the slit information. Each slit
        // is only listed once, not in multiples. The slit length reflects
        // if the slit is a single, double, triple, etc.
        int slitArray3Length = 0;
        for (int i = 0; i < slitArray2.length; i++)
        	if (slitArray2[i].getSlitNumber() > 0) {
        		slitArray3Length++;
        	}
        
        Slit[] slitArray3 = new Slit[slitArray3Length];
        int slitArray3Index = 0;
        for (int i = 0; i < slitArray2.length; i++) {
        	if (slitArray2[i].getSlitNumber() > 0) {
        		System.arraycopy(slitArray2, i, slitArray3, slitArray3Index, 1);
        		slitArray3Index++;
        	}
        }
        // Give slitArray3 correct slitY values and "close" slitX values. Later,
        // the slitX will be modified to account for bar tilt so that the object
        // is placed in the horizontal center of the slit, regardless of its 
        // vertical displacement from center.
        for (int i = 0; i < slitArray3.length; i++) {
        	slitArray3[i].setSlitY(slitArray3[i].getSlitY() 
        			+ (slitArray3[i].getSlitMul() - 1) * 
        			(deadSpace + rowRegionHeight / 2));
        	slitArray3[i].setSlitX(slitArray3[i].getSlitObjX());
        }
        // Rotate the slits in the CSU plane backwards by the Position Angle.
        for (int i = 0; i < slitArray3.length; i++) {
        	double xOld = slitArray3[i].getSlitX();
        	double yOld = slitArray3[i].getSlitY();
        	double theta = -pa * Math.PI / 180;
        	slitArray3[i].setSlitX(xOld * Math.cos(theta) 
        			- yOld * Math.sin(theta));
        	slitArray3[i].setSlitY(xOld * Math.sin(theta) 
        			+ yOld * Math.cos(theta));
        }
        for (int i = 0; i < slitArray3.length; i++) {
        	double xOld = slitArray3[i].getSlitObjX();
        	double yOld = slitArray3[i].getSlitObjY();
        	double theta = -pa * Math.PI / 180;
        	slitArray3[i].setSlitObjX(xOld * Math.cos(theta) 
        			- yOld * Math.sin(theta));
        	slitArray3[i].setSlitObjY(xOld * Math.sin(theta) 
        			+ yOld * Math.cos(theta));
        }
        // Transform slitX and slitY back into wcs coodinates by adding the 
        // center position x and y. Also, calculate the Ra and Dec and reverse 
        // the numbering to conform to top-to-bottom descending convention.
        for (int i = 0; i < slitArray3.length; i++) {
        	slitArray3[i].setWcsX(slitArray3[i].getSlitX() 
        			+ centerPosition.xCoordinate);
        	slitArray3[i].setWcsY(slitArray3[i].getSlitY() 
        			+ centerPosition.yCoordinate);
        	slitArray3[i].setWcsX(slitArray3[i].getWcsX() - 
    				Math.tan(BAR_TILT * Math.PI / 180) 
        			* (slitArray3[i].getSlitObjY() - slitArray3[i].getSlitY()));
        	slitxyToRaDec(slitArray3[i], centerPosition);
        	slitArray3[i].setTargetLocation(slitArray3[i].getSlitObjY() 
        			- slitArray3[i].getSlitY());
        	slitArray3[i].setSlitObjWcsX(slitArray3[i].getSlitObjX() 
        			+ centerPosition.xCoordinate);
        	slitArray3[i].setSlitObjWcsY(slitArray3[i].getSlitObjY() 
        			+ centerPosition.yCoordinate);
        	slitObjxyToRaDec(slitArray3[i], centerPosition);
        	slitArray3[i].setSlitNumber(slitArray3.length 
        			- slitArray3[i].getSlitNumber() + 1);
        }
        return slitArray3;
    }
     
	// Verify that the CSU Field Data conforms to restrictions. If not,
    // throw an exception to notify the user and substitute a default.
    private static CSUParameters verifyCSUData(double xRange, double xCenter, 
    		double slitWidth, double ditherSpace) 
    			throws InvalidArgumentException {
    	java.text.DecimalFormat thirdPlace = 
    		new java.text.DecimalFormat("0.000");
        if (xRange <= 0) {
        	xRange = 3;
        	System.out.println("Error: x range cannot be <= 0 " +
			" arc min. \nSubstituting default x range = 3.");
        }
        if (xRange > CSU_WIDTH / 60) {
        	xRange = 3;
        	System.out.println("Error: x range cannot be > the " +
        			"CSU width, which is " + thirdPlace.format(CSU_WIDTH / 60) +
        			" arc min. \nSubstituting default x range = 3.");
        }
        if (-xCenter <= -CSU_WIDTH / 120) {
        	xCenter = 0;
        	System.out.println("Error: x center cannot be <= " +
        			"-(CSU width) / 2, which is " + 
        			thirdPlace.format(-CSU_WIDTH / 120) + " arc min. " +
        			"\nSubstituting default x center = 0.");
        }
        if (-xCenter >= CSU_WIDTH / 120) {
        	xCenter = 0;
        	System.out.println("Error: x center cannot be >= " +
        			"(CSU width) / 2, which is " + 
        			thirdPlace.format(CSU_WIDTH / 120) + " arc min." +
					"\nSubstituting default x center = 0.");
        }
        if (slitWidth <= 0) {
        	slitWidth = 0.7;
        	System.out.println("Error: slit width cannot be <= 0 " +
        			" arc sec. " +
					"\nSubstituting default slit width = 0.7.");
        }
        if (slitWidth > CSU_WIDTH) {
        	slitWidth = 0.7;
        	System.out.println("Error: slit width cannot be > CSU " +
        			"width, which is " + thirdPlace.format(CSU_WIDTH) + 
        			" arc sec. \nSubstituting default slit width = 0.7.");
        }
        if (ditherSpace < 0) {
        	ditherSpace = 2;
        	System.out.println("Error: dither space cannot be < 0 " +
        			"arc sec." + 
			"\nSubstituting default dither space = 2.");
        }
        if (ditherSpace > SINGLE_SLIT_HEIGHT / 2) {
        	ditherSpace = 2;
        	System.out.println("Error: dither space cannot be >= " +
        			"half of the single slit height, which is " + 
        			thirdPlace.format(SINGLE_SLIT_HEIGHT / 2) + " arc sec." + 
        			"\nSubstituting default dither space = 2.");
        }
        // Calculate some CSU parameters which are more useful.
        double minLegalX = 60 * (xCenter - xRange / 2) + CSU_WIDTH / 2;
        double maxLegalX = 60 * (xCenter + xRange / 2) + CSU_WIDTH / 2;
        double deadSpace = ditherSpace + OVERLAP / 2;
        double rowRegionHeight = SINGLE_SLIT_HEIGHT - 2 * ditherSpace;
        return new CSUParameters(minLegalX, 
        		maxLegalX, deadSpace, rowRegionHeight, xRange, xCenter, 
        		slitWidth, ditherSpace);
        
    }
    
    // Convert the coordinates of an RaDec point from Ra/Dec to x/y.
    private static void raDecToXY(RaDec p) {
		p.setYCoordinate(p.getDecDeg() / Math.abs(p.getDecDeg()) * 
				(3600 * Math.abs(p.getDecDeg()) + 60 * p.getDecMin() + 
				p.getDecSec()));
		p.setXCoordinate(Math.cos(p.getYCoordinate() * Math.PI / 180 / 3600) * 
				15 * 
				(p.getRaHour() * 3600 + 
				p.getRaMin() * 60 +  
				p.getRaSec()));
	}
    
    // Convert the coordinates of an RaDec point from x/y to Ra/Dec.
    private static void xyToRaDec(RaDec p) {
		p.setDecDeg((int) (Math.abs(p.getYCoordinate()) / p.getYCoordinate()) * 
				(int) Math.floor(Math.abs(p.getYCoordinate() / 3600)));		
		p.setDecMin((int) Math.floor((Math.abs(p.getYCoordinate() / 3600) - 
				Math.floor(Math.abs(p.getYCoordinate() / 3600))) * 60));
		p.setDecSec(60 * ((Math.abs(p.getYCoordinate() / 3600) - 
				Math.floor(Math.abs(p.getYCoordinate() / 3600))) * 
				60 - p.getDecMin()));
		
		double decrad = (p.getYCoordinate() * Math.PI / 180 / 3600);
		
		p.setRaHour((int) Math.floor(p.getXCoordinate() / 3600 / 15 / Math.cos(decrad)));
		p.setRaMin((int) Math.floor(60 * ((p.getXCoordinate() / 3600 / 15 / Math.cos(decrad)) - p.getRaHour())));
		p.setRaSec(60 * (((p.getXCoordinate() / 3600 / 15 / Math.cos(decrad)) - p.getRaHour()) * 60 - p.getRaMin()));
	}
	
    // Convert the coordinates of an AstroObj from Ra/Dec to x/y.
    private static void astroObjRaDecToXY(AstroObj obj, RaDec cp) {
		obj.setWcsY(obj.getDecDeg() / Math.abs(obj.getDecDeg()) * 
				(3600 * Math.abs(obj.getDecDeg()) + 60 * obj.getDecMin() + 
				obj.getDecSec()));
		obj.setWcsX(Math.cos(cp.getYCoordinate() * Math.PI / 180 / 3600) * 
				15 * 
				(obj.getRaHour() * 3600 + 
				obj.getRaMin() * 60 +  
				obj.getRaSec()));
	}
	
    // Convert the coordinates of a Slit from x/y to Ra/Dec.
    private static void slitxyToRaDec(Slit slit, RaDec cp) {
		slit.setDecDeg((int) (Math.abs(slit.getWcsY()) / slit.getWcsY()) * 
				(int) Math.floor(Math.abs(slit.getWcsY() / 3600)));
		slit.setDecMin((int) Math.floor((Math.abs(slit.getWcsY() / 3600) - 
				Math.floor(Math.abs(slit.getWcsY() / 3600))) * 60));
		slit.setDecSec(60 * ((Math.abs(slit.getWcsY() / 3600) - 
				Math.floor(Math.abs(slit.getWcsY() / 3600))) * 
				60 - slit.getDecMin()));
		double decrad = (cp.getYCoordinate() * Math.PI / 180 / 3600);
		slit.setRaHour((int) Math.floor(slit.getWcsX() / 3600 / 15 / Math.cos(decrad)));
		slit.setRaMin((int) Math.floor(60 * 
				((slit.getWcsX() / 3600 / 15 / Math.cos(decrad)) 
						- slit.getRaHour())));
		slit.setRaSec(60 * (((slit.getWcsX() / 3600 / 15 / Math.cos(decrad)) 
				- slit.getRaHour()) * 60 - slit.getRaMin()));
	}

    // Convert the coordinates of a Slit from Ra/Dec to x/y.
    private static void slitObjxyToRaDec(Slit slit, RaDec cp) {
    	slit.setObjDecDeg((int) 
    			(Math.abs(slit.getSlitObjWcsY()) / slit.getSlitObjWcsY()) * 
				(int) Math.floor(Math.abs(slit.getSlitObjWcsY() / 3600)));
		slit.setObjDecMin((int) 
				Math.floor((Math.abs(slit.getSlitObjWcsY() / 3600) - 
				Math.floor(Math.abs(slit.getSlitObjWcsY() / 3600))) * 60));
		slit.setObjDecSec(60 * ((Math.abs(slit.getSlitObjWcsY() / 3600) - 
				Math.floor(Math.abs(slit.getSlitObjWcsY() / 3600))) * 
				60 - slit.getObjDecMin()));
		double decrad = (cp.getYCoordinate() * Math.PI / 180 / 3600);
		slit.setObjRaHour((int) Math.floor
				(slit.getSlitObjWcsX() / 3600 / 15 / Math.cos(decrad)));
		slit.setObjRaMin((int) Math.floor
				(60 * ((slit.getSlitObjWcsX() / 3600 / 15 / Math.cos(decrad)) -
						slit.getObjRaHour())));
		slit.setObjRaSec(60 * 
				(((slit.getSlitObjWcsX() / 3600 / 15 / Math.cos(decrad)) 
						- slit.getObjRaHour()) * 60 - slit.getObjRaMin()));
	}

    // Write the outputRegFile from a slitArray.
    private static void printRegionsFromSlitArray(Slit[] slitArray, RaDec cp, 
    		CSUParameters csuData, double positionAngle, String outputRegFile) {
        RaDec oldRedBoxCoord = 
        	new RaDec((cp.getXCoordinate() + csuData.getXCenter() * 60), 
        			(cp.getYCoordinate()));
        double redBoxCSUx = oldRedBoxCoord.xCoordinate - cp.getXCoordinate();
        double redBoxCSUy = oldRedBoxCoord.yCoordinate - cp.getYCoordinate();
    	double theta = -positionAngle * Math.PI / 180;
    	double redBoxCSUxRotated = redBoxCSUx * Math.cos(theta) 
    			- redBoxCSUy * Math.sin(theta);
    	double redBoxCSUyRotated = redBoxCSUx * Math.sin(theta) 
    			+ redBoxCSUy * Math.cos(theta);
    	double redBoxWCSxRotated = redBoxCSUxRotated + cp.getXCoordinate();
        double redBoxWCSyRotated = redBoxCSUyRotated + cp.getYCoordinate();
    	FileOutputStream out;
    	PrintStream p;
    	try
    	{
    		out = new FileOutputStream(outputRegFile);
    		p = new PrintStream( out );
    		java.text.DecimalFormat wholeNum = new java.text.DecimalFormat("0");
    		java.text.DecimalFormat secondPlace = 
    			new java.text.DecimalFormat("0.00");
    		java.text.DecimalFormat thirdPlace = 
    			new java.text.DecimalFormat("0.000");
    		java.text.DecimalFormat fifthPlace = 
    			new java.text.DecimalFormat("0.00000");
    		p.println("global color=green font=\"helvetica 10 normal\" " +
    				"select=1 highlite=0 edit=0 move=0 delete=1 " +
    				"include=1 fixed=0 source \nfk5");
            p.println("circle(" + 
        			wholeNum.format(cp.getRaHour()) + 
        			":" + wholeNum.format(cp.getRaMin()) +
        			":" + secondPlace.format(cp.getRaSec()) +
        			"," + wholeNum.format(cp.getDecDeg()) +
        			":" + wholeNum.format(cp.getDecMin()) +
        			":" + secondPlace.format(cp.getDecSec()) + "," + 
        			CSU_FP_RADIUS + "\""
        			+ ")\t# color=yellow font=\"helvetica 18 normal\" " +
        			"text={Keck Focal Plane}");
            p.println("box(" 
            		+ fifthPlace.format((cp.getXCoordinate() / 
            				Math.cos(cp.getYCoordinate() *
            						Math.PI / 180 / 3600) / 3600))
            		+ "," + fifthPlace.format(cp.getYCoordinate() / 3600)
            		+ "," +
        			fifthPlace.format(CSU_WIDTH) + "\"" + "," + 
        			fifthPlace.format(CSU_HEIGHT) + "\"" 
        			+ "," + thirdPlace.format(positionAngle) + 
        			")\t# color=magenta font=\"helvetica 18 normal\" " +
        			"text={CSU Plane}");
            p.println("box(" + 
            		fifthPlace.format(redBoxWCSxRotated / 
            				Math.cos(cp.getYCoordinate() * 
            						Math.PI / 180 / 3600) / 3600) +
        			"," + fifthPlace.format(redBoxWCSyRotated / 3600) + 
        			"," + fifthPlace.format((csuData.getXRange() * 60)) + "\""
        			+ "," + fifthPlace.format(CSU_HEIGHT) + "\"" +
        			"," + thirdPlace.format(positionAngle) + 
        			")\t# color=red");
            for (int i = 0; i < slitArray.length; i++) {
            	p.println("circle(" + 
            			fifthPlace.format((slitArray[i].getSlitObjWcsX() / 
            					Math.cos(cp.getYCoordinate() * 
            							Math.PI / 180 / 3600) / 3600)) + 
            			"," + fifthPlace.format(
            					(slitArray[i].getSlitObjWcsY() / 3600)) +
            			",0.5\")\t# " +
            			"text={" + slitArray[i].getSlitObjName() + 
            			"}");
            	p.println("box(" + 
            			fifthPlace.format((slitArray[i].getWcsX() / 
            					Math.cos(cp.getYCoordinate() * 
            							Math.PI / 180 / 3600) / 3600)) + 
            			"," + fifthPlace.format(
            					(slitArray[i].getWcsY() / 3600)) +
            			"," + fifthPlace.format(csuData.getSlitWidth()) + "\""
            			+ "," + fifthPlace.format(slitArray[i].getSlitLength())
            			+ "\"," + thirdPlace.format((BAR_TILT + positionAngle)) 
            			+ ")\t# text={Slit# " + slitArray[i].getSlitNumber() + 
            			"}");
            } 
    		p.close();
    	}
    	catch (Exception e)
    	{
    		System.err.println ("Error writing to file");
    	}
    	System.out.println("The corresponding SAOImage Ds9 region file is: " 
    			+ outputRegFile);
    }
    
    private static void writeOutSlitList(String outputSlitFile, 
    		Slit[] slitArray3, double slitWidth) {
		java.text.DecimalFormat secondPlace = 
			new java.text.DecimalFormat("0.00");
		FileOutputStream out;
		PrintStream p;
		try
		{
			out = new FileOutputStream(outputSlitFile);
			p = new PrintStream(out);
			
	    	java.text.DecimalFormat wholeNum = new java.text.DecimalFormat("0");
	    	
	    	for (int i = 0; i < slitArray3.length; i++) {
	    			p.print(slitArray3[slitArray3.length - 1 - i].
	    											getSlitNumber() + "\t" + 
	    					wholeNum.format(slitArray3[slitArray3.
	    					                           length - 1 - i].
	    					                           getRaHour()) + "\t" + 
	    					wholeNum.format(slitArray3[slitArray3.
	    					                           length - 1 - i].
	    					                           getRaMin()) + "\t" + 
	    					secondPlace.format(slitArray3[slitArray3.
	    					                              length - 1 - i].
	    					                              getRaSec()) + "\t" + 
	    					wholeNum.format(slitArray3[slitArray3.
	    					                           length - 1 - i].
	    					                           getDecDeg()) + "\t" + 
	    					wholeNum.format(slitArray3[slitArray3.
	    					                           length - 1 - i].
	    					                           getDecMin()) + "\t" + 
	    					secondPlace.format(slitArray3[slitArray3.
	    					                              length - 1 - i].
	    					                              getDecSec()) + "\t" + 
	    					secondPlace.format(slitWidth) + "\t" + 
	    					secondPlace.format(slitArray3[slitArray3.
	    					                              length - 1 - i].
	    					                              getSlitLength()) 
	    					                              + "\t" + 
	    					slitArray3[slitArray3.length - 1 - i].
	    					getSlitObjName() + "\t" + 
	    					secondPlace.format(slitArray3[slitArray3.
	    					                              length - 1 - i].
	    					                              getSlitObjPriority())
	    					                              + "\t" + 
	    					secondPlace.format(slitArray3[slitArray3.
	    					                              length - 1 - i].
	    					                              getTargetLocation()) 
	    					                              + "\t" +
	    					wholeNum.format(slitArray3[slitArray3.
	    					                           length - 1 - i].
	    					                           getObjRaHour()) + "\t" + 
	    					wholeNum.format(slitArray3[slitArray3.
	    					                           length - 1 - i].
	    					                           getObjRaMin()) + "\t" + 
	    					secondPlace.format(slitArray3[slitArray3.
	    					                              length - 1 - i].
	    					                              getObjRaSec()) 
	    					                              + "\t" + 
	    					wholeNum.format(slitArray3[slitArray3.
	    					                           length - 1 - i].
	    					                           getObjDecDeg()) + "\t" + 
	    					wholeNum.format(slitArray3[slitArray3.
	    					                           length - 1 - i].
	    					                           getObjDecMin()) + "\t" + 
	    					secondPlace.format(slitArray3[slitArray3.
	    					                              length - 1 - i].
	    					                              getObjDecSec()));
	    			p.print("\n");
	    	}
			p.close();
		}
		catch (Exception e)
		{
			System.err.println ("Error writing to file");
		}
    }
    
    /** Inner Classes **/
    // Comparator for finding the highest priority object. Also factors in 
    // proximity to CSU center X-position when two or more objects have the same
    // maximum priority.
    private class PriorityAndXOrderCompare implements Comparator<AstroObj> {
    	double xCenter;
    	public PriorityAndXOrderCompare(double xCenter) {
    		this.xCenter = xCenter;
    	}
    	public int compare(AstroObj obj1, AstroObj obj2) {
    		if (obj1.getObjPriority() < obj2.getObjPriority())
    			return 1;
    		if (obj1.getObjPriority() > obj2.getObjPriority())
    			return -1;
    		if (obj1.getObjPriority() == obj2.getObjPriority() && 
    				(Math.abs((obj1.getObjX() - CSU_WIDTH / 2) -
    						xCenter * 60) < 
    				Math.abs((obj2.getObjX() - CSU_WIDTH / 2)) - 
    				xCenter * 60))
    			return -1;
    		if (obj1.getObjPriority() == obj2.getObjPriority() && 
    				(Math.abs((obj1.getObjX() - CSU_WIDTH / 2) - 
    						xCenter * 60) > 
    				Math.abs((obj2.getObjX() - CSU_WIDTH / 2)) - 
    					xCenter * 60))
    			return 1;
    		else
    			return 0;
    	}
    };
	
    /** MAIN **/
	public static void main(String[] args) { 
		// Start the timer.
		Timer howLong = new Timer();
		howLong.start();
		
		/** Process the command line input arguments and print the help text
		 * if needed. **/
		// Print the help text if there is an invalid number of arguments or if
		// the first argument is "help".
		String help = "help";
		if ((args.length != 21) && (args.length != 15)) 
		{
			PrintHelpText();
		}
		if (args[0].equals(help)) {
			PrintHelpText();
		}
		
		/** Read in file and create astroObjArrayOriginal. **/
		astroObjArrayOriginal = readInFile(args[0]);
        
        int highObjRaHour = (int) astroObjArrayOriginal[0].getRaHour();
        int lowObjRaHour = (int) astroObjArrayOriginal[0].getRaHour();
        int highObjDecDeg = (int) astroObjArrayOriginal[0].getDecDeg();
        int lowObjDecDeg = (int) astroObjArrayOriginal[0].getDecDeg();
        for (int i = 0; i < astroObjArrayOriginal.length; i++) {
        	if ((int) astroObjArrayOriginal[i].getRaHour() > highObjRaHour) {
        		highObjRaHour = (int) (astroObjArrayOriginal[i].getRaHour());
        	}
        	if ((int) astroObjArrayOriginal[i].getRaHour() < lowObjRaHour) {
        		lowObjRaHour = (int) (astroObjArrayOriginal[i].getRaHour());
        	}
           	if ((int) astroObjArrayOriginal[i].getDecDeg() > highObjDecDeg) {
           		highObjDecDeg = (int) (astroObjArrayOriginal[i].getDecDeg());
        	}
        	if ((int) astroObjArrayOriginal[i].getDecDeg() < lowObjDecDeg) {
        		lowObjDecDeg = (int) (astroObjArrayOriginal[i].getDecDeg());
        	}
        	
        }
        
        
        if (highObjRaHour - lowObjRaHour > 1) {
        	if ((highObjRaHour == 23) && (lowObjRaHour == 0)) {
        	}
        	else {
        		System.err.println("Dummy! You can't expect me to run on an input " +
            			"object list that spans \nmore than one hour in Right " +
            			"Ascension. Exiting . . . ");
            	System.exit(1);
        	}
        }
        if (highObjDecDeg - lowObjDecDeg > 1) {
        	System.err.println("Dummy! You can't expect me to run on an input " +
        			"object list that spans \nmore than one degree in  " +
        			"Declination. Exiting . . . ");
        	System.exit(1);
        }
        boolean raCoordWrap = false;
        if ((highObjRaHour == 23) && (lowObjRaHour == 0)) {
            for (int i = 0; i < astroObjArrayOriginal.length; i++) {
            	if (astroObjArrayOriginal[i].getRaHour() == 0) {
            		astroObjArrayOriginal[i].setRaHour(12);
            		raCoordWrap = true;
            	}
            	if (astroObjArrayOriginal[i].getRaHour() == 23) {
            		astroObjArrayOriginal[i].setRaHour(11);
            		raCoordWrap = true;
            	}
            }
        }
        
        
		// Instnatiate a new RaDec variable to store the input field center.
		RaDec fieldCenter = new RaDec();
		
		// Construct CSUFieldData from the input args array.
		CSUFieldData csuFieldData = 
			new CSUFieldData(args, astroObjArrayOriginal);
		
		fieldCenter = csuFieldData.fieldCenter;
		
		if (raCoordWrap) {
			if (csuFieldData.fieldCenter.raHour == 0) {
				csuFieldData.fieldCenter.raHour = 12;
			}
			if (csuFieldData.fieldCenter.raHour == 23) {
				csuFieldData.fieldCenter.raHour = 11;
			}
		}
		
		
		// Instnatiate a new CSUParameters variable and give it the inputs.
		// Verify that the inputs are permitted.
        CSUParameters verifiedCSUData = new CSUParameters();
        try {
			verifiedCSUData = verifyCSUData(csuFieldData.xRange, 
					csuFieldData.xCenter, csuFieldData.slitWidth,
					csuFieldData.ditherSpace);
		} catch (InvalidArgumentException e) {
			System.out.println(e.getMessage());
			e.printStackTrace();
		}

		// Compute the wcs x and y coordinates of the field center from Ra/Dec.
		raDecToXY(fieldCenter);
        // Instantiate a new Mascgen_beta.
		Mascgen_Delta myMascgen_Delta = new Mascgen_Delta();
		double totalPriority = 0;
		int runNum = 0; // Keep track of the number of optimization runs.
		RaDec tempFieldCenter = new RaDec();
		double tempPA;
		// Convert the Ra/Dec coordinates of the input object list into wcs x 
		// and y in ar seconds.
        for (int i = 0; i < astroObjArrayOriginal.length; i++) {
        	astroObjRaDecToXY(astroObjArrayOriginal[i], fieldCenter);
        }
        
        /** Create the correct RowRegions and OverlapRegions for the 
         * user-supplied CSU Field Data. **/
        // Divide the CSU_HEIGHT arc seconds of CSU height into 46 RowRegions, 
        // 45 OverlapRegions, and 2 "half" OverlapRegions on the edges.
   		orArray[0] = new Region(0, verifiedCSUData.getDeadSpace());
        for (int i = 0; i < 46; i++) {
        	if (i > 0) {
        		orArray[i] = new Region(verifiedCSUData.getDeadSpace() + 
        				i * verifiedCSUData.getRowRegionHeight()
        				+ (i - 1) * 2 * verifiedCSUData.getDeadSpace(), 
        				verifiedCSUData.getDeadSpace() + i * 
        				verifiedCSUData.getRowRegionHeight() + 
        				(i - 1) * 2 * verifiedCSUData.getDeadSpace() + 
        				2 * verifiedCSUData.getDeadSpace());
        	}
        	rrArray[i] = new Region((verifiedCSUData.getDeadSpace() + 
        			2 * verifiedCSUData.getDeadSpace() * i) + 
        			i * verifiedCSUData.getRowRegionHeight(), 
        			(verifiedCSUData.getDeadSpace() + 
        					2 * verifiedCSUData.getDeadSpace() * i) + 
        			(1 + i) * verifiedCSUData.getRowRegionHeight());
        }
        orArray[46] = new Region(CSU_HEIGHT - verifiedCSUData.getDeadSpace(), 
        		CSU_HEIGHT);
        
        // Now, run the three-level for loop over position angle, field center
        // y coordinate, and field center x coordinate. Count the total number
        // of loops (runNum).
		for (int j = -csuFieldData.stepsX; j < csuFieldData.stepsX + 1; j++) {
			tempFieldCenter.setXCoordinate(fieldCenter.getXCoordinate() - 
					j * csuFieldData.stepSizeX);
			for (int k = -csuFieldData.stepsY; k < csuFieldData.stepsY + 1; k++){
				tempFieldCenter.setYCoordinate(fieldCenter.getYCoordinate() -
						k * csuFieldData.stepSizeY);
		        for (int i = 0; i < astroObjArrayOriginal.length; i++) {
		        	astroObjRaDecToXY(astroObjArrayOriginal[i], tempFieldCenter);
		        }
				for (int m = -csuFieldData.paSteps; m < csuFieldData.paSteps + 1; m++) {
					double tempTotalPriority;
					tempPA = csuFieldData.positionAngle + m * csuFieldData.paStepSize;
					AstroObj[] tempAOArray = myMascgen_Delta.optimize(
							astroObjArrayOriginal, 
							verifiedCSUData, tempFieldCenter, tempPA);
					tempTotalPriority = prioritySum(tempAOArray);
					runNum++;
					if (tempTotalPriority > totalPriority) {
						finalFieldCenter.setXCoordinate(
								tempFieldCenter.getXCoordinate());
						finalFieldCenter.setYCoordinate(
								tempFieldCenter.getYCoordinate());
						finalPA = tempPA;
						totalPriority = tempTotalPriority;
						optimumAstroObjArray = tempAOArray;
					}
				}	
			}
		}
        xyToRaDec(finalFieldCenter);
        
        for (int i = 0; i < optimumAstroObjArray.length; i++) {
        	astroObjRaDecToXY(optimumAstroObjArray[i], finalFieldCenter);
        }

		// Make a new array of Slits.
		Slit[] slitArray3;
		// Run the slitConfigurationGenerator on the optimized AstroObj array,
		// which was returned by the above optimize call.
		slitArray3 = myMascgen_Delta.slitConfigurationGenerator(
				optimumAstroObjArray,
				verifiedCSUData.getSlitWidth(), verifiedCSUData.getDeadSpace(), 
				verifiedCSUData.getRowRegionHeight(),
				finalFieldCenter, finalPA, csuFieldData.barPositionList);
		
        if (raCoordWrap) {
            for (int i = 0; i < slitArray3.length; i++) {
            	if (slitArray3[i].getRaHour() == 12) {
            		slitArray3[i].setRaHour(0);
            	}
            	if (slitArray3[i].getObjRaHour() == 12) {
            		slitArray3[i].setObjRaHour(0);
            	}
            	if (slitArray3[i].getRaHour() == 11) {
            		slitArray3[i].setRaHour(23);
            	}
            	if (slitArray3[i].getObjRaHour() == 11) {
            		slitArray3[i].setObjRaHour(23);
            	}
            }
            if (finalFieldCenter.getRaHour() == 12) {
        		finalFieldCenter.setRaHour(0);
        	}
        	if (finalFieldCenter.getRaHour() == 11) {
        		finalFieldCenter.setRaHour(23);
        	}
        } 
        
	    /** The final step is to print the slit configuration. **/
		java.text.DecimalFormat secondPlace = 
			new java.text.DecimalFormat("0.00");
		System.out.println("The optimized slit configuration has been " +
				"found after " + runNum + " runs." +
				"\n\tTotal Priority = " + totalPriority + 
				"\n\tNumber of Slits = " + slitArray3.length +
				"\n\tCSU Center Position = ");
		finalFieldCenter.printRaDecCoordTabs();
		System.out.println("\tPosition Angle = " + 
				secondPlace.format(finalPA) + "¡\n");
		// Write the slit list out.
		writeOutSlitList(csuFieldData.outputSlitFile, slitArray3, 
				verifiedCSUData.getSlitWidth());
		System.out.println("The slit list file is: " + csuFieldData.outputSlitFile);
		// Write out the region file.
		printRegionsFromSlitArray(slitArray3, finalFieldCenter, verifiedCSUData, 
				finalPA, csuFieldData.outputRegFile);
		System.out.println("The corresponding bar position list is: " 
				+ csuFieldData.barPositionList);
		// Report the program execution time.
		howLong.end();
		howLong.printTime();
	}// Close main
} // Close class.