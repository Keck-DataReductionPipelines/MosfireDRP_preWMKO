
public class Slit {

	private int slitNumber; // The number of the slit.
	private double slitWidth; // Slit width.
	private double slitLength; // Slit length.
	private double slitX = -1; // Slit center x-coordinate.
	private double slitY = -1; // Slit center y-coordinate.
	private String slitObjName; // Object name.
	private double slitObjPriority; // Object priority.
	private double objMag;
	private double raHour;
	private double raMin;
	private double raSec;
	private double decDeg;
	private double decMin;
	private double decSec;
	private double epoch;
	private double equinox;
	private double wcsX;
	private double wcsY;
	private double slitObjX; // Object x coodinate in arc seconds.
	private double slitObjY; // Object y coordinate in arc seconds.
	private double slitObjWcsX; // Object x coodinate in arc seconds.
	private double slitObjWcsY; // Object y coordinate in arc seconds.
	private double objRaHour;
	private double objRaMin;
	private double objRaSec;
	private double objDecDeg;
	private double objDecMin;
	private double objDecSec;
	private int slitMul = -1; // Slit length multiple (double, triple, etc)
	private double targetLocation; // Distance in arcsecs vertically (along the 
								   // slit) from center of target to center of
	                               // slit.
	
	
	/** Constructor to initialize a slit with given values. **/
    public Slit(int num, double width, double length, String name, 
    		double priority, double x, double y) {
        this.slitNumber = num;
        this.slitWidth = num;
        this.slitLength = length;
    	this.slitObjName = name;
        this.slitObjPriority = priority;
        this.slitObjX = x;
        this.slitObjY = y;
    }
	
    public Slit() {
    	this(-1, -1, -1, "blank", 0, -1, -1);
    }
    
    /** Print the slit information. **/
    public void printSlit(){
    	java.text.DecimalFormat secondPlace = new java.text.DecimalFormat("0.00");
    	java.text.DecimalFormat wholeNum = new java.text.DecimalFormat("0");
    	
    	System.out.print(slitNumber + "\t" + 
    			wholeNum.format(raHour) + "\t" + 
    			wholeNum.format(raMin) + "\t" + 
    			secondPlace.format(raSec) + "\t" + 
    			wholeNum.format(decDeg) + "\t" + 
    			wholeNum.format(decMin) + "\t" + 
    			secondPlace.format(decSec) + "\t" + 
    			secondPlace.format(slitWidth) + "\t" + 
    			secondPlace.format(slitLength) + "\t" + 
    			slitObjName + "\t" + 
    			secondPlace.format(slitObjPriority) + "\t" + 
    			secondPlace.format(targetLocation) + "\t" +
    			wholeNum.format(objRaHour) + "\t" + 
    			wholeNum.format(objRaMin) + "\t" + 
    			secondPlace.format(objRaSec) + "\t" + 
    			wholeNum.format(objDecDeg) + "\t" + 
    			wholeNum.format(objDecMin) + "\t" + 
    			secondPlace.format(objDecSec));
    	System.out.print("\n");
    }
    
    public void printSlit2(){
    	java.text.DecimalFormat secondPlace = new java.text.DecimalFormat("0.00");
    	java.text.DecimalFormat wholeNum = new java.text.DecimalFormat("0");
    	
    	System.out.print(slitNumber + "\t" + 
    			wholeNum.format(raHour) + "\t" + 
    			wholeNum.format(raMin) + "\t" + 
    			secondPlace.format(raSec) + "\t" + 
    			wholeNum.format(decDeg) + "\t" + 
    			wholeNum.format(decMin) + "\t" + 
    			secondPlace.format(decSec) + "\t" + 
    			secondPlace.format(slitWidth) + "\t" + 
    			secondPlace.format(slitLength) + "\t" + 
    			slitObjName + "\t" + 
    			secondPlace.format(slitObjPriority) + "\t" + 
    			secondPlace.format(targetLocation) + "\t" +
    			wholeNum.format(objRaHour) + "\t" + 
    			wholeNum.format(objRaMin) + "\t" + 
    			secondPlace.format(objRaSec) + "\t" + 
    			wholeNum.format(objDecDeg) + "\t" + 
    			wholeNum.format(objDecMin) + "\t" + 
    			secondPlace.format(objDecSec) + "\t" +
    			slitMul);
    	System.out.print("\n");
    }
    public void printSlitShort(){
    	System.out.print(slitObjX + " " + slitY + " " + 
    			slitWidth + " " + slitLength);
    	System.out.print("\n");
    }
    public void printSlitShort2(){
    	System.out.print(slitObjX + " " + slitY + " " + 
    			slitWidth / 2 + " " + slitLength / 2);
    	System.out.print("\n");
    }
    
    /** Set the Slit x coordinate. **/
    public void setSlitX(double num){
    	this.slitX = num;
    }
    
    /** Return the slit x coordinate. **/
    public double getSlitX() {
    	return slitX;
    }
    
    /** Set the Slit y coordinate. **/
    public void setSlitY(double num){
    	this.slitY = num;
    }
    
    /** Return the slit y coordinate. **/
    public double getSlitY() {
    	return slitY;
    }
    
    /** Set the Slit Number. **/
    public void setSlitNumber(int num){
    	this.slitNumber = num;
    }
    
    /** Return the slit number. **/
    public int getSlitNumber() {
    	return slitNumber;
    }
    
    /** Return the slit object name. **/
    public String getSlitObjName() {
    	return slitObjName;
    }
    
    /** Set the slit object name. **/
    public void setSlitObjName(String name) {
    	this.slitObjName = name;
    }
    
    /** Set the slit object priority. **/
    public void setSlitObjPriority(double num) {
    	this.slitObjPriority = num;
    }
    
    /** Return the slit object priority. **/
    public double getSlitObjPriority() {
    	return slitObjPriority;
    }
    
    /** Set the slit width. **/
    public void setSlitWidth(double num) {
    	this.slitWidth = num;
    }
    
    /** Set the slit object x coordinate. **/
    public void setSlitObjX(double num) {
    	this.slitObjX = num;
    }
    
    /** Return the slit object x coordinate. **/
    public double getSlitObjX() {
    	return slitObjX;
    }
    
    /** Return the slit object y coordinate. **/
    public double getSlitObjY() {
    	return slitObjY;
    }
    
    /** Set the slit object y coordinate. **/
    public void setSlitObjY(double num) {
    	this.slitObjY = num;
    }
    /** Set the slit object x coordinate. **/
    public void setSlitObjWcsX(double num) {
    	this.slitObjWcsX = num;
    }
    
    /** Return the slit object x coordinate. **/
    public double getSlitObjWcsX() {
    	return slitObjWcsX;
    }
    
    /** Return the slit object y coordinate. **/
    public double getSlitObjWcsY() {
    	return slitObjWcsY;
    }
    
    /** Set the slit object y coordinate. **/
    public void setSlitObjWcsY(double num) {
    	this.slitObjWcsY = num;
    }
    /** Set the slit length multiple. **/
    public void setSlitMul(int num) {
    	this.slitMul = num;
    }
    
    /** Return the slit length multiple. **/
    public int getSlitMul() {
    	return slitMul;
    }
    
    /** Set the slit length. **/
    public void setSlitLength(double num) {
    	this.slitLength = num;
    }
    
    public void setRaHour(double num) {
		raHour = num;
	}
	
	public double getRaHour() {
		return raHour;
	}
	
	public void setRaMin(double num) {
		raMin = num;
	}
	
	public double getRaMin() {
		return raMin;
	}
	
	public void setRaSec(double num) {
		raSec = num;
	}
	
	public double getRaSec() {
		return raSec;
	}
	public void setDecDeg(double num) {
		decDeg = num;
	}
	
	public double getDecDeg() {
		return decDeg;
	}
	
	public void setDecMin(double num) {
		decMin = num;
	}
	
	public double getDecMin() {
		return decMin;
	}
	
	public void setDecSec(double num) {
		decSec = num;
	}
	
	public double getDecSec() {
		return decSec;
	}
	
	public void setWcsX(double num) {
		wcsX = num;
	}
	
	public double getWcsX() {
		return wcsX;
	}
	
	public void setWcsY(double num) {
		wcsY = num;
	}
	
	public double getWcsY() {
		return wcsY;
	}
	
	public void setEquinox(double num) {
		equinox = num;
	}
	
	public double getEquinox() {
		return equinox;
	}
	
	public void setEpoch(double num) {
		epoch = num;
	}
	
	public double getEpoch() {
		return epoch;
	}
	
	public void setObjMag(double num) {
		objMag = num;
	}
	
	public double getObjMag() {
		return objMag;
	}
	
    /** Return true if the SlitObject's name is not "blank". **/
    public boolean isNotBlank() {
    	if (slitObjName == "blank")
    		return false;
    	else
    		return true;
    }
    
    /** Return true if the SlitObject's name is "blank". **/
    public boolean isBlank() {
    	if (slitObjName == "blank")
    		return true;
    	else
    		return false;
    }
    
    public double getTargetLocation() {
    	return targetLocation;
    }
    
    public void setTargetLocation(double num) {
    	targetLocation = num;
    }
    
    
    public void setObjRaHour(double num) {
		objRaHour = num;
	}
	
	public double getObjRaHour() {
		return objRaHour;
	}
	
	public void setObjRaMin(double num) {
		objRaMin = num;
	}
	
	public double getObjRaMin() {
		return objRaMin;
	}
	
	public void setObjRaSec(double num) {
		objRaSec = num;
	}
	
	public double getObjRaSec() {
		return objRaSec;
	}
	public void setObjDecDeg(double num) {
		objDecDeg = num;
	}
	
	public double getObjDecDeg() {
		return objDecDeg;
	}
	
	public void setObjDecMin(double num) {
		objDecMin = num;
	}
	
	public double getObjDecMin() {
		return objDecMin;
	}
	
	public void setObjDecSec(double num) {
		objDecSec = num;
	}
	
	public double getObjDecSec() {
		return objDecSec;
	}

    public double getSlitLength() {
    	return slitLength;
    }
    
    
}
