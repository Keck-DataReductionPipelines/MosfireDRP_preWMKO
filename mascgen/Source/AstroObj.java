
public class AstroObj {

	private String objName; // Object name.
	private double objPriority; // Object priority.
	private double objMag;
	private double objX; // Object x coodinate in arc seconds.
	private double objY; // Object y coordinate in arc seconds.
	private int objRR = -1; // RowRegion occupied by Object.
	private int objOR = -1; // OverlapRegion occupied by Object.
	private double raHour;
	private double raMin;
	private double raSec;
	private double decDeg;
	private double decMin;
	private double decSec;
	private double epoch;
	private double equinox;
	private double junk1;
	private double junk2;
	private double wcsX;
	private double wcsY;
	
	/** Constructor to initialize an AstroObj with given values. **/
    public AstroObj(String name, double priority, double mag, double raHour,
    		double raMin, double raSec, double decDeg, double decMin, 
    		double decSec, double epoch, double equinox, double junk1, 
    		double junk2, int RR, int OR, double objX, double objY) {
        this.objName = name;
        this.objPriority = priority;
        this.objMag = mag;
        this.raHour = raHour;
        this.raMin = raMin;
        this.raSec = raSec;
        this.decDeg = decDeg;
        this.decMin = decMin;
        this.decSec = decSec;
        this.epoch = epoch;
        this.equinox = equinox;
        this.junk1 = junk1;
        this.junk2 = junk2;
        this.objRR = RR;
        this.objOR = OR;
        this.objX = objX;
        this.objY = objY;
    }
    
    public AstroObj(String name, double priority, double mag, double raHour,
    		double raMin, double raSec, double decDeg, double decMin, 
    		double decSec, double epoch, double equinox, double junk1, 
    		double junk2) {
        this.objName = name;
        this.objPriority = priority;
        this.objMag = mag;
        this.raHour = raHour;
        this.raMin = raMin;
        this.raSec = raSec;
        this.decDeg = decDeg;
        this.decMin = decMin;
        this.decSec = decSec;
        this.epoch = epoch;
        this.equinox = equinox;
        this.junk1 = junk1;
        this.junk2 = junk2;
    }
    
    public AstroObj() {
    	this("blank", 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
    }
    
    public void printObj(){
    	System.out.print(objName + "\t\t" + objPriority + "\t\t" + objX + 
    			           "\t\t" + objY + "\t\t" +
    			           	"RR: " + objRR + 
    			           "\t\tOR: " + objOR + 
    			           "\n");
    }
    
    public void printObjRaDec() {
    	System.out.print(objName + "\t\t" + objPriority + "\t\t" + raHour + 
		           "\t\t" + raMin + "\t\t" + raSec + "\t\t" + decDeg + "\t\t" +
		           decMin + "\t\t" + decSec +
		           "\n");
    }
    
    /** Return the name of the AstroObj. **/
    public String getObjName() {
        return objName;
    }
    
    public void setObjName(String name) {
        objName = name;
    }
    
    /** Return the priority of the AstroObj. **/
    public double getObjPriority() {
        return objPriority;
    }
    
    public void setObjPriority(double num) {
        objPriority = num;
    }
    
	public void setObjMag(double num) {
		objMag = num;
	}
	
	public double getObjMag() {
		return objMag;
	}
    
    /** Return the x coordinate of the AstroObj. **/
    public double getObjX() {
        return objX;
    }
    
    /** Set the x coordinate of the AstroObj. **/
    public void setObjX(double num) {
    	this.objX = num;
    }
    
    /** Return the y coordinate of the AstroObj. **/
    public double getObjY() {
        return objY;
    }
    
    /** Set the y coordinate of the AstroObj. **/
    public void setObjY(double num) {
    	this.objY = num;
    }
    
    /** Set the Object's RowRegion number. **/
    public void setObjRR(int num) {
    	this.objRR = num;
    }
    
    /** Return the Object's RowRegion number. **/
    public int getObjRR() {
    	return objRR;
    }
    
    /** Set the Object's OverlapRegion number. **/
    public void setObjOR(int num) {
    	this.objOR = num;
    }
    
    /** Return the Object's OverlapRegion number. **/
    public int getObjOR() {
    	return objOR;
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
    
	public void setEpoch(double num) {
		epoch = num;
	}
	
	public double getEpoch() {
		return epoch;
	}
	
	public void setEquinox(double num) {
		equinox = num;
	}
	
	public double getEquinox() {
		return equinox;
	}
	
	public void setJunk1(double num) {
		junk1 = num;
	}
	
	public double getJunk1() {
		return junk1;
	}
	
	public void setJunk2(double num) {
		junk2 = num;
	}
	
	public double getJunk2() {
		return junk2;
	}
	
    /** Return true if the Object's name is not "blank". **/
    public boolean isNotBlank() {
    	if (objName == "blank")
    		return false;
    	else
    		return true;
    }
    
    /** Return true if the Object's name is "blank". **/
    public boolean isBlank() {
    	if (objName == "blank")
    		return true;
    	else
    		return false;
    }

}
