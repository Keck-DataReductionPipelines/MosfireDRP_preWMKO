


public class RaDec {
	
	double raHour;
	double raMin;
	double raSec;
	
	double decDeg;
	double decMin;
	double decSec;
	
	double epoch;
	double equinox;
	
	double xCoordinate;
	double yCoordinate;
	
	public RaDec(double raHour, double raMin, double raSec, double decDeg, 
			double decMin, double decSec) {
		this.raHour = raHour;
		this.raMin = raMin;
		this.raSec = raSec;
		this.decDeg = decDeg;
		this.decMin = decMin;
		this.decSec = decSec;
	}
	
	public RaDec(double x, double y) {
		this.xCoordinate = x;
		this.yCoordinate = y;
	}
	
	public RaDec() {
		
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
	
	public void printRaDecCoord() {
		System.out.println("RA: " + raHour + ", " + raMin + ", " + raSec);
		System.out.println("DEC: " + decDeg + ", " + decMin + ", " + decSec);
	}
	
	public void printRaDecCoordTabs() {
		java.text.DecimalFormat secondPlace = new java.text.DecimalFormat("0.00");
		System.out.println("\t\tRA: " + raHour + "h, " + raMin + "m, " + secondPlace.format(raSec) + "s");
		System.out.println("\t\tDEC: " + decDeg + "¡, " + decMin + "', " + secondPlace.format(decSec) + "\"");
	}
	
	public void setXCoordinate(double num) {
		xCoordinate = num;
	}
	
	public double getXCoordinate() {
		return xCoordinate;
	}
	
	public void setYCoordinate(double num) {
		yCoordinate = num;
	}
	
	public double getYCoordinate() {
		return yCoordinate;
	}
	
}
