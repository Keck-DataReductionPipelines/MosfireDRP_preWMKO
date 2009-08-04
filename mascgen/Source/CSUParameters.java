

public class CSUParameters {
	private double minLegalX;
	private double maxLegalX;
	private double deadSpace;
	private double rowRegionHeight;
	private double xRange;
	private double xCenter;
	private double slitWidth;
	private double ditherSpace;
	
	public CSUParameters(double minLegalX, double maxLegalX, double deadSpace, 
			double rowRegionHeight, double xRange, double xCenter, double 
			slitWidth, double ditherSpace) {
		this.minLegalX = minLegalX;
		this.maxLegalX = maxLegalX;
		this.deadSpace = deadSpace;
		this.rowRegionHeight = rowRegionHeight;
		this.xRange = xRange;
		this.xCenter = xCenter;
		this.slitWidth = slitWidth;
		this.ditherSpace = ditherSpace;
	}
	
	public CSUParameters() {
		
	}
	
	public double getMinLegalX() {
		return minLegalX;
	}
	
	public double getMaxLegalX() {
		return maxLegalX;
	}
	
	public double getDeadSpace() {
		return deadSpace;
	}
	
	public double getRowRegionHeight() {
		return rowRegionHeight;
	}
	
	public double getXRange() {
		return xRange;
	}
	
	public double getXCenter() {
		return xCenter;
	}
	
	public double getSlitWidth() {
		return slitWidth;
	}
	
	public double getDitherSpace() {
		return ditherSpace;
	}
	
}
