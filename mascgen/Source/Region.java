

public class Region {
	private double minY;
	private double maxY;
	
	/** Constructor to create a OverlapRegion with given values. **/
    public Region(double min, double max) {
        this.minY = min;
        this.maxY = max;
    }
    
    public Region() {
    	this(0, 0);
    }
    
    /** Return the minimum y value of the OverlapRegion. **/
    public double getMinY() {
    	return minY;
    }
    
    /** Return the maximum y value of the OverlapRegion. **/
    public double getMaxY() {
    	return maxY;
    }
    
    /** Return the height of the OverlapRegion. **/
    public double getHeight() {
    	return maxY - minY;
    }
}
