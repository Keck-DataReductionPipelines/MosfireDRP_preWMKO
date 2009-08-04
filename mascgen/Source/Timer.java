

public class Timer {
	  private long start;
	  private long end;

	  public Timer() {
	    reset();
	    }

	  public void start() {
	    System.gc();
	    start = System.currentTimeMillis();
	    }

	  public void end() {
	    System.gc();
	    end = System.currentTimeMillis();
	    }

	  public long duration(){
	    return (end-start);
	    }

	  public void reset() {
	    start = 0;  end   = 0;
	    }
	  
	  public void printTime() {
		  System.out.println("\nProgram execution time: " + 
					(end-start) / 1000.00 + " seconds.");
	  }
	  
	}
