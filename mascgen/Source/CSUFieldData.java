
public class CSUFieldData {
//	 Instnatiate a new RaDec variable to store the input field center.
	RaDec fieldCenter = new RaDec();
	double xRange;
	double xCenter;
	double slitWidth;
	double ditherSpace;
	int stepsX;
	double stepSizeX;
	int stepsY;
	double stepSizeY;
	double positionAngle;
	int paSteps;
	double paStepSize;
	String outputSlitFile;
	String outputRegFile;
	String barPositionList;
	
	
	CSUFieldData(String[] args, AstroObj[] astroObjArrayOriginal) {
		this.xRange = Double.parseDouble(args[1]); // in arc min
	    this.xCenter = -Double.parseDouble(args[2]); // in arc min
	    this.slitWidth = Double.parseDouble(args[3]); // in arc sec
	    this.ditherSpace = Double.parseDouble(args[4]); // in arc sec
		if (args.length == 21) {
			this.fieldCenter.setRaHour(Double.parseDouble(args[5]));
			this.fieldCenter.setRaMin(Double.parseDouble(args[6]));
			this.fieldCenter.setRaSec(Double.parseDouble(args[7]));
			this.fieldCenter.setDecDeg(Double.parseDouble(args[8]));
			this.fieldCenter.setDecMin(Double.parseDouble(args[9]));
			this.fieldCenter.setDecSec(Double.parseDouble(args[10]));
			this.stepsX = Integer.parseInt(args[11]);
			this.stepSizeX = Double.parseDouble(args[12]); // in arc sec
			this.stepsY = Integer.parseInt(args[13]);
			this.stepSizeY = Double.parseDouble(args[14]); // in arc sec
			this.positionAngle = Double.parseDouble(args[15]); // in degrees
			this.paSteps = Integer.parseInt(args[16]);
			this.paStepSize = Double.parseDouble(args[17]); // in degrees
			this.outputSlitFile = args[18];
			this.outputRegFile = args[19];
			this.barPositionList = args[20];
		}
		else {
			this.stepsX = Integer.parseInt(args[5]);
			this.stepSizeX = Double.parseDouble(args[6]); // in arc sec
			this.stepsY = Integer.parseInt(args[7]);
			this.stepSizeY = Double.parseDouble(args[8]); // in arc sec
			this.positionAngle = Double.parseDouble(args[9]); // in degrees
			this.paSteps = Integer.parseInt(args[10]);
			this.paStepSize = Double.parseDouble(args[11]); // in degrees
			this.outputSlitFile = args[12];
			this.outputRegFile = args[13];
			this.barPositionList = args[14];
			double centerRAsecs = raWeightedSum(astroObjArrayOriginal) / 
					prioritySum(astroObjArrayOriginal);
			double centerDECsecs = decWeightedSum(astroObjArrayOriginal) / 
					prioritySum(astroObjArrayOriginal);
			this.fieldCenter.setRaHour(Math.floor(centerRAsecs / 3600));
			this.fieldCenter.setRaMin(Math.floor((centerRAsecs - 3600 * 
					(Math.floor(centerRAsecs / 3600))) / 60));
			this.fieldCenter.setRaSec(centerRAsecs - 
					60 * Math.floor((centerRAsecs - 3600 * 
							(Math.floor(centerRAsecs / 3600))) / 60)
					- 3600 * Math.floor(centerRAsecs / 3600));
			this.fieldCenter.setDecDeg(Math.floor(centerDECsecs / 3600));
			this.fieldCenter.setDecMin(Math.floor((centerDECsecs - 3600 * 
					(Math.floor(centerDECsecs / 3600))) / 60));
			this.fieldCenter.setDecSec(centerDECsecs - 
					60 * Math.floor((centerDECsecs - 3600 * 
							(Math.floor(centerDECsecs / 3600))) / 60)
					- 3600 * Math.floor(centerDECsecs / 3600));
		}
	}
	
	private double decWeightedSum(AstroObj[] array) {
		double result = 0;
		for (int i = 0; i < array.length; i++)
			result += (array[i].getDecDeg() * 3600 + 
					array[i].getDecMin() * 60 +
					array[i].getDecSec()) * array[i].getObjPriority();
		return result;
	}

	private double prioritySum(AstroObj[] array) {
		double result = 0;
		for (int i = 0; i < array.length; i++)
			result += array[i].getObjPriority();
		return result;
	}

	private double raWeightedSum(AstroObj[] array) {
		double result = 0;
		for (int i = 0; i < array.length; i++)
			result += (array[i].getRaHour() * 3600 + 
					array[i].getRaMin() * 60 +
					array[i].getRaSec()) * array[i].getObjPriority();
		return result;
	}

}
