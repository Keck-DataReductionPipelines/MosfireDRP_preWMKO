import java.awt.Frame;
import java.awt.GridLayout;
import java.awt.Label;
import java.awt.Panel;
import java.awt.TextField;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Scanner;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;



class RegGen extends JFrame {
	private static final long serialVersionUID = 1L;
	String argObjectList = null;
	TextField field1 = new TextField(argObjectList);
	String argRegionFile = "RegionFile.reg";
	TextField field2 = new TextField(argRegionFile);
	String radius = "1";
	TextField field3 = new TextField(radius);
    public RegGen() {
    	    
    	    
        // 1... Create/initialize components
        JButton runBtn = new JButton("Run");
        runBtn.addActionListener(new RunBtnListener());
        
        JPanel panel0 = new JPanel();
        JLabel regGenTitle = new JLabel("RegGen - SAOImage Ds9 Region " +
        		"File Generator");
        panel0.add(regGenTitle);


        // 2... Create content panel, set layout
        JPanel panel1 = new JPanel();
        panel1.setLayout(new GridLayout(3,3,0,5));
        panel1.add(new Label("Input Object List"));
        field1.setText(argObjectList);
        panel1.add(field1);
        panel1.add(new Label(" "));
        
        panel1.add(new Label("Output Region File"));
        field2.setText(argRegionFile);
        panel1.add(field2);
        panel1.add(new Label("end with .reg"));
        
        panel1.add(new Label("Region Circle Radius"));
        field3.setText(radius);
        panel1.add(field3);
        panel1.add(new Label("in arc seconds"));
        
        JPanel panel2 = new JPanel();
        panel2.add(runBtn);

        this.add(panel0,"North");
        this.add(panel1, "Center");
        this.add(panel2, "South");
        this.setSize(420,170);
        this.addWindowListener(new Terminate());
    }
	

	class RunBtnListener implements ActionListener {
        public void actionPerformed(ActionEvent e) {
        	AstroObj[] objectList = readInFile(field1.getText());
        	printRegionsFromObjectArray(objectList, field2.getText(), 
        			field3.getText());
    	}
    }
    
    class Terminate extends WindowAdapter{
  	  public void windowClosing(WindowEvent e){
  	    //terminate the program when the window is closed  
  	    System.exit(0);
  	  }//end windowClosing
  	}//end class Terminate
    
    // Print a region file from an astroObjArray.
    private static void printRegionsFromObjectArray(AstroObj[] astroObjArray, 
    		String regFile, String radius) {
    	
        // Calculate average object priority.
        double cumPri = 0;
        for (int i = 0; i < astroObjArray.length; i++) {
        	cumPri += astroObjArray[i].getObjPriority();
        }
        double avgPri = cumPri / astroObjArray.length;
        
        double sigmaSum = 0;
        for (int i = 0; i < astroObjArray.length; i++) {
        	sigmaSum += (Math.pow(astroObjArray[i].getObjPriority() - avgPri, 2));
        }
        double sigma = Math.sqrt(sigmaSum / astroObjArray.length);
    	
    	FileOutputStream out1;
		PrintStream p1;
		try
		{
			out1 = new FileOutputStream(regFile);
			p1 = new PrintStream(out1);
			
			p1.println("global font=\"helvetica 12 normal\"");
			java.text.DecimalFormat wholeNum = new java.text.DecimalFormat("0");
	        for (int i = 0; i < astroObjArray.length; i++) {
	        	if (astroObjArray[i].getObjPriority() <= avgPri - sigma){
	        		p1.println("fk5;circle(" + 
		        			wholeNum.format(astroObjArray[i].getRaHour()) + 
		        			":" + wholeNum.format(astroObjArray[i].getRaMin()) +
		        			":" + astroObjArray[i].getRaSec() +
		        			"," + wholeNum.format(astroObjArray[i].getDecDeg()) +
		        			":" + wholeNum.format(astroObjArray[i].getDecMin()) +
		        			":" + astroObjArray[i].getDecSec() + "," +
		        			radius + "\")\t# " +
		        			"color=red text={" + astroObjArray[i].getObjName() + 
		        			", p = " + astroObjArray[i].getObjPriority() + "}");
	        	}
	        	if ((astroObjArray[i].getObjPriority() > avgPri - sigma) && 
	        			(astroObjArray[i].getObjPriority() <= avgPri)){
	        		p1.println("fk5;circle(" + 
		        			wholeNum.format(astroObjArray[i].getRaHour()) + 
		        			":" + wholeNum.format(astroObjArray[i].getRaMin()) +
		        			":" + astroObjArray[i].getRaSec() +
		        			"," + wholeNum.format(astroObjArray[i].getDecDeg()) +
		        			":" + wholeNum.format(astroObjArray[i].getDecMin()) +
		        			":" + astroObjArray[i].getDecSec() + "," +
		        			radius + "\")\t# " +
		        			"color=yellow text={" + astroObjArray[i].getObjName() + 
		        			", p = " + astroObjArray[i].getObjPriority() + "}");
	        	}
	        	if ((astroObjArray[i].getObjPriority() > avgPri) && 
	        			(astroObjArray[i].getObjPriority() <= avgPri + sigma)){
	        		p1.println("fk5;circle(" + 
		        			wholeNum.format(astroObjArray[i].getRaHour()) + 
		        			":" + wholeNum.format(astroObjArray[i].getRaMin()) +
		        			":" + astroObjArray[i].getRaSec() +
		        			"," + wholeNum.format(astroObjArray[i].getDecDeg()) +
		        			":" + wholeNum.format(astroObjArray[i].getDecMin()) +
		        			":" + astroObjArray[i].getDecSec() + "," +
		        			radius + "\")\t# " +
		        			"color=magenta text={" + astroObjArray[i].getObjName() + 
		        			", p = " + astroObjArray[i].getObjPriority() + "}");
	        	}
	        	if (astroObjArray[i].getObjPriority() > avgPri + sigma){
	        		p1.println("fk5;circle(" + 
		        			wholeNum.format(astroObjArray[i].getRaHour()) + 
		        			":" + wholeNum.format(astroObjArray[i].getRaMin()) +
		        			":" + astroObjArray[i].getRaSec() +
		        			"," + wholeNum.format(astroObjArray[i].getDecDeg()) +
		        			":" + wholeNum.format(astroObjArray[i].getDecMin()) +
		        			":" + astroObjArray[i].getDecSec() + "," +
		        			radius + "\")\t# " +
		        			"color=blue text={" + astroObjArray[i].getObjName() + 
		        			", p = " + astroObjArray[i].getObjPriority() + "}");
	        	}
	        }

			p1.close();
		}
		catch (Exception er)
		{
			System.err.println ("Error writing to file");
		}
    }
    
    
//  Read in file and return array of AstroObjs.
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
	    	{
				Frame errorFrame = new Frame("Error!");
				Panel errorPanel = new Panel();
				Label errorMsg = new Label("Object List not found. " +
						"\nDisregard the output from this Mascgen execution.");
				errorPanel.add(errorMsg);
				errorFrame.add(errorPanel);
				errorFrame.setSize(440,75);
				errorFrame.setLocation(640,0);
	    		errorFrame.addWindowListener(new WindowAdapter() {
			        public void windowClosing(WindowEvent evt) {
			            Frame frame = (Frame)evt.getSource();
			            frame.setVisible(false);
			            frame.dispose();
			        }
			    });
	    		errorFrame.setVisible(true);
			}
		}   
	    AstroObj[] astroObjArray = 
	        astroObjArrayList.toArray(new AstroObj[astroObjArrayList.size()]);
	    return astroObjArray;
	}
	
	
	public static void main(String[] args) {
        RegGen window = new RegGen();
        window.setVisible(true);
	}

}
