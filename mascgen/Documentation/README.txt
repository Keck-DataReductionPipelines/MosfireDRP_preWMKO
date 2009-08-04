MOSFIRE Automatic Slit Mask Configuration Generator 
Author: ChristopherKlein 
Advisor: Charles Steidel
Summer 2007 
Contact: cklein@caltech.edu

Please read the user's manual, MascgenPaper.pdf, before attempting to
use Mascgen productively.

There are three versions of the program bundled into three different JAR
files:

Mascgen_Delta.jar 
	Command Line version. Must be run in Terminal with the command 
		$ java -jar Mascgen_Delta.jar <arguments> 
	Definitely read the help text (put "help" as the only argument or run 
	with no arguments) before using this version.

Mascgen_Delta_GUI.jar
	GUI version tailored to Mac OS X. Just double click.

Mascgen_Delta_GUI_WinXP.jar
	GUI version tailored to Windows XP. Should be able to run with double 
	click in Windows XP.

If you want to run Mascgen under Linux or Solaris, use the WinXP version. It
is more approriately formatted for those default GUI themes. If the main 
Mascgen window looks a little squished, just resize it. Also, you may have 
to start Mascgen with the Terminal command
	$ java -jar Mascgen_Delta_GUI_WinXP.jar
when running it under Linux or Solaris.

Both GUI versions will run on all four operating systems, but the GUI
formatting will be a little off if you try to run the wrong version
under the wrong OS. The only difference in performance is that the WinXP
version does not update the status window after each new highest
priority slit configuration is found during the iteration run. (The
reason is that the GUI components had to be modified to work under
Windows and Labels had to be switched out for JTextAreas.) If you're
using the WinXP version, you'll just have to wait for the program to
complete before getting a status update.


If you have any problems or any suggestions/comments, please feel free
to contact the author.

