# Robert Lasenby 2009

Mascgen_Delta_GUI_WinXP.class : Mascgen_Delta_GUI_WinXP.java
	javac $<
# see http://java.sun.com/docs/books/tutorial/deployment/jar/appman.html
# for details of jar stuff
Mascgen_Delta_GUI_WinXP.jar : Mascgen_Delta_GUI_WinXP.class
	jar cfe $@ Mascgen_Delta_GUI_WinXP *.class

Mascgen_Delta.class : Mascgen_Delta.java
	javac $<
Mascgen_Delta.jar : Mascgen_Delta.class
	jar cfe $@ Mascgen_Delta *.class

all : Mascgen_Delta_GUI_WinXP.jar

clean : 
	rm *.class
