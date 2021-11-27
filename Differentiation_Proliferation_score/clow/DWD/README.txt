README

Introduction

Distance Weighted Discrimination (DWD) can perform systematic bias adjustment 
in microarray data (https://genome.unc.edu/pubsup/dwd/). 
The current version is open source and runs in Java environment.  
You can download JRE at http://java.sun.com/j2se/1.4.2/download.html. JRE must 
be 1.4 or above.  

It is absolutely necessary to verify that the default JRE is the one you 
downloaded, rather than something else, such as from Oracle (1.3.1 or older version). 

You can start DOS and type java -version to verify it.


This document will focus on the following areas: Installation and Running of 
the software.
 

Installation

The current version of DWD is only tested in Win XP and 2000.  The software 
should be directly installed under c: drive, or else you would have some 
configuration work to do.  There are two subdirectories under the DWD directory:
lib and DWDdata.  There are also eight important files in the DWD directory: 
README.txt,RunDWD.bat,DWD.jar and other five jar files.

The lib subdirectory contains all libraries and converted c executable. 

The DWDdata subdirectory can house all data files, including your original 
ones and output from DWD.  

Note: You do not have to put your input and output files in this directory, 
but this folder is critical.  If you delete this folder, an error will 
pop up to reminder you this. 


Running of DWD

Two methods:
Method 1: Start DOS command line window, cd to DWD, type RunDWD.bat.
Method 2: In window explorer, double click RunDWD.bat.

I tend to use the method 1, because an error can be retained after you close DWD.


Some Notes Related to Data Files:

1.Data set format  
	The current version of DWD can take Stanford-like text delimited file 
	and MAGE-ML format file as input.
	
   A. For Stanford-like text delimited file.
    	The first column of the Stanford-like text delimited file is identification, 
	the second column contains some annotation, and data start from the 
	third column.  The gene identification must be unique (duplicates are
	not allowed right now). The first row contains the information of samples.  
	This kind of data sets has been extensively tested.  Two sample files 
	are also included in the DWDdata subdirectory. 
	The software can also detect the missing elements and non-digit values 
	where a digit value is supposed to be and report each error one 
	by one each time. 
	
	
	Missing values should first be imputed. see 
	http://bioinformatics.oupjournals.org/cgi/reprint/17/6/520.pdf for a good 
	introduction to the main ideas of, and methods typically used for, imputation, 
	and see http://www.scripps.edu/researchservices/dna_array/new/Data_Analysis_SAM.htm 
	for imputation software. When you edit the file with MS Excel, please delete any 
	blank lines. We suggest using TextPad editor to modify your files.
	Here is the website to download TextPad: http://www.textpad.com/download/
	
   B. For MAGE-ML format file.
    	We have tested files generated from Agilent and Affymatrix, either with internal 
    	data or external data files. 
    	At this time, validation with MAGE-ML.dtd does not work. Validation itself works, 
    	but it massed up with the reading internal data generated with Agilent software. 
    	The work-around is to delete those lines related to MAGE-ML.dtd in the .xml input 
    	files.  This is very important. 
	

2. Parameter selection
   A. For Stanford-like text delimited file.
   	It is straight forward as the GUI displayed.
   	
   B. For MAGE-ML format file.	
	When MAGE-ML format files are used as input, there are a few criteria we adopt.
	
     (1). If there are both raw data (MeasuredBioAssayData) and derived data 
	(DerivedBioAssayData), only derived data will be retrieved for data adjustment.
     (2). If there are more than one raw data/derived data in the one MAGE-ML file, they 
     	must be adjusted first.
     (3). If you want to re-run the data adjustment,you have to restart DWD. The button 
     	"Reset"	in both BioAssayData selection and Quantitation selection will only reset
     	those two parameters. 


3. Output
   A. For Stanford-like text delimited file.
	The output files include: DWD_input.txt, DWD_Vec.txt, and 
	DWD_Non_Std_Output.txt/ DWD_Std_Output.txt (if you use default output).  
	DWD_Non_Std_Output.txt/ DWD_Std_Output.txt  is the final output 
	corresponding to the two DWD types.  But other files (DWD_input.txt and DWD_Vec.txt)
	will also be used in the visual diagnostics analysis.  Please do not delete them.  
	They are automatically overwritten from one run of DWD to another.
	
	
   B. For MAGE-ML format file.
   	In addition to the same files generated as above, there are two extra output files, 
   	when MAGE-ML format files are used as input.  They are ExternalAdjustedDataFile.txt 
   	and DWD_Non_Std_Output.xml/ DWD_Std_Output.xml (if you use default output). The
   	adjusted data will be stored in an external text file named ExternalAdjustedDataFile.txt.
   	The text file, DWD_Non_Std_Output.txt/ DWD_Std_Output.txt (if you use default output)
   	will be used for the visualization.
   	
