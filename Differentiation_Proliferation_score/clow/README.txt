
The contens of this zip file should include:

	README.txt
	DWD.zip
	Differentiation Predictor/DifferentiationPredictor.R
	Differentiation Predictor/DifferentiationPredictor_functions.R
	Differentiation Predictor/14unsortedSamples_noReplicates_forPaper.txt
	Differentiation Predictor/GSE16997_series_entrez.txt
	Differentiation Predictor/Jason-122-Mouse.txt
	Differentiation Predictor/UNC337arraydata_imputedCollapsed.txt
	

This collection provides all necessary data and scripts required to calculate
the differentiation score model (using GSE16997), and assign the differentiation score to test cases
(UNC337, MMS and Jason122Mouse) as shown in Figure 4 of Prat et al. (2010).

This code has been validated in R version 2.10.1.  The current impute library is required for
implementation. Please see www.r-project.org and www.bioconductor.org for installion of R and impute,
respectively.

The next requirement is installation of the DWD executable. This file has only been compiled for
MS Windows operating systems. The DWD.zip file should be extracted into a directory labeled 'DWD'.  
This directory must be placed in the root 'C' patition such that the executable 'BatchAdjustSM.exe'
exists as C:/DWD/lib/BatchAdjustSM.exe.  The directories C:/DWD/DWDdata and C:/DWD/bin should also
exist if all has been correctly placed.

The step before running the code requires editing the file DifferentiationPredictor.R.  
The line that appears as

setwd("E:/somedir/DifferentiationPredictor/")

should be edited to indicate the location of the 'Differentiation Predictor' directory.  
One example would be

setwd("E:/Documents and Settings/APrat/Desktop/DifferentiationPredictor/")

Save this file and open the R software package.  In R click File->Source R code and select the
DifferentiationPredictor.R file.  The program will load the processed datasets, calculate the 
differentiation score model from the Lim et al. data, write the model to a file, assign the 
differentiation score to each data set, and summarize these results in a pdf figure.




