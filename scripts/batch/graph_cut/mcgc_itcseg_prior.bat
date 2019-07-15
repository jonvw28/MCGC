@echo off
GOTO HEADER

This batch file takes a LAS file and creates a prior using the variable window
size maxima finder from ITCsegment methodolgy and saves this out in a format
appropriate for the mcgc pipeline

Call: "mcgc_itcseg_prior.bat data_folder out_folder las_file allom grid_size 
		height_threshold tau"

The script takes the data directory as its first argument, the location to save
the output as its second and the name of the raw LAS file as its third argument.
The fourth argument is then the allometry lookup table suitable for ITCsegment.
It also takes the grid size for preparing the CHM grid, a minimum height
threshold for each prior tree to be included in the prior. The final argument
is which quantile of allometry to use in look up table.

The script will then use the variable window size maxima finder from ITCsegment.
This will be from CHM created from the LAS file. This is then cleaned to a csv 
containing only the x, y and top height columns. The priors are also filtered to
include only 'trees' where their height is above the threshold 

* The cleaned prior goes into the folder specified by out_folder
* Typical values should be 0.5 for grid_size, 10 for height thresholds

Arguments:

* data_folder:		Folder where LAS file is stored
* out_folder:		Folder in which to store the resulting prior
* las_file:			LAS file name
* allom:			filename of .csv lookup of a,b for quantile regression
					of height-radius allometry for R = exp(a) H^{b}
* grid_size:		Width for each pixel in rasterised CHM, in units of LAS file
* height_threshold:	Minimum height in the CHM for a tree to be kept
* tau:				Quantile for use in lookup table of H vs R table

NOTE The directories must be referenced with respect to the directory where
	 this script is called
NOTE data must already have noise filtered, and the returns levelled with 
	 respect to the ground
NOTE This script requires LASTools and R to be added to the system path variable					
					
Example: "mcgc_itcseg_prior.bat data anaylsis data.las lut.csv 0.5 10 50"


Jonathan Williams
jonvw28@gmail.com
05/02/2018

:HEADER

IF EXIST %2\%~n3_prior.csv (
	@echo Error: %2\%~n3_prior.csv output already exists, remove or rename this file and try again
	GOTO ERREND
)

@echo Computing Prior...
IF NOT EXIST %2 mkdir %2
::Get a Canopy Height Model
lascanopy -i %1\%~3 -step %5 -max -o %2\%~n3_CHM.tif
::Generate Prior
Rscript %~dp0..\..\R\itcseg_prior.r %2\%~n3_CHM_max.tif %5 %6 %4 %7

del %2\*.tif
del %2\*.tfw
ren %2\%~n3_CHM_max_prior.csv %~n3_prior.csv 
@echo Done!

:ERREND