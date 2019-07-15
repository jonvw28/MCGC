@echo off

GOTO HEADER

This is a batch file used to pre-process a folder of LAS files for use with the
graph-cut ITC method

Call: "gc_prepoc.bat data_folder output_folder n_cores"

The script takes the data directory as its first argument, the location where 
the results are to be placed as its second, and number of cores to use for its
third

The script will then create and populate various directories as relevant for
the analysis.
* The denoised data goes into \2_filtered
* The ground-subtracted data goes into \3_levelled

Arguments:

* data_folder:		Folder where LAS file is stored
* las_file:			LAS file name
* output_folder:	Folder in which to store the results in their subfolders

NOTE The directories must be referenced with respect to the directory where
	 this script is called
NOTE Raw data must have ground points pre-labelled as class 2 and noise as 
	 class 7
NOTE The raw data will be moved into a sub-folder \data\1_raw so as to preserve 
	 this for the user
NOTE This script requires LASTools to be added to the system path variable

Example: "gc_prepoc.bat data analysis 3"

Jonathan Williams
jonvw28@gmail.com
12/06/2017

:HEADER

::This section of the code code uses LAStools to pull out ground points and scale the LiDAR returns down to the ground

::Back up raw data
@echo Backing up Raw Data...
IF NOT EXIST %2 mkdir %2
IF NOT EXIST %2\1_raw mkdir %2\1_raw
FOR %%f in (%1\*.las) do (
	IF NOT EXIST %2\1_raw\%%~nf.las cp %1\%%~nf.las %2\1_raw\%%~nf.las
)

::Generate header file as text
IF NOT EXIST %2\1_raw\info mkdir %2\1_raw\info
lasinfo -i %2\1_raw\*.las -odir %2\1_raw\info -otxt -cores %3
@echo Done!
@echo ~~~~~~~~~~~~~~~~~~~~~~

::Filter for noise
@echo Filtering out noise...(expect warning about dropped points)
IF NOT EXIST %2\2_filtered mkdir %2\2_filtered
IF NOT EXIST %2\2_filtered\info mkdir %2\2_filtered\info
lasnoise -i %2\1_raw\*.las -drop_class 7 -remove_noise -odir %2\2_filtered -odix _filtered -olas -cores %3
lasinfo -i %2\2_filtered\*.las -odir %2\2_filtered\info -otxt -cores %3
@echo Done!
@echo ~~~~~~~~~~~~~~~~~~~~~~

::rescale points to subtract off ground
@echo Subtracting Ground...
IF NOT EXIST %2\3_levelled mkdir %2\3_levelled
IF NOT EXIST %2\3_levelled\info mkdir %2\3_levelled\info
lasheight -i %2\2_filtered\*.las -odir %2\3_levelled\ -odix _level -replace_z -olas -cores %3
lasinfo -i %2\3_levelled\*.las -odir %2\3_levelled\info -otxt -cores %3
@echo Done!