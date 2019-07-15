@echo off

GOTO HEADER

This is a batch file used to pre-process a LAS file for use with the graph-cut
ITC method

Call: "gc_prepoc.bat data_folder las_file output_folder"

The script takes the data directory as its first argument, the name of the raw
LAS file as its second argument, and the location to save the results as its 
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

Example: "gc_prepoc.bat data P4_buffer.las analysis"

Jonathan Williams
jonvw28@gmail.com
12/06/2017

:HEADER

::This section of the code code uses LAStools to pull out ground points and scale the LiDAR returns down to the ground

::Back up raw data
@echo Backing up Raw Data...
IF NOT EXIST %3 mkdir %3
IF NOT EXIST %3\1_raw mkdir %3\1_raw
IF NOT EXIST %3\1_raw\%2 cp %1\%2 %3\1_raw\%2

::Generate header file as text
IF NOT EXIST %3\1_raw\info mkdir %3\1_raw\info
lasinfo -i %3\1_raw\%2 -odir %3\1_raw\info -otxt
@echo Done!
@echo ~~~~~~~~~~~~~~~~~~~~~~

::Filter for noise
@echo Filtering out noise...(expect warning about dropped points)
IF NOT EXIST %3\2_filtered mkdir %3\2_filtered
IF NOT EXIST %3\2_filtered\info mkdir %3\2_filtered\info
lasnoise -i %3\1_raw\%2 -drop_class 7 -remove_noise -odir %3\2_filtered -odix _filtered -olas 
lasinfo -i %3\2_filtered\%~n2_filtered.las -odir %3\2_filtered\info -otxt 
@echo Done!
@echo ~~~~~~~~~~~~~~~~~~~~~~

::rescale points to subtract off ground
@echo Subtracting Ground...
IF NOT EXIST %3\3_levelled mkdir %3\3_levelled
IF NOT EXIST %3\3_levelled\info mkdir %3\3_levelled\info
lasheight -i %3\2_filtered\%~n2_filtered.las -odir %3\3_levelled\ -odix _level -replace_z -olas 
lasinfo -i %3\3_levelled\%~n2_filtered_level.las -odir %3\3_levelled\info -otxt 
@echo Done!
@echo ~~~~~~~~~~~~~~~~~~~~~~