# Script to generate table of local maxima of CHM in tiff format as provided 
# by LASTools
#
# Syntax:
#
#		>> RScript prior_cleaner.R  data_file scale min_h allom_lut/window_size [tau]
#
# Example:
#
#		>> Rscript itcseg_prior.r data.las 0.5 10 allom_lut.csv 90
#
# Method:
#
#       This script takes the CHM in the tiff file and finds local maxima using
#       the apporach found in ITCsegment package. If an allom lookup table is
#       supplied then an allometrically scaling search window is used, otherwise
#       a fixed window size is used as set by the same argument.
#
#       here a maxima is simply a point which has no point higher than it within 
#       its search window. The output is a csv file with x,y,z cooridnates of 
#       local maxima of the CHM
#
# Inputs:
#
#	data_file:              Name of CHM in TIFF format (with raw height values)
#       scale:                  Cell size of CHM
#       min_h:                  Minimum height of a local maxima to be kept
#       allom_lut/window_size:  If you want an allometrically scaling search
#                               window then put the name of a csv file 
#                               containing a,b for the 100 percentiles of
#                               quantile regression for R = exp(a)*H^b
#                               Otherwise, put the fixed window width in terms
#                               of pixel numbers
#       [tau]:                  optional choice of percentile to set boundary
#                               near edge of plot to reject
#
# Output:
#
#	A csv file of the local maxima co-ordinates
#
#	Jonathan Williams
#	jonvw28@gmail.com	
#	06/11/2017		

# Catch arguments to the script
args <- commandArgs(trailingOnly = FALSE)
idx <- grep('--file=',args)
script.path <- dirname(normalizePath(sub('--file=', "", args[idx])))
args<- args[-(1:(idx+1))]

# Check that these make sense
if (length(args) == 0){ 
        # Has Input CHM Been Specified
        stop("Input not specified")
        
} else if(length(args) == 1){ 
        # Check if CHM scale is set
        warning("CHM scale has not been specified, defaulting to 1, and minimum height not specified, defaulting to 10")
        args <- c(args,1,10)
        
} else if(length(args)>=2 && is.na(as.numeric(args[2]))){
        # Check CHM scale is valid
        stop("Scale must be a number")
        
} else if(length(args)>=3 && is.na(as.numeric(args[3]))){
        # Check minimum height is valid
        stop("minimum height must be a number")
        
} else if(length(args)==4 && is.na(as.numeric(args[4]))){
        # Check if Tau is set
        warning("Tau not set, defaulting to 90")
        args <- c(args,90)
} else if(length(args)==5 && is.na(as.numeric(args[5]))){
        # Check Tau is valid
        stop("tau must be a number")        
}
        
# Check folder is unix friendly
args[1] <- gsub("\\\\","/",args[1])
args[4] <- gsub("\\\\","/",args[4])

# get package dependencies
if(!require(raster)){
        install.packages('raster')
}
library(raster)
library(rgdal)

# Get required scripts
source(paste(script.path,'/tom_itc_max/utils.r',sep=""))
source(paste(script.path,'/tom_itc_max/detect_maxima.r',sep=""))
source(paste(script.path,'/tom_itc_max/gen_ITCsegmax_prior.r',sep=""))

# Find maxima based on lookup table or not
if(length(args)==3){
        maxima <- gen_ITCSegmax_prior(args[1],
                                      scale=as.numeric(args[2]),
                                      min_h=as.numeric(args[3]))
} else if(length(args)==4){
        maxima <- gen_ITCSegmax_prior(args[1],
                                      scale=as.numeric(args[2]),
                                      min_h=as.numeric(args[3]),
                                      window_size = as.numeric(args[4]))
} else {
        maxima <- gen_ITCSegmax_prior(args[1],lut=args[4],
                                      scale=as.numeric(args[2]),
                                      min_h=as.numeric(args[3]),
                                      tau=as.numeric(args[5]))
}
colnames(maxima) <- c('X','Y','Z')
write.table(maxima,file=paste(substr(args[1],1,nchar(args[1])-4),"_prior.csv",sep=""),
            sep=",",col.names = T, row.names = F)
