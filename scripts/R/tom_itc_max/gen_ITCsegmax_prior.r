#' Function which reads in a CHM TIFF file and produces a table of (x,y,z) of 
#' local maxima based on methods from ITCsegment package
#' 
#' @param chm_tiff the filename of the CHM as a TIFF file
#' @param scale The resolution of the raster CHM in meters
#' @param lut_csv An optional lookup table for H to R allometry based on quantile
#' quantile regression. Can be generated by rq_lut() method
#' @param window_size The fixed window size (in terms of pixels) to be used if 
#' lut not supplied
#' @param min_h The minimum height for data to be included
#' @return A matrix of coordinates of local maxima (x,y,z)
#' @export
#' @author Jon Williams
#' @details
#'
#' Created 17-11-06
#
# NB: requires both utils.R and detect_maxima.r from the work of Tom Swinfield
#
gen_ITCSegmax_prior <- function(chm_tiff,scale,lut_csv=NULL,window_size=3,tau=90,min_h=0){
        
        # get package dependencies
        if(!require(raster)){
                install.packages('raster')
        }
        library(raster)
        
        # read in CHM and get cooridnates of [0,0] on grid
        chm <- raster(chm_tiff)
        xmin <- chm@extent@xmin
        ymax <- chm@extent@ymax
        chm_mat <- as.matrix(chm)
        
        # If supplying a look-up use it, otherwise used a fixed window size
        if(class(lut_csv) != 'NULL'){
                lut <- read.csv(lut_csv)
                chm_max <- detect_maxima(chm_mat,scale=scale,lut=lut,tau=tau,
                                         min_h=min_h)
        } else {
                chm_max <- detect_maxima(chm_mat,scale=scale,
                                         lm.searchwin = window_size,tau=tau,
                                         min_h=min_h)
        }
        
        # Convert pixels to coordinates- centre on middle of grid
        # Note: matrix indexed from top left, rows for y, cols for x
        
        if(is.character(chm_max)){
                # if no maxima found then link to a near ground point
                maxima = matrix(nrow=1,ncol=3, data=0)
                maxima[1,] = c(xmin,ymax,1)
                
        } else {
                # We have maxima from algorithm
                maxima <- matrix(nrow=nrow(chm_max),ncol=3, data=0)
                for (i in 1:nrow(chm_max)){
                        ths_x <- scale*(chm_max[i,2]-1 + 0.5) + xmin
                        ths_y <- ymax- scale*(chm_max[i,1]-1+0.5)
                        ths_z <- chm_max[i,3]
                        maxima[i,1:3] <- c(ths_x,ths_y,ths_z)
                }
        }
        
        return(maxima)
}