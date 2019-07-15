#' Find local maxima of a raster CHM
#'
#' @param x A Matrix encoding a canopy height model raster
#' @param scale The resolution of the raster in meters
#' @param lm.searchwin The search window for finding maxima. If this is set to NULL the window will scale
#' allometrically according to tree height
#' @param tau The percentile to use in calculating the tree maximum size near egdes of image
#' @param min_h The minimum height for data to be included
#' @return A matrix of coordinates of tree maxima locations in pixels and tree heights
#' @export
#' @author Tom Swinfield (edited by Jon Williams)
#' @details
#'
#' Created 16-10-10
#' Edited 17-11-02

#x<-img
#scale = res(imagery)[1]
#tau = 90

detect_maxima <- function(x, scale, lut = NULL, lm.searchwin = NULL, tau = 90, min_h = 0) {
        
        # Creates an empty data frame to store the maxima locations:
        Max <- x
        Max[, ] <- 0    # Sets max to 0s.
        
        y <- which(x >= min_h, arr.ind = T) # Extracts the locations of pixels which are bigger than the min tree height.
        z<-x[y[,1]+nrow(x)*(y[,2]-1)] # extracts the tree heights from the matrix
        
        # If using a variable window, set this up and work out distance to edge
        # to ignore
        if (class(lm.searchwin) == 'NULL'){
                WinSize<-allom.dist(z, scale=scale, lut, tau = tau)
                pix2edge<-edge.dist(z, scale=scale, lut, tau = tau)
        } else {
                pix2edge <- floor(lm.searchwin/2)
        }
                
        # Extracts only the pixels far enough from the image edge.
        y.ind<-y[, 1] > pix2edge &
                y[, 1] < (nrow(x) - pix2edge) &
                y[,2] > pix2edge &
                y[,2] < (ncol(x) - pix2edge)
        
        y<-y[y.ind,]
        z<-z[y.ind]
        
        if (class(lm.searchwin) == 'NULL'){
                WinSize<-WinSize[y.ind]
        }
                
        length.total<-length(z)
        
        # Throw an error if there simply aren't enough point to generate a maxima
        if(length.total<10){
                stop('Not enough pixels in CHM to generate maxima')
        }
        
        coordSeeds <-
                matrix(NA,
                       nrow = length.total,
                       ncol = 2,
                       dimnames = list(NULL, c('row', 'col')))
        treeHeights <- vector('numeric', length = length.total)
        
        Ntrees <- 0 # the number of trees encountered
        
        # Works through each pixel one by one:
        for (i in 1:nrow(y))
        {
                r = as.numeric(y[i, 1]) # Extracts the row pos.
                k = as.numeric(y[i, 2]) # Extracts the column pos.
                # Sets the search window size for the local maxima:
                if (class(lm.searchwin) == 'NULL')
                        searchWinSize <-
                        WinSize[i] # a search window scaled by max tree height.
                else
                        searchWinSize <- lm.searchwin # a fixed size search window.
                
                hc.sWS <-
                        ceiling(searchWinSize / 2) # half the search window size rounded to floor
                hf.sWS <-
                        floor(searchWinSize / 2) # half the search window size rounded to ceiling
                FIL <- matrix(searchWinSize, searchWinSize, data = NA)
                # Extracts the search window for the pixel:
                FIL <- x[(r - hf.sWS):(r + hf.sWS),
                         (k - hf.sWS):(k + hf.sWS)]
                # Extracts the window from Max indicating whether a tree has already being designated or not:
                Max.chk <- Max[(r - hf.sWS):(r + hf.sWS),
                               (k - hf.sWS):(k + hf.sWS)]
                
                # If the focal pixel has the greatest value in the window & there is no tree already assigned in the output matrix within the window & the max value is no 0...
                # because the order is from tallest to shortest, large trees will always suppress the designation of small trees.
                if (FIL[hc.sWS, hc.sWS] == max(FIL, na.rm = T) &
                    max(Max.chk, na.rm = T) == 0 &
                    max(FIL, na.rm = T) != 0)
                {
                        Ntrees <- Ntrees + 1 # increments
                        Max[r, k] <- 1 # A logical assignment of tallest tree
                        
                        coordSeeds[Ntrees,] <-c(r, k) # Assigns the sequence in which the trees were found.
                        treeHeights[Ntrees] <- x[r, k] # Assigns the height of the trees
                }
        }
        cat(Ntrees, 'trees detected; creating spatial objects\n')
        if(Ntrees>1){
                coordSeeds <-
                        coordSeeds[1:Ntrees,] # Extracts the tree crowns, which are already in order.
                coordSeeds <- cbind(coordSeeds, height = treeHeights[1:Ntrees]) # Extracts the tree crowns, which are already in order.
                return(coordSeeds)
        } else if(Ntrees==1){
				outCoords <- matrix(NA,
                       nrow = 1,
                       ncol = 3,
                       dimnames = list(NULL, c('row', 'col', 'height')))
                outCoords[1,] <- c(coordSeeds[1,],treeHeights[1]) # Extracts the tree crowns, which are already in order.
                return(outCoords)
		} else {
                return('Notrees')
        }
     
}