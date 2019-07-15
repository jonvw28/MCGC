Pts2TreeShp <- function(datafile,shape_dir,shape_name,class_col=8,x_col=1,
                        y_col=2,z_col=3){
#
# A function that takes a pointcloud dataset in csv format and returns
# a ShapeFile of the convex hulls of each tree labelled by a given 
# number in class_col. Each shape also has associated x,y,z data of
# the top point of each tree in the dbf file
#
# Data is then saved out to the ShapeFile layer as set by the arguments
# shape_dir and shape_name
#
# Arguments
#
#       datafile:       .csv file of points
#       shape_dir:      name of directory to be created to store ShapeFile data
#       shape_name:     name of ShapeFile layer to be created
#       class_col:      column number containing numerical label for tree 
#                       membership of each point (use 0 if not a member of tree)
#       x_col:          column number containing x co-cordinates of each point
#       y_col:          column number containing y co-cordinates of each point
#       z_col:          column number containing z co-cordinates of each point        

# Some input checking
if(!is.character(shape_dir)){
        stop('Argument shape_dir must be a string')
}
if(!is.character(shape_name)){
        stop('Argument shape_name must be a string')
}        

# get dependencies
if(!require(sp)){
        install.packages('sp')
}
library(sp)
if(!require(rgdal)){
        install.packages('rgdal')
}
library(rgdal)

# read data
data <- read.csv(datafile,header=F)

# remove unclassified points
zero_idx <- which(data[,class_col]==0)
if(length(zero_idx!=0)){
        data <- data[-zero_idx,]   
}

# If no points classified then end here
if(nrow(data)==0){
        return(NA)
}

# Get tree indices
tr_ids <- sort(unique(data[,class_col]))

# Set up objects to collect data
all_poly <- list()
x <- numeric(length(tr_ids))
y <- numeric(length(tr_ids))
z <- numeric(length(tr_ids))

# get each tree
for(i in 1:length(tr_ids)){
        
        #get tree id
        j <- tr_ids[i]
        
        # Get convex hull indices
        ths_tr <- data[data[,class_col]==j,]
        ths_hl <- chull(ths_tr[,c(x_col,y_col)])
        ths_hl <- c(ths_hl,ths_hl[1])
        
        # convert to spatial polygon, as an island not hole
        ths_poly <- sp::Polygon(cbind(ths_tr[ths_hl,1],ths_tr[ths_hl,2]),hole=F) 
        ths_polys <- sp::Polygons(list(ths_poly),as.character(i))
        
        # Add to our list
        all_poly[[length(all_poly)+1]] <- ths_polys
        
        # Get location of tree top
        ths_mxid <- which.max(ths_tr[,z_col])
        x[i] <- max(ths_tr[ths_mxid,x_col])
        y[i] <- max(ths_tr[ths_mxid,y_col])
        z[i] <- max(ths_tr[ths_mxid,z_col])
        
        # clean up
        rm(ths_tr,ths_hl,ths_poly,ths_polys,ths_mxid)
}

# Compiles all polygons and data in a spdf
all_polys <- sp::SpatialPolygons(all_poly,1:length(tr_ids))
final_polys <- sp::SpatialPolygonsDataFrame(all_polys,
                                            data=data.frame(x=x, y=y,z=z,
                                                            row.names=row.names(all_polys)))

# Save out the result
rgdal::writeOGR(obj=final_polys,dsn=shape_dir,layer=shape_name,driver="ESRI Shapefile")
}