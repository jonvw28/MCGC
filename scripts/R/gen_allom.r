gen_allom <- function(data,tau,fileName,zone=NULL,biome=NULL){
#
# Function that generates qunatile regression allometry for the data from
# Jucker et al 2016 paper on allometric relations for different biomes
#
# The function takes the Height and CD and performs logarithmic quantile
# regression on these. The model fits CD as a*H^b. The parameter tau sets which
# quantile of the data is to be used (eg 0.95 for 95-th qunatile). This is
# designed to be an upper bound for allometry. fileName then sets the prefix to
# be included in the output file eg "neotropical_tropic".
#
# When set, zone and biome subset the data to only use those data from a 
# specific biogeographical zone and biome respectively.
#
# Arguments
#
#       data:           .csv file of data from paper
#       tau:            Quantile to be used (eg 0.95 for 95-th)
#       fileName:       prefix for file name of allometry table to be saved
#       zone:           Name of biogeographic zone as set out in the data
#       biome:          Name of biome as set out in the the data
       
        
        
        # Load Dependancies
        if(!require("quantreg")){
                install.packages("quantreg")
        }
        library(quantreg)
        
        # Subset relevant data if requested
        alo_data <- read.csv(data,stringsAsFactors = F)
        if(!is.null(zone)){
                reg_data <- subset(alo_data,alo_data$Biogeographic_zone == zone)
        } else {
                reg_data <- all_data
        }
        if(!is.null(biome)){
                biome_data <- subset(reg_data,reg_data$Biome == biome)
        } else {
                biome_data <- reg_data
        }
        
        # Build log-log quantile regression model
        mdl <- rq(formula = log(CD) ~ log(H), data = biome_data, tau = tau)
        max_H = max(biome_data$H)
        X <- 0:ceiling(max_H)
        Y = exp(mdl$coefficients[1])*X^(mdl$coefficients[2])
        
        # Save the results out as Crown Diameter and Radius respectively vs H
        write.table(cbind(X,Y),file=paste(fileName,'_HD_',as.character(round(100*tau)),".csv",sep=""),col.names = F,row.names = F,sep=",")
        write.table(cbind(X,Y/2),file=paste(fileName,'_HR_',as.character(round(100*tau)),".csv",sep=""),col.names = F,row.names = F,sep=",")
}