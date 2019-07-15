mcgc_cleaner <- function(mcgc_dir){
#
# A function that clears outr all .las and .csv files in a given layer of the
# mcgc gridded algorithm. It should leave any error reports along with the
# .txt files for the output.
#
# mcgc_dir should be the directory containing the directories 'pass_1', 'pass_2'
# etc and needs to be supplied as an absolute path
#
# Jonathan Williams 16/03/2018
#

# move to directory to be cleared
setwd(mcgc_dir)
        
# remove las files if any are here (from pts left)
lases <- list.files(pattern = "\\.las$")
file.remove(lases)

# Clear all but error reports and final output for each pass
passes <- list.dirs(recursive=F)

for(dir in passes){
        if(dir=="."){
                next       
        }
        setwd(paste(mcgc_dir,'/',dir,sep=''))
        csvs <- list.files(pattern = "\\.csv$")
        file.remove(csvs)
        lases <- list.files(pattern = "\\.las$")
        file.remove(lases)
}
}

mcgc_2_layer_cleaner <- function(mcgc_dir){
#
# A function that clears out all .las and .csv files in a given double layer run
# of the mcgc gridded algorithm. It should leave any error reports along with the
# .txt files for the output.
#
# mcgc_dir should be the top level directory (ie the one containing the
# directories 'layer_1' and 'layer_2' and needs to be supplied as an absolute path
#
# Jonathan Williams 16/03/2018
#

# move to directory to be cleared
setwd(mcgc_dir)

layer_1 <- paste(mcgc_dir,'/layer_1',sep='')
mcgc_cleaner(layer_1)

layer_2 <- paste(mcgc_dir,'/layer_2',sep='')
if(file.exists(layer_2)){
        mcgc_cleaner(layer_2)
}
}