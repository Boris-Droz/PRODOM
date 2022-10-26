####################################################################################################################################
###
### Tcheck UV-absorbance in the range of our interest
####
#####============================================================
## Historic - version- authors
##############################
# v1.0 - Boris Droz - June 2022

# PRODOM project - School of Biological, Earth and Environmental Sciences
# Environmental Research Institute (ERI)
# University College Cork
#################################################################################
## Descriptif
##############
# Check the values of Abs for all sample within the range of interest
###############################################################################

## NEEDS ###
############
# input folder:
# - /input/ABS

## Load library version and requirement.
#######################################
# work under R-4.2.0 [64bits]

## set parameter
################
# Setting working directory
workdir <- "C:/Users/bodro/Documents/Actual/EPA 2019-W-MS43 - PRODOM/Data_analysis/PARAFAC analysis"

# range of abs interest
abs_range = c(250, 600) 

#################################################################################
########################### Preprocessing ABS spectra ###########################
#################################################################################
#################################################################################
setwd(workdir)

# make folder list (should be similar in EEM and ABS folder)
list.folder <- list.dirs(path = paste(workdir,"/input/EEM",sep=""), full.names = FALSE, recursive = FALSE)

d.out <- NULL
for ( k in 1:length(list.folder))
    {
    ################################################################################
    ## Read absorbance data ####
    absorbance_dir = paste(workdir,"/input/ABS/",list.folder[k], sep="")
    abs_fns <- list.files(absorbance_dir)
    
    for (q in 1:length(abs_fns))
        {
        abs_data <- read.table(paste(absorbance_dir,"/",abs_fns[q], sep="")
                               , header=F, sep=",")
        
        abs_data <- abs_data[abs_data[,1]>=abs_range[1]&
                               abs_data[,1]<=abs_range[2] ,]
          
        d.out <- rbind(d.out,data.frame(sample=abs_fns[q],
                                        min_abs=min(abs_data[,2]),
                              max_abs= max(abs_data[,2]) ) )
        }
    } # end list folder  

d.out <- rbind(data.frame(sample="TOTAL",
                                min_abs=min(d.out$min_abs),
                                max_abs= max(d.out$max_abs) ),d.out )

# data table of all values
write.table(d.out, file=paste(workdir,"/output/ABS_values.txt",sep=""),
            sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)

###########################################################################################
###########################################################################################
###########################################################################################
###### GO FOR COFFE -- no time for that!!!
#######################
#######################

