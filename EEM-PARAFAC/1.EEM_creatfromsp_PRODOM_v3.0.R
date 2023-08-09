####################################################################################################################################
###
### Excitation-emission matrice fluorescence spectra
##  convert sp file from perkin-elmer into a EEM to be further analyzed with staRdom package
####
#####============================================================
## Historic - version- authors
##############################
# v1.0 - Elena Fernandez December 2021
# v2.0 - Modif by Boris Droz - June 2022
## v3.0 - Integrate automatic read of sp data (PErkin Elmer) - Stefan Seeger

# PRODOM project - School of Biological, Earth and Environmental Sciences
# Environmental Research Institute (ERI)
# University College Cork
#################################################################################
## NEEDS ###
############
# input data are in two folder at the workdir roots
# - /input/EEM_raw
# self-creat an input folder with EEM for further calculation

## Load library version and requirement.
#######################################
# work under R-4.2.0 [64bits]

## set parameter
################
# Setting working directory
workdir <- "C:/Users/bodro/Documents/Actual/EPA 2019-W-MS43 - PRODOM/Data_analysis/Vlux_validation"

# set parameter similar to header and column in matrix() min-max range and increment)
ex.para <- c(250,590,10)
em.para <- c(300, 647.5,0.5)
#################################################################################################################################################
# FUNCTION check and produced subDir folder
###########################################
#February 2017 - B.Droz

creat.subDir <- function (mainDir,subDir)
{
  if (file.exists(subDir)){
    setwd(file.path(mainDir, subDir))
  } else {
    dir.create(file.path(mainDir, subDir))
    setwd(file.path(mainDir, subDir))
  }
}

################################################################################
###### R-Toolbox to open spectroflourometer (sp files) data obtained with 
##  FLWInlab software ######
#################################################
# authors: Stefan Seeger & Elena Fernandez Pascual
# Created: 2015-11-18
# Last revision: 2023-07-17
source(paste(workdir,"/script/sp_convertion_toolbox.r",sep=""))

#################################################################################
########################### Preprocessing EEM spectra ###########################
#################################################################################
#################################################################################
setwd(workdir)

# make folder list (should be similar in EEM and ABS folder)
list.folder <- list.dirs(path = paste(workdir,"/input/EEM_raw",sep=""), full.names = FALSE, recursive = FALSE)

#d.out <- NULL
for ( k in 1:length(list.folder))
    {
  
  # self-creat input-outpu folder
  creat.subDir (paste(workdir,"/input/",sep=""),"EEM")
  creat.subDir (paste(workdir,"/input/EEM/",sep=""),list.folder[k])
  
  cat("PROCESS",list.folder[k], "\n",append = F) 
  ################################################################################
  #### Read and convert sp data ####
  sample_dir <- paste(workdir,"/input/EEM_raw/",list.folder[k],sep="")
  
  sampl_names <- unique( unlist(sapply( strsplit(list.files(sample_dir),"#") ,"[[",1) ) )
  
  for (l in 1: length(sampl_names))
      {
      # function to load all FLWinlab sp files and combine them into one matrix
      filenameList <-dir(sample_dir)[ grep(sampl_names[l], dir(sample_dir)) ]
      filenameList <- filenameList[grep("\\.sp$", filenameList)]
      eem <- loadSpList(sprintf("%s/%s", sample_dir, filenameList))
      
      ## Add the excitation and emission wavelengths in the created matrix:
      row_names <- c(NA , seq(ex.para[1], ex.para[2], by=ex.para[3]) )
      col_names <- seq(em.para[1], em.para[2], by=em.para[3])
      
      # for each sample
      Sample1_col <- eem[["dataMatrix"]]
      Sample1_row <- cbind (col_names, Sample1_col)
      eem <- rbind(row_names, Sample1_row)
      write.table(eem, file= paste(workdir,"/input/EEM/",list.folder[k],"/",sampl_names[l],".csv",sep=""),
                     sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
}

###########################################################################################
###########################################################################################
###########################################################################################
###### GO FOR COFFE
#######################
#######################

