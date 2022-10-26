#################################################################
### Compute similarity between two set of corrected sample
#################################################################
## 
## Historic - version- authors
##############################
# v1.0 - Boris Droz - October 2022

# PRODOM project - School of Biological, Earth and Environmental Sciences
# Environmental Research Institute (ERI)
# University College Cork
#################################################################################
## Descriptif
##############
## compute similarity analysis with staRdom
## Please see details in: U.J. Wünsch, R. Bro, C.A. Stedmon, 
## P. Wenig, K.R. Murphy, Emerging patterns in the global 
## distribution of dissolved matter fluorescence, Anal. Methods, 11 (2019), 
## pp. 888-893

## Load library version and requirement.
#######################################
# work under R-4.2.0 [64bits]

library(staRdom) # 1.1.25
library(eemR) # 1.0.1

##################################################################################
## set parameter
################
# Setting working directory
workdir <- "C:/Users/bodro/Documents/Actual/EPA 2019-W-MS43 - PRODOM/Data_analysis/PARAFAC analysis"

folder <- "/output/EEM_corr/"

# gave a list of samples to compare
list.rep <- list(c("DR3_subsoil_220621","DR3_subsoil_220621_R"),
                 c("DY12_010721","DY12_010721_R"),
                 c("DY12_150721","DY12_150721_R"),
                 c("DY12_190521", "DY12_190521_R"),
                 c("DY12_210421", "DY12_210421_R"),
                 c("DY12_210721","DY12_210721_R"),
                 c("DY12_220621","DY12_220621_R"))

#########################################################################
##########################################################################
###########################  ANALYSIS ###########################
#################################################################################
#################################################################################
setwd(workdir)

# get all file to compare
output <- NULL
for (i in 1:length(list.rep))
{
  
  #### Read EEM data ####
  sample_dir1 <- paste(workdir,folder,list.rep[[i]][1],".csv",sep="")
  eem1 <- eem_read(sample_dir1, recursive = TRUE, import_function = eem_csv)
  
  sample_dir2 <- paste(workdir,folder,list.rep[[i]][2],".csv",sep="")
  eem2 <- eem_read(sample_dir2, recursive = TRUE, import_function = eem_csv)
  
  #ssc(eem1[[1]]$x, eem2[[1]]$x, tcc = FALSE) # return ssc
  
  output <- c(output, ssc(eem1[[1]]$x, eem2[[1]]$x, tcc = TRUE))# return TCC
  
  }

mean(output,na.rm=T)
sd(output,na.rm=T)


