#################################################################
### ROOT MEAN SQUARE ERROR OF ALL BLANK
########################################
## 
# Historic
##########
# v1.0 - Boris Droz - June 2022

# PRODOM project - School of Biological, Earth and Environmental Sciences
# Environmental Research Institute (ERI)
# University College Cork
##################################################################
## Descriptif
##############
## root mean square error (RMSE) between the matrices 
## can be estimated (Burdick  and Tu, 1989). If, All matrix correlations had RMSE 
## alues below 0.06 indicating a high degree of similarity between Milli-Q® 
## spectra and supporting instrument and hardware stability during the 
## experimental period. (Peleato, N. M.; et al. J. Environ. Sci. 2015, 
## 27, 159-167. https://doi.org/10.1016/j.jes.2014.04.014.)

## Load library version and requirement.
#######################################
# work under R-4.2.0 [64bits]

##################################################################################
## set parameter
################
# Setting working directory
workdir <- "C:/Users/bodro/Documents/Actual/EPA 2019-W-MS43 - PRODOM/Data_analysis/PARAFAC analysis"

#########################################################################
##########################################################################
###########################  ANALYSIS ###########################
#################################################################################
#################################################################################
setwd(workdir)

f.dir <- list.dirs(paste(workdir,"/input/EEM",sep=""),recursive=FALSE)

# get all file with blank name
fns <- NULL
for (i in 1:length(f.dir))
{fns <- c(fns,list.files(f.dir[i], full.names=TRUE, pattern="blank") )}

sel.id <- c(1:16,18:22)

# data list
# make the mean
for (i in sel.id) ## dont know why data i==17 as less row??
  {m <- read.csv(fns[i], sep=",",header=F)
    m <-as.matrix(m[2:nrow(m),2:ncol(m)])
    if (i==1) 
      {m.mean <- m}else{m.mean <- m+m.mean}
  }
m.mean <- m.mean/length(sel.id)

for (i in sel.id) ## dont know why data i==17 as less row??
  {m <- read.csv(fns[i], sep=",",header=F)
  m <-as.matrix(m[2:nrow(m),2:ncol(m)])
  m <- (m-m.mean)^2
  if (i==1) 
  {m.sum <- m}else{m.sum <- m+m.sum}
  }

# sum all pixel
sqrt(sum(m.sum)/(length(sel.id)*696*35))

rmse.m <- sqrt(m.sum/length(sel.id))
mean(rmse.m)
sd(rmse.m)


