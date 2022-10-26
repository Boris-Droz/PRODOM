####################################################################################################################################
###
### Compare two fluorescence excitation−emission matrix (EEM) together.
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
# Two EEM csv data

## Load library version and requirement.
#######################################
# work under R-4.2.0 [64bits]

##############################################################################
######################### START RUN #########################################
##############################################################################
# compare two run for comparing indice matrix------
d0 <- read.table(file.choose(), header=TRUE, sep = ",")
d1 <- read.table(file.choose(), header=TRUE, sep = ",")

# in case of not exactly same matrix to compare
###############################################

d0 <- d0[,1:17] # for comparing indice matrix
d1 <-d1[,1:17]
d1<-na.omit(d1)

d0 <- d0[match(d1[,1],d0[,1]), ] # do since here to compare parafac resut
d0 <- na.omit(d0)

d1 <- d1[match(d0[,1],d1[,1]), ]
d1 <- na.omit(d1)

#################################################
# compare result
dcomp <- (d1[,2:17]-d0[,2:17])/d0[,2:17] # for indice matrix
dcomp <- (d1[,2:7]-d0[,2:7])/d0[,2:7] # for parafac matrix

boxplot(dcomp)

summary(dcomp)
