################################################################################
###############################################################################
## Clean and Correct time serie from sensor
############################################
## Used in the case of clean optical sensor who creat some data shift
## + some glitch on the data (saturation and trend in baseline)
### .......
################################################################################
#
################################################################################
################################################################################
## Historic
## --------
## v1.0 August 2023 B.DROZ

# PRODOM project - School of Biological, Earth and Environmental Sciences
# Environmental Research Institute (ERI)
# University College Cork
################################################################################
# DESCRIPTION
###############
## RUN under R-4.3.1
################################################################################
## Load library
################
library(bcp)
library(strucchange)
library(segmented)
library(tree)
library(tidyr)
################################################################################
## SCRIPT PARAMETER  
##  --> NEED TO GO TROUGH AND MAKE THE APPROPRIATE MODIFICATION !!!
#####################
# Set folder 
## =========
workdir <- "C:/Users/bodro/Documents/Actual/EPA 2019-W-MS43 - PRODOM/Data_analysis/Vlux_continuous_data"

f.sensor <- "/input/VLux_data_20211001-20220908.csv"

# f.sensor <- "/input/VLux_data_20211001-20220908_clean.csv"
  
## select only part of the data
###############################
t.start <- "2021-10-05 07:00:00 BST"
t.end <- "2022-07-31 07:00:00 BST"

################################################################################
## SCRIPT STAR HERE
################################################################################
# read data
data.in <- read.table (paste(workdir,f.sensor, sep=""),
                         header=TRUE, sep =","); head(data.in)   

# prepare time data
data.in <- cbind(data.in, iso_time = as.POSIXct(paste(data.in$DATE, data.in$TIME), 
                                                format="%d-%m-%y %H:%M:%S"))
  
#select part of the data to test 
data.in <- data.in[data.in$iso_time>t.start& data.in$iso_time<=t.end, ]

################################################################################
## Start data analysis
##################
### CLEAN OUTLIER -- conservative version
##################
#correct/remove extrem value data -- base on limit of the sensor
# fill NA values with previous data
data.in$TRYPTOPHAN[data.in$TRYPTOPHAN<=3.2] <- 0
data.in$TRYPTOPHAN[data.in$TRYPTOPHAN>600] <- NA
data.in <- data.in %>% fill(TRYPTOPHAN) # fill NA

data.in$FDOM[data.in$FDOM<=3.2] <- 0
data.in$FDOM[data.in$FDOM>600] <- NA
data.in <- data.in %>% fill(FDOM) # fill NA

data.in$ABSORBANCE[data.in$ABSORBANCE<=0.25] <- 0
data.in$ABSORBANCE[data.in$ABSORBANCE>2] <- NA
data.in <- data.in %>% fill(ABSORBANCE)

data.in$TURBIDITY[data.in$TURBIDITY<=17] <- 0
data.in$TURBIDITY[data.in$TURBIDITY>1000] <- NA
data.in <- data.in %>% fill(TURBIDITY)

####################################################
data.in <- data.in[!duplicated(data.in$iso_time), ] # DELET duplicate data in time
##########################
# smooth data 
#############
## Simple Moving Average
df <- zoo(data.in[,3:7],data.in$iso_time)
df <- rollmean(df,50)
df <- as.data.frame(df)

df <- data.frame(TIME=as.POSIXct(rownames(df)),df )

# write table out 
write.table(df, file=paste(workdir,"/input/VLux_data_20211001-20220908_clean.csv", sep=""),
            sep=",", append=FALSE, row.names=FALSE,col.names=TRUE, quote=FALSE)



