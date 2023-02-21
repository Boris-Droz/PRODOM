####################################################################################################################################
####################################################################################################################################
###					
###  TEST NORMALITY
### ===========================
###
###-->  variable	transformation to have normal data  
####################################################################################################################################
####################################################################################################################################
## v1.1, June 2015 - Boris DROZ & Gerrad Jones, ETHZ & EAWAG - Projet Se around the world
## v2.0 November 2016 - Boris Droz - repair bug from v.1.1 and add transf
## v3. July 2022 modified for the PRODOM project

# PRODOM project - School of Biological, Earth and Environmental Sciences
# Environmental Research Institute (ERI)
# University College Cork
####################################################################################################################################
## DESCRIPTION
###############
# Performed given transformation for a given data set
# Additionnally performed a box cox transformation
# Nomality for each transformation are asses by Shapiro test, hist plot and Q-Q plot
##                            shapiro W value close to 1 indicate close to normality
# Additionnaly two test are also reported, where 0 correspond at symmetric and unimodal distributions for both.
#              skew value: asymetrie test: neg. asymetrie on left and pos. asymetrie on right
#              kurtosis value:  = "peakedness" of the distribution and the heaviness of its 
#              tail (in some point see wiki exemple)
# Input data:
#############
# - data set
# - var.list --> defined all variable within the data
# - remove list --> defined during the PARAFAC analysis -- just copie it on the appropriate folder


# Output file:
##############
# - Table with all parameter to chooose the best transformation
# - histogramme and q-q plot for each transformation
#
####################################################################################################################################
## Load library
##################
# work under R-4.2.0 [64bits]
## Set library path
# .libPaths("C:/Users/Public/Documents/R/win-library/3.1")

library(rms)  
library(pracma)
library(moments)

####################################################################################################################################
## SCRIPT PARAMETER
#####################
## FIRST creat an input and output folder on the workdir
workdir <- "C:/Users/bodro/Documents/Actual/EPA 2019-W-MS43 - PRODOM/Data_analysis/ML_pred"

# input folder
f.dat <- "input_data_MLv3.csv"# File name of data points
f.var <-"var.list.txt"# File with the variable list
f.rem <- "exclude.list.txt" # File with the list of the sample to remove

REM.ZERO <- "YES" # option to remove zero -- > needed if you use SVC in the model.proj script fo y var.

####################################################################################################################################
####################################################################################################################################
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
####################################################################################################################################
####################################################################################################################################
## SCRIPT STAR HERE
####################################################################################################################################

# Set folder --> 
date <- Sys.Date()
inpath <- paste(workdir,"/input/", sep="")
outpath <-paste(workdir,"/output/",date,"_norm.transf/",sep="")
creat.subDir( paste(workdir,"/output",sep=""), paste("/",date,"_norm.transf", sep="") )

## set your work space
setwd(inpath)

## TRANSFORMATION list
#######################
## Follow recommendation of Velleman & Hoaglin (1981)
############
transf <- list( alist(x=,(x)), alist(x=,log10((abs(x)))),
                alist(x=,log10((abs(x)+1))), # cool if a lot of data with zero
                alist(x=,(x^2)),alist(x=,sqrt(abs((x)))), 
                # inverse also good to test
                alist(x=,(1/x)),alist(x=,(1/log10(abs(x)))), alist(x=,(1/x^2)), alist(x=,(1/sqrt(abs(x))))
                ) ; length(transf)

## text of the transformation sim. as transf.
text.transf <- c( "x","log(x)","log(x+1)","x^2","sqrt(x)", 
                  "inv(x)","inv(log(x))", "inv(x^2)","inv(sqrt(x))" ) 
              length(text.transf)

##############################################
# (1) Loading  datasets #
##############################################
# Open the data set and the unknow point around the world
data.in <- read.table(paste(inpath,f.dat,sep=""),header = TRUE ,
                      na.strings = "NA", sep=",")

var.list <- read.table(paste(inpath,f.var,sep=""), header= TRUE)
rem.list <- read.table(paste(inpath,f.rem,sep=""), header= FALSE)

########################################
## Creat Organized data set for the test
## --> remove outlier defined during parafac analysis
## --> select only variable present in the var.list 
var.list <-na.omit(var.list)
data.in <- data.in[!(data.in$sample %in% rem.list$V1),]
pos <-match(var.list$var.name, names(data.in))
data.in <- data.in[,pos]

# delet NA
data.in <- na.omit(data.in) # delet row with NA value
nrow(data.in)

#########################################
## APPLIED TRANSFORMATION ON THE DATA SET
#########################################

test.table <- NULL

for (t in 1:length(transf))
    {
      shapw <- NULL
      pval <- NULL
      skew <- NULL
      kurt <- NULL
      
        for ( n in 1:ncol(data.in) ){
            ## FOR EACH VARIABLE
            # TRANSFORMATION
          if (REM.ZERO =="YES")
          {
            data.sel <-data.in[,n]
            data.sel <- data.sel[data.sel!=0]
            data.sel <- as.function(transf[[t]])(data.sel)
          }else{
            data.sel <- as.function(transf[[t]])(data.in[,n])
          }
          
            # change Inf by NA
            is.na(data.sel) <- sapply(data.sel, is.infinite)
          
		        # Select data and DEALET NA values
      	    data.sel <- na.omit(data.sel)

            ## PLOT HIST AND Q-QPLOT FOR EACH DATA
            jpeg(paste(outpath,"/",names(data.in)[n],"_hist_qplot_",text.transf[t],".jpg",sep=""),width = 60, height = 30,units="cm",res=150)
            
              par(mfrow=c(1,2)); hist(data.sel, breaks="Sturges", freq = FALSE); qqnorm(data.sel); qqline(data.sel)
            
            dev.off() 
            
            ## TEST NORMALITY FOR ALL VARIABLE
            shapw <- c(shapw,shapiro.test(data.sel)$statistic)
            pval <- c(pval,shapiro.test(data.sel)$p.value)
            
            ## TEST SKEW AND KURT
            skew <- c(skew,skewness(data.sel,na.rm =TRUE))
            kurt <- c(kurt,kurtosis(data.sel,na.rm =TRUE))
            
            }
      test.table <- rbind(test.table,numb = length(data.sel), Shap = shapw, p_val = pval, skew = skew, kurt = kurt )
      
      }
      
dimnames(test.table) <- list(paste(rep(text.transf,each=5),c("_numb","_Shap_W","_P_value","_skew","_kurt")),
                             dimnames(data.in)[[2]])

### max shapiro -- best transf..
pos <- seq(from=2, to =(length(text.transf)*5)-3, by=5)
#pos <- seq(from=2, to =(10*5)-3, by=5)

max.shapi <- apply(test.table[pos,],2,max)
col <- match(max.shapi,test.table[pos,])- seq(from=0,to=(ncol(test.table)-1),by=1) * nrow(test.table[pos,])

test.table <- rbind(best.shapi = text.transf[col], test.table)
    
write.table(test.table,file = paste(outpath,"/AA_Normality_test.txt",sep=""),
            sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)      


