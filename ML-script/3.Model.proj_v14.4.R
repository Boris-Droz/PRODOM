####################################################################################################################################
####################################################################################################################################
###					                ###
###  MODEL PROJECTION		    ###   
###                         ###					
####################################################################################################################################
####################################################################################################################################
## Historic
## --------
## ...
## v7.0 September 2016 - Boris DROZ & Gerrad Jones, ETHZ & EAWAG
## v8.0 Januar 2017 - Boris DROZ --> more ML technic implemented
## v9.0 include xboost and ranger
## v10.0 include extraction outlier (option)
## v12.0 January 2018 - inclu stop criteria --> decrease over-fitting
## v13 August 2022 adapt for non spatial data + update for R4.2.0 and all new package
## v14 October 2022 - support vector classification option

# PRODOM project - School of Biological, Earth and Environmental Sciences
# Environmental Research Institute (ERI)
# University College Cork
####################################################################################################################################
## DESCRIPTION
###############
## RUN under R-4.2.0 
##
## Performed separately several modelling technics (choose by user)
##        Machine learning model
##        ----------------------
##                  - Neural Networks (nnet -- nnet package v.7.3-17)
##                  - Stuttgart Neural Network Simulator (SNNS -- RSNNS package v.0.4-14) 
##                  - Extrem learning machine (ELM -- elmNNRcpp package v. 1.0.4)
##                  - Random Forest (RF -- randomForest package v.4.7-1.1) 
##                  - Case-specific Random Forest (CSRF --- ranger package v 0.14.1)
##                  - Bagging trre method (BAG -- ipred package v.0.9-13)
##                  - Extreme Gradient Boosting (EGB --- xgboost v.1.6.0.1)
##                  - Generalized Boosted Regression Models (GBM -- gbm package v.2.1.8)
##            
##        Simple linear model
##        -------------------
##                  - Generalized Linear Models (GLM)
##                  - GLM step-wise (GLM_sw)
##                  - Generalized Additive Models (GAM --package gam v. 1.20.2)
##             
##
# Input: - data.table : var. to predict and cooordinate for each data point
#=======  
#
# Output: - model projection
#=======  - result from the cross validation
#         - Importance are compute using garson for nnet and rsnns, mean decrease in node impurity for tree method
#             and permutation procedure for linear model
# 
####################################################################################################################################
## Load library
##################

## Set library path
#.libPaths("C:/Users/Public/Documents/R/win-library/3.1")
# Cheak and download or releaod package
## v1.0 Emmanuel Rey 2013 EAWAG
########################################
# check library version and load library
check.lib <- function(package.list) {
  
  for(pkg in package.list){
    
    libTest <- try(library(pkg,character.only=TRUE),silent=TRUE)
    
    if(class(libTest)=='try-error'){
      
      updteTest <- try(install.packages(pkg))   
      
      if(class(updteTest)=='try-error'){update.packages(pkg)}
      
      else{install.packages(pkg)}
      
      libTest <- try(library(pkg,character.only=TRUE),silent=TRUE)
      
      print(pkg)
      
    }
  }
}

## library list
lib.list <- c("ade4","gam","nnet","neuralnet","NeuralNetTools","glmnet",
              "randomForest","AICcmodavg", "RSNNS","elmNNRcpp","gbm",
              "ipred","ranger","xgboost","e1071",
              "rms","foreign","ggplot2","gplots","plyr","Hmisc","corrplot","MASS",
              "AICcmodavg","stringr","boot", "pracma","party","snowfall","parallel")

## ""glmulti", prob --> prob with Rstudio

check.lib(lib.list)
####################################################################################################################################
## SCRIPT PARAMETER  --> NEED TO GO TROUH AND MAKE THE APPROPRIATE MODIFICATION !!!
#####################
# Set folder 
## =========
workdir <- "C:/Users/bodro/Documents/Actual/EPA 2019-W-MS43 - PRODOM/Data_analysis/ML_pred"

# input folder
f.dat <- "input_data_MLv3.csv"# File name of data points
f.var <-"var.list.txt"# File with the variable list
f.rem <- "exclude.list.txt" # File with the list of the sample to remove

## choose the name of the variable selection you will use -- var.list column name with 
## 0 for not use and 1 for use
var.sele <- c("ML_opt3","ML_UV")

## NUMBER OF MODEL BUILDING
# ========================
n.rep <- 100
## SAMPLING of each proportion of data at each replic
#====================================================
sampl <- 80 ## percent kept for analysis --usual 90

# transform the data according to the var.list file
TRANSF <-"YES" ## YES or NO

# rescalling to center the data
RES <- "YES" ## YES or NO

# do K-fold cross validation procedure 
## --> if too loo number of data CV is not recommended
CR.VAL <- "NO"

# suport vector classification
SVC <- "YES"
k.nel <- "radial" # kernel type

# ratio occurrence -- 
## threshold ratio of occurence of non zero data 
# if proportion of zero in cali data set bigger stop it
RA.O <- 0.5

## -- CHOOSE MODELLING TECHN:
## YES or NO
#############################
# machine learning method
NNET <- "YES"
RSNNS <-"NO" 
ELM <- "NO"

RF <- "NO" 
CSRF <- "NO"

BAG <- "YES" 
EGB <- "NO"
GBM <- "YES"

# linear additive
GLM <- "NO"
GLM_sw <- "NO"
GAM <- "NO"

## FUNCTION PARAMETER
## ==================
## PERMUTATION FOR var. opt. Calc. (glm, gam, gbm model)
n.perm <- 100

# xx fold cross validation calibration
fold.cv <- 10

## NNET PARAMETER
##################
## Weight decay Folow recomendation of B. D. Ripley: "Pattern Recognition and Neural Networks", Cambridge, 1996.
## between 0.1 and 0.01
deca <- 0.1

UNIT.Max <- 15 # NUMBER OF HIDDEN UNITS --> between the number of input nodes and number of output nodes 

it.max <- 10000 # Need to be 10'000 to be MONTE-CARLO Permutation

# RANDOM FOREST PARAMETER
## number of tree --> RF and BAG, GBM model (100 good compromise)
################
nb.tree <- 100

overfit="YES" # additionnal restriction used for RF and CSRF to avoid overfitting
alpha.lim = 0.1
##### .....

####################################################################################################################################
####################################################################################################################################
## ADD IN FUNCTION
#########################################################
#######################################################################################################
## lm paramter table
####################

asses.lm. <- function (x,y)
          {
            lin.corr <- lm(y ~ x) #Matrix ."slope","Int","R2","RMSE"
            ## if only one parameter
            if (nrow(summary(lin.corr)$coefficients)==1){
              if (row.names(summary(lin.corr)$coefficients)=="(Intercept)")
                  {
                    slope <- NA #slope
                    std.slope <- NA
                    int. <- summary(lin.corr)$coefficients[1,1] #Intercept
                    std.int <- summary(lin.corr)$coefficients[1,2]
                  }else{
                    slope <- summary(lin.corr)$coefficients[1,1] #slope
                    std.slope <- summary(lin.corr)$coefficients[1,2]
                    int. <- NA #Intercept
                    std.int <- NA
                  }
              
            }else{
            slope <- summary(lin.corr)$coefficients[2,1] #slope
            std.slope <- summary(lin.corr)$coefficients[2,2]
            int. <- summary(lin.corr)$coefficients[1,1] #Intercept
            std.int <- summary(lin.corr)$coefficients[1,2]
            }
            rsquare <- summary(lin.corr)$r.squared #r squared
            MAE <- mean( abs((y - x)),na.rm = TRUE )# mean abs error
            MSE <- mean((y - x)^2,na.rm = TRUE) # mean square error
            RMSE <- sqrt(mean((y - x)^2,na.rm = TRUE)) #root mean square error
            
            lm.vector <- c(slope=slope, std.slope=std.slope, int.=int., std.int=std.int, 
                           rsquare=rsquare, MAE=MAE, MSE=MSE, RMSE=RMSE)
            
            return(lm.vector) 
          }  
  
#################################################################################################################################################
# FUNCTION check and produced subDir folder
###########################################
#February 2017

# modified on the 18/8/2022 to 
creat.subDir <- function (mainDir,subDir)
{
  if ( dir.exists(paste(mainDir,"/",subDir, sep="") ) ){
    
    i <- 1
    while( file.exists( paste(mainDir,"/",subDir,"_",i, sep="") ) )
    {i <-i+1}
    
    dir.create(file.path(mainDir, paste(subDir,"_",i, sep="") ))
    outpath <- file.path(mainDir, paste(subDir,"_",i, sep=""))
    
  } else {
    dir.create(file.path(mainDir, subDir))
    outpath <- file.path(mainDir, subDir)
  }
  
  return(outpath)
}
################################################################################
## TRANSFORMATION list
#######################
## Follow recommendation of Velleman & Hoaglin (1981)
############
transf <- list( alist(x=,(x)), # 1 
                alist(x=,log10((abs(x)))), # 2
                alist(x=,log10((abs(x)+1))), # 3 cool if a lot of data with zero
                alist(x=,(x^2)), # 4
                alist(x=,sqrt(abs((x)))), # 5
                # inverse also good to test
                alist(x=,(1/x)), # 6
                alist(x=,(1/log10(abs(x)))), # 7
                alist(x=,(1/x^2)), # 8
                alist(x=,(1/sqrt(abs(x))))# 9
) ; length(transf)

inv.transf <- list( alist(x=,(x) ), # 1 
                 alist(x=,(10^x) ), # 2
                 alist(x=,(10^x-1) ), # 3 cool if a lot of data with zero
                 alist(x=,(sqrt(x)) ), # 4
                 alist(x=,(x^2) ), # 5
                 # inverse also good to test
                 alist(x=,(1/x) ), # 6
                 alist(x=,(1/10^x) ), # 7
                 alist(x=,(1/sqrt(x) ) ), # 8
                 alist(x=,(1/x^2) )# 9
                ) 

#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################
##########   ---- CROSS VALIDATION FUNCTIONS   ##################################################################################################
#################################################################################################################################################
## All fucntion are adapted from Ecospat package for continuous data
## All data with same weights are considered

## Historic
## --------
## CV-NNET Droz. B - 12.6.2015 -modified 7.9.2015 / 10.2.2017
## CV-RSNNS Droz. B - February 2017
## CV-ELM Droz. B - February 2017
## CV- RF Droz. B - 12.6.2015
## CV - CSRF Droz B. 17.3.2017
## CV- BAG Droz. B - February 2017
## CV- GBM Droz. B - February 2017

## CV- GLM Droz. B - 12.6.2015 --> work for GAM too
## CV- GLM step Droz. B - 12.6.2015
## 

#################################################################################################################################################
## CV-NNET ##
CV.NNET. <- function(data.cv, nnet.unit, nnet.WT=0.01, it.max =1000, K=10, cv.lim = 10, name.sp)
{
  
  n <- nrow(data.cv)
  
  # FORCE AS continuous row name as 1.2.3....
  rownames(data.cv) <- seq(from=1, to= n, by=1 ) 
  
  id <- as.vector(row.names(data.cv), mode = "numeric")
  
  K.lst <- K
  
  cat("K has been finally set to",K.lst, "\n",append = F)
  
  f <- ceiling(n/K.lst)
  s <- sample(rep(1:K.lst,f), n)
  
  # RESPONSE PREDICTION
  
  for (i in 1:K.lst)
  {
    j.out <- id[(s == i)]
    j.out <- sort(j.out)
    j.in <- id[(s != i)]
    j.in <- sort(j.in)
    
    data.cal.cv <- data.cv[j.in,]
    data.test.cv<- data.cv[j.out,]
    
    x.cal <- data.cal.cv[,colnames(data.cal.cv)!=name.sp]
    y.cal <- data.cal.cv[,colnames(data.cal.cv)==name.sp]
    
    x.test <- data.test.cv[,colnames(data.test.cv)!=name.sp]
    
    nnet.cal <-  nnet(x.cal,y.cal, data=data.cal.cv, decay = nnet.WT, size = nnet.unit, 
                      linout = TRUE, maxit = it.max, Hess = TRUE, trace = FALSE)  
    
    nnet.val <- predict(nnet.cal, newdata = x.test , type = "raw",na.rm = TRUE)
    
    if (i == 1)
    {
      vect.id <- j.out
      vect.predicted <- as.vector(nnet.val)
    } else if (i > 1)
    {  
      vect.id <- append(vect.id, j.out, after=length(vect.id))
      vect.predicted <- append(vect.predicted, as.vector(nnet.val), after=length(vect.predicted))
    }  
  }     
  
  df.tmp.res <-  data.frame( cbind( id=round(as.numeric(vect.id),digit=0),
                                    predictions=as.numeric(vect.predicted)) [order(vect.id),])
  
  df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(data.cv[,name.sp]),predictions=df.tmp.res[,2])
  
   return(df.res)    
  }
################################################################################################################################################# 
## CV-RSNNS ##
CV.RSNNS. <- function(data.cv, nral.hidden= round(ncol(data.cv)*2/3) , it.max =1000, K=10, cv.lim = 10, name.sp)
{
  n <- nrow(data.cv)
  
  # FORCE AS continuous row name as 1.2.3....
  rownames(data.cv) <- seq(from=1, to= n, by=1 ) 
  
  id <- as.vector(row.names(data.cv), mode = "numeric")
  
  K.lst <- K
  
  cat("K has been finally set to",K.lst, "\n",append = F)
  
  f <- ceiling(n/K.lst)
  s <- sample(rep(1:K.lst,f), n)
  
  # RESPONSE PREDICTION
  
  for (i in 1:K.lst)
  {
    j.out <- id[(s == i)]
    j.out <- sort(j.out)
    j.in <- id[(s != i)]
    j.in <- sort(j.in)
    
    data.cal.cv <- data.cv[j.in,]
    data.test.cv<- data.cv[j.out,]
    
    x.cal <- data.cal.cv[,colnames(data.cal.cv)!=name.sp]
    y.cal <- data.cal.cv[,colnames(data.cal.cv)==name.sp]
    
    x.test <- data.test.cv[,colnames(data.test.cv)!=name.sp]
    
    nnet.cal <- mlp(x.cal,y.cal, size= nral.hidden, maxit= it.max, learnFunc = "Std_Backpropagation")
    
    nnet.val <- predict(nnet.cal,x.test)
    
    if (i == 1)
    {
      vect.id <- j.out
      vect.predicted <- as.vector(nnet.val)
    } else if (i > 1)
    {  
      vect.id <- append(vect.id, j.out, after=length(vect.id))
      vect.predicted <- append(vect.predicted, as.vector(nnet.val), after=length(vect.predicted))
    }  
  }     
  
  df.tmp.res <-  data.frame( cbind( id=round(as.numeric(vect.id),digit=0),
                                    predictions=as.numeric(vect.predicted)) [order(vect.id),])
  
  df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(data.cv[,name.sp]),predictions=df.tmp.res[,2])
  
  return(df.res)    
}
#################################################################################################################################################
# CV-ELM for elmNNRcpp package 1.04

CV.ELM. <- function(data.cv, elm.obj, K=10, cv.lim = 10, name.sp)
{
  n <- nrow(data.cv)
  
  # FORCE AS continuous row name as 1.2.3....
  rownames(data.cv) <- seq(from=1, to= n, by=1 ) 
  
  id <- as.vector(row.names(data.cv), mode = "numeric")
  
  K.lst <- K
  
  cat("K has been finally set to",K.lst, "\n",append = F)
  
  f <- ceiling(n/K.lst)
  s <- sample(rep(1:K.lst,f), n)
  
  # RESPONSE PREDICTION
  for (i in 1:K.lst)
  {
    j.out <- id[(s == i)]
    j.out <- sort(j.out)
    j.in <- id[(s != i)]
    j.in <- sort(j.in)
    
    data.cal.cv <- data.cv[j.in,]
    data.test.cv<- data.cv[j.out,]
    
    x.cal <- data.cal.cv[,colnames(data.cal.cv)!=name.sp]
    y.cal <- data.cal.cv[,colnames(data.cal.cv)==name.sp]
    
    x.test <- data.test.cv[,colnames(data.test.cv)!=name.sp]
    
    elm.cal <- elm_train(x=as.matrix(x.cal), y=as.matrix(y.cal),
                      nhid=elm.obj$nhid, actfun= elm.obj$actfun  )
    
    elm.val <- elm_predict(elm.cal,as.matrix(x.test))
    
    if (i == 1)
    {
      vect.id <- j.out
      vect.predicted <- as.vector(elm.val)
    } else if (i > 1)
    {  
      vect.id <- append(vect.id, j.out, after=length(vect.id))
      vect.predicted <- append(vect.predicted, as.vector(elm.val), after=length(vect.predicted))
    }  
  }     
  
  df.tmp.res <-  data.frame( cbind( id=round(as.numeric(vect.id),digit=0),
                                    predictions=as.numeric(vect.predicted)) [order(vect.id),])
  
  df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(data.cv[,name.sp]),predictions=df.tmp.res[,2])
  
  return(df.res)    
}

#################################################################################################################################################
## CV-RF ##
CV.RF. <- function(RF.obj, data.cv, K=10, cv.lim = 10, name.sp)
{
  
  n <- nrow(data.cv)
  
  # FORCE AS continuous row name as 1.2.3....
  rownames(data.cv) <- seq(from=1, to= n, by=1 )
  
  # CONTROL IF THE NUMBER OF OBSERVATIONS IS SUFFICIENT
  
  if ((K > n) || (K <= 1))
    stop("K outside allowable range")
  
  id <- as.vector(row.names(data.cv), mode = "numeric")
  K <- K
  K.lst <- round(K)
  kvals <- unique(round(n/(1:floor(n/2))))
  temp <- abs(kvals - K.lst)
  if (!any(temp == 0))
  {
    K.lst <- kvals[temp == min(temp)][1]
  }
  if (K.lst != K)
  {
    cat("K has been set to", K.lst, "\n",append = F)
  }
  
  cat("K has been finally set to",K.lst, "\n",append = F)
  K.lst <- K.lst
  
  f <- ceiling(n/K.lst)
  s <- sample(rep(1:K.lst,f), n)
  
  # RESPONSE PREDICTION
  
  for (i in 1:K.lst)
  {
    j.out <- id[(s == i)]
    j.out <- sort(j.out)
    j.in <- id[(s != i)]
    j.in <- sort(j.in)
    
    data.cal.cv <- data.cv[j.in,]
    data.test.cv<- data.cv[j.out,]
    
    rf.cal <- update(RF.obj, data = data.cal.cv, weights=rep(1,nrow(data.cal.cv)))
    
    rf.val <- predict(rf.cal, data.test.cv , type = "response",weights=rep(1,nrow(data.test.cv)))
    
    if (i == 1)
    {
      vect.id <- j.out
      vect.predicted <- as.vector(rf.val)
    } else if (i > 1)
    {  
      vect.id <- append(vect.id, j.out, after=length(vect.id))
      vect.predicted <- append(vect.predicted, as.vector(rf.val), after=length(vect.predicted))
    }  
  }     
  
  df.tmp.res <-  data.frame( cbind( id=round(as.numeric(vect.id),digit=0) ,predictions=as.numeric(vect.predicted)) [order(vect.id),])
  
  df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(data.cv[,name.sp]),predictions=df.tmp.res[,2])
  
  return(df.res)    
}

################################################################################################################################################
## CV-CSRF ##
CV.CSRF. <- function(CSRF.obj, data.cv, K=10, cv.lim = 10, name.sp)
{
  n <- nrow(data.cv)
  
  # FORCE AS continuous row name as 1.2.3....
  rownames(data.cv) <- seq(from=1, to= n, by=1 )
  
  # CONTROL IF THE NUMBER OF OBSERVATIONS IS SUFFICIENT
  
  if ((K > n) || (K <= 1))
    stop("K outside allowable range")
  
  id <- as.vector(row.names(data.cv), mode = "numeric")
  K <- K
  K.lst <- round(K)
  kvals <- unique(round(n/(1:floor(n/2))))
  temp <- abs(kvals - K.lst)
  if (!any(temp == 0))
  {
    K.lst <- kvals[temp == min(temp)][1]
  }
  if (K.lst != K)
  {
    cat("K has been set to", K.lst, "\n",append = F)
  }
  
  cat("K has been finally set to",K.lst, "\n",append = F)
  K.lst <- K.lst
  
  f <- ceiling(n/K.lst)
  s <- sample(rep(1:K.lst,f), n)
  
  # RESPONSE PREDICTION
  for (i in 1:K.lst)
  {
    j.out <- id[(s == i)]
    j.out <- sort(j.out)
    j.in <- id[(s != i)]
    j.in <- sort(j.in)
    
    data.cal.cv <- data.cv[j.in,]
    data.test.cv<- data.cv[j.out,]
    
    csrf.cal <- update(CSRF.obj, data = data.cal.cv) 
    csrf.val <- predict(csrf.cal, data.test.cv , type = "response") 
    
    if (i == 1)
    {
      vect.id <- j.out
      vect.predicted <- as.vector(csrf.val)$predictions
    } else if (i > 1)
    {  
      vect.id <- append(vect.id, j.out, after=length(vect.id))
      vect.predicted <- append(vect.predicted, as.vector(csrf.val)$predictions, after=length(vect.predicted))
    }  
  }     
  
  df.tmp.res <-  data.frame( cbind( id=round(as.numeric(vect.id),digit=0) , predictions=as.numeric(vect.predicted)) [order(vect.id),])
  
  df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(data.cv[,name.sp]),predictions=df.tmp.res[,2])
  
  return(df.res)    
}

#################################################################################################################################################
## CV-BAG ##
CV.BAG. <- function(data.cv, K=10, cv.lim = 10, nb.tree= 100, name.sp)
{
  n <- nrow(data.cv)
  
  # FORCE AS continuous row name as 1.2.3....
  rownames(data.cv) <- seq(from=1, to= n, by=1 )
  
  # CONTROL IF THE NUMBER OF OBSERVATIONS IS SUFFICIENT
  
  if ((K > n) || (K <= 1))
    stop("K outside allowable range")
  
  id <- as.vector(row.names(data.cv), mode = "numeric")
  K <- K
  K.lst <- round(K)
  kvals <- unique(round(n/(1:floor(n/2))))
  temp <- abs(kvals - K.lst)
  if (!any(temp == 0))
  {
    K.lst <- kvals[temp == min(temp)][1]
  }
  if (K.lst != K)
  {
    cat("K has been set to", K.lst, "\n",append = F)
  }
  
  cat("K has been finally set to",K.lst, "\n",append = F)
  K.lst <- K.lst
  
  f <- ceiling(n/K.lst)
  s <- sample(rep(1:K.lst,f), n)
  
  # RESPONSE PREDICTION
  
  for (i in 1:K.lst)
  {
    j.out <- id[(s == i)]
    j.out <- sort(j.out)
    j.in <- id[(s != i)]
    j.in <- sort(j.in)
    
    data.cal.cv <- data.cv[j.in,]
    data.test.cv<- data.cv[j.out,]
    
    x.cal <- data.cal.cv[,colnames(data.cal.cv)!=name.sp]
    y.cal <- data.cal.cv[,colnames(data.cal.cv)==name.sp]
    
    bag.cal <- ipredbagg (y.cal , x.cal, nbag=nb.tree)
    
    bag.val <- predict(bag.cal, newdata = data.test.cv , type = "raw" , na.rm = TRUE)
    
    if (i == 1)
    {
      vect.id <- j.out
      vect.predicted <- as.vector(bag.val)
    } else if (i > 1)
    {  
      vect.id <- append(vect.id, j.out, after=length(vect.id))
      vect.predicted <- append(vect.predicted, as.vector(bag.val), after=length(vect.predicted))
    }  
  }     
  
  df.tmp.res <-  data.frame( cbind( id=round(as.numeric(vect.id),digit=0) ,predictions=as.numeric(vect.predicted)) [order(vect.id),])
  
  df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(data.cv[,name.sp]),predictions=df.tmp.res[,2])
  
  return(df.res)    
}

#################################################################################################################################################
## CV-EGB ##
CV.EGB. <- function(data.cv, K=10, cv.lim = 10, it.max= 1000, name.sp)
{
  n <- nrow(data.cv)
  
  # FORCE AS continuous row name as 1.2.3....
  rownames(data.cv) <- seq(from=1, to= n, by=1 )
  
  # CONTROL IF THE NUMBER OF OBSERVATIONS IS SUFFICIENT
  
  if ((K > n) || (K <= 1))
    stop("K outside allowable range")
  
  id <- as.vector(row.names(data.cv), mode = "numeric")
  K <- K
  K.lst <- round(K)
  kvals <- unique(round(n/(1:floor(n/2))))
  temp <- abs(kvals - K.lst)
  if (!any(temp == 0))
  {
    K.lst <- kvals[temp == min(temp)][1]
  }
  if (K.lst != K)
  {
    cat("K has been set to", K.lst, "\n",append = F)
  }
  
  cat("K has been finally set to",K.lst, "\n",append = F)
  K.lst <- K.lst
  
  f <- ceiling(n/K.lst)
  s <- sample(rep(1:K.lst,f), n)
  
  # RESPONSE PREDICTION
  for (i in 1:K.lst)
  {
    j.out <- id[(s == i)]
    j.out <- sort(j.out)
    j.in <- id[(s != i)]
    j.in <- sort(j.in)
    
    # CREAT DATA SET FOR x-FOLD CV
    data.cal.cv <- data.cv[j.in,]
    data.test.cv<- data.cv[j.out,]
    
    x.cal <- data.cal.cv[,colnames(data.cal.cv)!=name.sp]
    y.cal <- data.cal.cv[,colnames(data.cal.cv)==name.sp]
    
    x.test <- data.test.cv[,colnames(data.test.cv)!=name.sp]
    
    ## -- MODEL --
    # set the parameter of the booster  (default)
    para <- list(booster= "gbtree", max_depth = 6, eta = 0.3, silent = 1, nthread = 2,
                 objective = "reg:linear", eval_metric = "rmse")
    
    # run model
    egb.cal <- xgboost(data = as.matrix(x.cal), label = t(as.vector(y.cal)), missing = NA, weight = NULL,
                         params = para , nrounds = it.max, verbose = 0 )
    
    egb.val <-  predict(egb.cal, as.matrix(x.test)) 
    
    if (i == 1)
    {
      vect.id <- j.out
      vect.predicted <- egb.val
    } else if (i > 1)
    {  
      vect.id <- append(vect.id, j.out, after=length(vect.id))
      vect.predicted <- append(vect.predicted, egb.val, after=length(vect.predicted))
    }  
  }     
  
  df.tmp.res <-  data.frame( cbind( id=round(as.numeric(vect.id),digit=0) , predictions=as.numeric(vect.predicted)) [order(vect.id),])
  
  df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(data.cv[,name.sp]),predictions=df.tmp.res[,2])
  
  return(df.res)    
}
#################################################################################################################################################
## CV-GBM ##
CV.GBM. <- function(gbm.obj, data.cv, nb.tree= 100, K=10, cv.lim = 10, name.sp)
{
  data.cv <- as.data.frame(data.cv)
  
  n <- nrow(data.cv)
  
  # FORCE AS continuous row name as 1.2.3....
  rownames(data.cv) <- seq(from=1, to= n, by=1 ) 
  
  id <- as.vector(row.names(data.cv), mode = "numeric")
  
  K.lst <- K
  
  cat("K has been finally set to",K.lst, "\n",append = F)
  
  f <- ceiling(n/K.lst)
  s <- sample(rep(1:K.lst,f), n)
  
  # RESPONSE PREDICTION
  
  for (i in 1:K.lst)
  {
    j.out <- id[(s == i)]
    j.out <- sort(j.out)
    j.in <- id[(s != i)]
    j.in <- sort(j.in)
    
    data.cal.cv <- data.cv[j.in,]
    data.test.cv<- data.cv[j.out,]
    
    x.cal <- data.cal.cv[,colnames(data.cal.cv)!=name.sp]
    y.cal <- data.cal.cv[,colnames(data.cal.cv)==name.sp]
    
    x.test <- data.test.cv[,colnames(data.test.cv)!=name.sp]
    
    gbm.cal <- update(gbm.obj, data = data.cal.cv, 
                      distribution = "gaussian") 
    
    gbm.val <- predict.gbm(gbm.cal, newdata=as.data.frame(x.test), 
                           n.trees=nb.tree, type="response",na.rm = TRUE)
    
    if (i == 1)
    {
      vect.id <- j.out
      vect.predicted <- as.vector(gbm.val)
    } else if (i > 1)
    {  
      vect.id <- append(vect.id, j.out, after=length(vect.id))
      vect.predicted <- append(vect.predicted, as.vector(gbm.val), after=length(vect.predicted))
    }  
  }     
  
  df.tmp.res <-  data.frame( cbind( id=round(as.numeric(vect.id),digit=0),
                                    predictions=as.numeric(vect.predicted)) [order(vect.id),])
  
  df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(data.cv[,name.sp]),predictions=df.tmp.res[,2])
  
  return(df.res)    
}

###############################################################################################################
## CV-GLM   ##
CV.glm. <- function(glm.obj, K=10, cv.lim = 10, name.sp)
{
  
  data.cv <<- glm.obj$data
  n <- nrow(data.cv)
  
  # FORCE AS continuous row name as 1.2.3....
  rownames(data.cv) <- seq(from=1, to= n, by=1 )
  
  # CONTROL IF THE NUMBER OF OBSERVATIONS IS SUFFICIENT
  
  if ((K > n) || (K <= 1))
    stop("K outside allowable range")
  
  id <- as.vector(row.names(data.cv), mode = "numeric")
  K <- K
  K.lst <- round(K)
  kvals <- unique(round(n/(1:floor(n/2))))
  temp <- abs(kvals - K.lst)
  if (!any(temp == 0))
  {
    K.lst <- kvals[temp == min(temp)][1]
  }
  if (K.lst != K)
  {
    cat("K has been set to", K.lst, "\n",append = F)
  }
  
  cat("K has been finally set to",K.lst, "\n",append = F)
  K.lst <- K.lst
  
  f <- ceiling(n/K.lst)
  s <- sample(rep(1:K.lst,f), n)
  
  # RESPONSE PREDICTION
  
  for (i in 1:K.lst)
  {
    j.out <- id[(s == i)]
    j.out <- sort(j.out)
    j.in <- id[(s != i)]
    j.in <- sort(j.in)
    
    data.cal.cv <- data.cv[j.in,]
    data.test.cv<- data.cv[j.out,]
    
    glm.cal <- update(glm.obj, data = data.cal.cv, weights=rep(1,nrow(data.cal.cv)))
    
    glm.val <- predict(glm.cal, data.test.cv , type = "response",weights=rep(1,nrow(data.test.cv)))
    
    if (i == 1)
    {
      vect.id <- j.out
      vect.predicted <- as.vector(glm.val)
    } else if (i > 1)
    {  
      vect.id <- append(vect.id, j.out, after=length(vect.id))
      vect.predicted <- append(vect.predicted, as.vector(glm.val), after=length(vect.predicted))
    }	
  }     
  
  df.tmp.res <-  data.frame( cbind( id=round(as.numeric(vect.id),digit=0) ,predictions=as.numeric(vect.predicted)) [order(vect.id),])
  
  df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(data.cv[,name.sp]),predictions=df.tmp.res[,2])
  
  return(df.res)    
}

###############################################################################################################
## CV-GLM stepwise
CV.glm.step <- function(glm.obj, K=10, cv.lim = 10, name.sp)
{
  
  data.cv <<- glm.obj$data
  n <- nrow(data.cv)
  
  # FORCE AS continuous row name as 1.2.3....
  rownames(data.cv) <- seq(from=1, to= n, by=1 )
  
  # CONTROL IF THE NUMBER OF OBSERVATIONS IS SUFFICIENT
  
  if ((K > n) || (K <= 1))
    stop("K outside allowable range")
  
  id <- as.vector(row.names(data.cv), mode = "numeric")
  K <- K
  K.lst <- round(K)
  kvals <- unique(round(n/(1:floor(n/2))))
  temp <- abs(kvals - K.lst)
  if (!any(temp == 0))
  {
    K.lst <- kvals[temp == min(temp)][1]
  }
  if (K.lst != K)
  {
    cat("K has been set to", K.lst, "\n",append = F)
  }
  
  cat("K has been finally set to",K.lst, "\n",append = F)
  K.lst <- K.lst
  
  f <- ceiling(n/K.lst)
  s <- sample(rep(1:K.lst,f), n)
  
  # RESPONSE PREDICTION
  
  for (i in 1:K.lst)
  {
    j.out <- id[(s == i)]
    j.out <- sort(j.out)
    j.in <- id[(s != i)]
    j.in <- sort(j.in)
    
    data.cal.cv <- data.cv[j.in,]
    data.test.cv<- data.cv[j.out,]
    
    #w.cal <- rep(1,nrow(data.cal.cv))
    #w.test <- rep(1,nrow(data.test.cv))
    
    glm.cal <- update(glm.obj, data = data.cal.cv ) #, weights= w.cal)
    
    glm.val <- predict(glm.cal, data.test.cv , type = "response") #,weights= w.test)
    
    if (i == 1)
    {
      vect.id <- j.out
      vect.predicted <- as.vector(glm.val)
    } else if (i > 1)
    {  
      vect.id <- append(vect.id, j.out, after=length(vect.id))
      vect.predicted <- append(vect.predicted, as.vector(glm.val), after=length(vect.predicted))
    }  
  }     
  
  df.tmp.res <-  data.frame( cbind( id=round(as.numeric(vect.id),digit=0) ,predictions=as.numeric(vect.predicted)) [order(vect.id),])
  
  df.res <- data.frame(id=df.tmp.res[,1],obs=as.vector(data.cv[,name.sp]),predictions=df.tmp.res[,2])
  
  return(df.res)    
}

#################################################################################################################################################
#################################################################################################################################################
#################################################################################################################################################
##########   ---- PERMUTATION VARIABLE FUNCTIONS   ##################################################################################################
##             --> calculate the relative importance of each variable !!!!!
#################################################################################################################################################

## Historic
## --------
## Original function write by C.Randin modified by B.Droz!
## var.imp.rf Droz. B - 12.6.2015
## var.imp.glm.step Droz. B - 12.6.2015 --> work for glm too

## Input :
## =======
### model :  model object
### cal : calibration dataset used to build the model
### names.pred: name of the predicting variables
### nperm : number of time each variable is permutated

####################################################################################################################################3333333
#########################
## VarImp RandomForest ##
#########################
var.imp.rf <- function(model,cal,names.pred,nperm=100)
{
	ref<-as.numeric(as.character(predict(model,cal)))
	VarImp<-vector()
	
	for (i in 1:length(names.pred))
	{
		print(names.pred[i])
		refi<-vector()
		
		for (j in 1:nperm)
		{
			if (j%%100==0)
			{
				cat("> Permutation ",j, "\n",append = FALSE)
			}
			
			cali<-cal
			cali[,names.pred[i]]<-cali[sample(1:nrow(cali),nrow(cali)),names.pred[i]]

			refi<-c(refi,1-cor(ref,as.numeric(as.character(predict(model,cali)))))
		}
		
		VarImp<-c(VarImp,round(mean(refi),3))
	}
	
	names(VarImp)<-names.pred
	return<-VarImp
}

###############################################################################
#########################
## VarImp elmNNRcpp - v1.0.4 ##
#########################
# Droz - August 2022

var.imp.elm <- function(model,cal,names.pred,nperm=100)
{
  require(elmNNRcpp)
  
  ref<-as.numeric(as.character(elm_predict(model,as.matrix(cal) )) )
  VarImp<-vector()
  
  for (i in 1:length(names.pred))
  {
    print(names.pred[i])
    refi<-vector()
    
    for (j in 1:nperm)
    {
      if (j%%100==0)
      {
        cat("> Permutation ",j, "\n",append = FALSE)
      }
      
      cali<-cal
      cali[,names.pred[i]]<-cali[sample(1:nrow(cali),nrow(cali)),names.pred[i]]
      
      refi<-c(refi,1-cor(ref,as.numeric(as.character(elm_predict(model,as.matrix(cali) )))) )
    }
    
    VarImp<-c(VarImp,round(mean(refi),3))
  }
  
  names(VarImp)<-names.pred
  return<-VarImp
}

## VarImp GLM, GLM Step, GAM
#############################
var.imp.glm.step <- function(model,cal,names.pred,nperm=100)
{
	
	ref<-predict(model,cal)
	VarImp<-vector()
	
	for (i in 1:length(names.pred))
	{
		print(names.pred[i])
		refi<-vector()
		
		for (j in 1:nperm)
		{
			if (j%%100==0)
			{
				cat("> Permutation ",j, "\n",append = F)
			}
			
			cali<-cal
			cali[,names.pred[i]]<-cali[sample(1:nrow(cali),nrow(cali)),names.pred[i]]
			refi<-c(refi,1-cor(ref,predict(model,cali)))
		}
		
		VarImp<-c(VarImp,round(mean(refi),3))
	}
	
	names(VarImp)<-names.pred
	return<-VarImp
}

####################################################################################################################################
####################################################################################################################################
## SCRIPT START HERE
####################################################################################################################################
# Set folder --> 
date <- Sys.Date()
inpath <- paste(workdir,"/input/", sep="")
folder <- paste("/",date,"_mod.proj", sep="")
outpath <- creat.subDir(paste(workdir,"/output",sep=""), folder)

## PARRALLEL CORE SETUP# 
##--------------------
# ptm <- proc.time()# ignite timer

# beginCluster(detectCores()) # ACTIVATE THIS MULTI CORE CALCULATION 
###############################################################################
### write the option of the run
################################
f.info <- paste(outpath,"/AA_INFO_PARA.txt",sep="")

cat( paste("*** BUILD AND RUN MODEL PREDICTION ---", Sys.Date()), file= f.info, sep="\n")
cat("###########################",file= f.info,append=TRUE, sep="\n")
cat("R-script model.proj_v14 - Droz 2022",file= f.info,append=TRUE, sep="\n")
cat(paste("RUN for the data set:", f.dat),file= f.info,append=TRUE, sep="\n")
cat(paste("name of variable selection do you use?", var.sele),file= f.info,append=TRUE, sep="\n")
cat(paste("remove sample list:", f.rem),file= f.info,append=TRUE, sep="\n")
cat(paste("number of replication?", n.rep),file= f.info,append=TRUE, sep="\n")
cat(paste("percent of sampling used for model cali:", sampl),file= f.info,append=TRUE, sep="\n")
cat(paste("Are variable transform?", TRANSF),file= f.info,append=TRUE, sep="\n")
cat(paste("Do you have rescaling?", RES),file= f.info,append=TRUE, sep="\n")
cat(paste("Do you use support vector classification ?", SVC),file= f.info,append=TRUE, sep="\n")

cat("###########################",file= f.info,append=TRUE, sep="\n")
cat( "MODEL TECHNIC USED" ,file= f.info,append=TRUE, sep="\n")
cat("###########################",file= f.info,append=TRUE, sep="\n")
cat(paste("Neural Networks?", NNET),file= f.info,append=TRUE, sep="\n")
cat(paste("Stuttgart Neural Network Simulator?", RSNNS),file= f.info,append=TRUE, sep="\n")
cat(paste("Extrem learning machine?", ELM),file= f.info,append=TRUE, sep="\n")
cat(paste("Random Forest?", RF),file= f.info,append=TRUE, sep="\n")
cat(paste("Case-specific Random Forest?", CSRF),file= f.info,append=TRUE, sep="\n")
cat(paste("Bagging tree method?", BAG),file= f.info,append=TRUE, sep="\n")
cat(paste("Extreme Gradient Boosting?", EGB),file= f.info,append=TRUE, sep="\n")
cat(paste("Generalized Boosted Regression Models?", GBM),file= f.info,append=TRUE, sep="\n")
cat(paste("Generalized Linear Models?", GLM),file= f.info,append=TRUE, sep="\n")
cat(paste("GLM step-wise?", GLM_sw),file= f.info,append=TRUE, sep="\n")
cat(paste("Generalized Additive Models?", GAM),file= f.info,append=TRUE, sep="\n")

cat("###########################",file= f.info,append=TRUE, sep="\n")
cat( "FUNCTION PARAMETER" ,file= f.info,append=TRUE, sep="\n")
cat("###########################",file= f.info,append=TRUE, sep="\n")            
cat(paste("numb of permutation: ", n.perm),file= f.info,append=TRUE, sep="\n")
cat(paste("X fold cross validation: ", fold.cv),file= f.info,append=TRUE, sep="\n")
cat(paste("NNET para - weight decay: ", deca),file= f.info,append=TRUE, sep="\n")
cat(paste("NNET para - number max of hidden units: ", UNIT.Max),file= f.info,append=TRUE, sep="\n")
cat(paste("NNET para - iteration max: ", it.max),file= f.info,append=TRUE, sep="\n")
cat(paste("RF, BAG, GBM para - numb tree: ", nb.tree),file= f.info,append=TRUE, sep="\n")
cat(paste("RF, CSRF - numb tree: ", nb.tree),file= f.info,append=TRUE, sep="\n")
cat(paste("RF, CSRF - overfit restr: ", overfit),file= f.info,append=TRUE, sep="\n")
cat(paste("RF, CSRF - numb tree: ", alpha.lim),file= f.info,append=TRUE, sep="\n")

###################################################################################################################################
# -- LOADING DATA SET -----
##############################################
# Open the data set
data.in <- read.table(paste(inpath,f.dat,sep=""),header = TRUE ,
                      na.strings = "NA", sep=",")

var.list <- read.table(paste(inpath,f.var,sep=""), header= TRUE)
rem.list <- read.table(paste(inpath,f.rem,sep=""), header= FALSE)

##########################################################################################################
## ---  BUILD MODEL FOR DIVERS DEPENDENT VARIABLE ---
#######################################################
var.y <- na.omit(var.list)
var.y <- var.y$var.name[var.y$var.type=="y"]
outpath ->  outpath.0  # keep the originate outpath
data.in -> data.0 ## keep the originate dataset
var.list -> var.list.0# keep the originate list
var.sele -> var.sele.0 ## keep original set list

# temporary value to compare data at the end
# --- AVERAGE TABLE ---
mod.name.perf <- NULL
modl.perf.av <- NULL
sv.perf.av <- NULL

for (tu in 1:length(var.y))# selected variable y to pred
{
  # creat a new folder for the considered y (dep.var)
  output <- paste(outpath.0,"/",var.y[tu],sep="")
  creat.subDir (outpath.0,var.y[tu])
  
  c.p. <- var.list.0[var.list.0$var.name==var.y[tu],
                   names(var.list.0) %in% var.sele.0] # column
  
  var.sele <- names(c.p.)[c.p.==1]
  
  if (length(var.sele)==0) {}else{# ? does the species exist
    
  for (zu in 1:length(var.sele)) ## list of selected variable
  {
  # creat a new folder for the considered y (dep.var)
  output <- paste(outpath.0,"/",var.y[tu],"/",var.sele[zu],sep="")
  creat.subDir (paste(outpath.0,"/",var.y[tu],sep=""),var.sele[zu]) 
  
  ########################################
  ## Creat Organized data set for the test
  ## --> remove selected sampl ( remove list)
  ## --> select only variable (x and y) who as one in the selected list 
  ## --> transformed to normalized it
  data.0 -> data.in # relaod originate data
  var.list <- var.list.0# reload originte list
  var.list <-na.omit(var.list)
  data.in <- data.in[!(data.in$sample %in% rem.list$V1),] # remove sample
  p. <- names(var.list)== var.sele[zu] # column
  pos <-match( c(var.list$var.name[var.list[,p.]==1 & 
                                  var.list$var.type=="x"],var.y[tu] ), names(data.in))
  data.in <- data.in[,pos]
  var.list <- var.list[var.list[,p.]==1,]
  
  ### SAMPLING CALI -VALI DATASET
  #################################
  # creat list for n repli of data for cali and validate model
  d.cali <- list(NULL)
  d.vali <- list(NULL)

  # creat list of sampling
  for (i in 1:n.rep) # list of sampling data
  {
    n.s <- sample(1:nrow(data.in),nrow(data.in)*sampl/100, replace=FALSE)
    d.cali[[i]] <- data.in[n.s,]
    d.vali[[i]] <- data.in[-n.s,]
    d.cali[[i]] <- na.omit(d.cali[[i]])
    d.vali[[i]] <- na.omit(d.vali[[i]])
    d.cali[[i]] <- d.cali[[i]][is.finite(rowSums(d.cali[[i]]) ),]
    d.vali[[i]] <- d.vali[[i]][is.finite(rowSums(d.vali[[i]]) ),]
  }
  
  ### APPLIED SVC - Support vector classification
  ###############################################
  if (SVC=="YES"){
    
    mod.tech <- "SVC"
    outpath <- paste(output,"/",mod.tech,sep="")
    creat.subDir (output,mod.tech)
    
    scv.cali <- list(NULL)
    scv.vali <- list(NULL)
    acc.svc <- NULL
    
    for (i in 1:n.rep)
    {
      y <- d.cali[[i]][,ncol(d.cali[[i]])] 
      y[y>0] <- 1 # classifier
      
      y.vali <-d.vali[[i]][,ncol(d.vali[[i]])] # data real
      y.vali[y.vali>0] <- 1 # classifier
      
      if (length(y)==sum(y)|length(y.vali)==sum(y.vali)) # if no zero
      {
        scv.cali[[i]] <- y
        scv.vali[[i]] <- y.vali
        
        acc.svc <-rbind(acc.svc,c(n.run=i,acc.cali=NA, acc.vali=NA))
        
      }else{
        d.slc <-  d.cali[[i]][,1:(ncol(d.cali[[i]])-1)] 
        d.slc <- data.frame(d.slc,y=as.factor(y))
        fit.svm <- svm(y ~ .,data=d.slc, kernel=k.nel, cost=10,scale = FALSE)
        #fit.svm$index #gives us the index of the Support Vectors
        scv.cali[[i]] <- fit.svm$fitted #so we have xx support vectors
        #table(Predicted=fit.svm$fitted,Actual=y)
        acc.cali <-mean(fit.svm$fitted==y)*100 #accuracy on Training Set
      
        scv.vali[[i]] <- predict(fit.svm, newdata=d.vali[[i]][,1:(ncol(d.vali[[i]])-1)],
                                 type="raw",na.rm = TRUE)
      
        acc.vali <- mean(scv.vali[[i]]==y.vali,na.rm=TRUE)*100 #accuracy on valid Set
      
        # Save Model - R file
        file.name<-paste(outpath,"/",mod.tech,"_model_",i,sep='')
        save(list=c('fit.svm'),file=file.name)
      
        acc.svc <- rbind(acc.svc,c(n.run=i,acc.cali=acc.cali, acc.vali=acc.vali))
      }
    
      write.table(acc.svc, file=paste(outpath,"/",mod.tech,"_accuracy_x_run.txt",sep=""),
                  sep="\t", append=FALSE, row.names=FALSE,col.names=TRUE, quote=FALSE)
    }
    
      # make the average accurency svc data
      tmp <- data.frame(species= var.y[tu], var.sel= var.sele[zu], technic= mod.tech, n.rep= n.rep,
                        m.acc.cali = mean(acc.svc[,2],na.rm=TRUE), sd.acc.cali = sd(acc.svc[,2],na.rm=TRUE),
                       m.acc.vali = mean(acc.svc[,3],na.rm=TRUE), sd.acc.vali = sd(acc.svc[,3],na.rm=TRUE) )
      
      sv.perf.av<- rbind(sv.perf.av,tmp)
    
  }else{
    y <- d.cali[[i]][,ncol(d.cali[[i]])] 
    y.vali <-d.vali[[i]][,ncol(d.vali[[i]])] # data real
    }
  
  ### NO enough occurence data stop zu loop her !!!!
  ###########################################
  if ( ( 1- sum(c(y.vali,y))/ length(c(y.vali,y)) )>=RA.O )
  {
  }else{
  
  #########################
  ### APPLIED DATA MODIFICATION
  ##########################
  if (TRANSF=="YES"){
    for (n in 1:length(d.cali))
        {
        for (i in 1:(ncol(d.cali[[n]]) ) )
            {
              f <- as.function(transf[[var.list$var.transf
                                       [var.list$var.name==names(d.cali[[n]][i]) ]  ]]) # define function from the list
              d.cali[[n]][,i] <- f(x=d.cali[[n]][,i])# applied modification
              d.vali[[n]][,i] <- f(x=d.vali[[n]][,i])
            }
        }
  }else{}
  
  #########################
  ## RESCALE - z-score
  #########################
  if (RES=="YES"){
    d.all <- rbind(d.cali[[1]],d.vali[[1]])
    d.all <- do.call(data.frame,lapply(d.all, function(x) replace(x, is.infinite(x),NA)))

    d.mean <- apply(d.all,2,FUN=mean,na.rm=TRUE)
    d.sd <- apply(d.all,2,FUN=sd, na.rm=TRUE)
    
    RES.data <- rbind(mean.var=d.mean,sd.var=d.sd)
    
    # write RESCALE table
    write.table(RES.data, file=paste(output,"/RESCALING_values.txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    for (n in 1:length(d.cali))
        {
      
      for (k in 1:ncol(d.cali[[n]]))
        {
          ## RESCALE THE VARIABLE with normal rescaling
          #f.res <- function(x) {(x-mean(x, na.rm=TRUE))/sd(x,na.rm=TRUE)}
          d.cali[[n]][,k] <- (d.cali[[n]][,k]- d.mean[k])/d.sd[k]
          d.vali[[n]][,k] <- (d.vali[[n]][,k]- d.mean[k])/d.sd[k]
        }
      }
    
  } else{}
  
  # remove NA/NAN/Inf and produced clean data set of cali vali
  for (n in 1:length(d.cali))
  {
    d.cali[[n]] <-d.cali[[n]][is.finite(rowSums(d.cali[[n]])),] # remove inf
    d.vali[[n]] <-d.vali[[n]][is.finite(rowSums(d.vali[[n]])),] 
    
    d.cali[[n]] <-na.omit(d.cali[[n]]) # remove na
    d.vali[[n]] <-na.omit(d.vali[[n]])
    
    if (SVC=="YES"){
      scv.cali[[n]] <- scv.cali[[n]][is.finite(rowSums(d.cali[[n]]))]# remove inf
      scv.vali[[n]] <- scv.vali[[n]][is.finite(rowSums(d.vali[[n]]))]
      
      scv.cali[[n]] <- scv.cali[[n]][apply(is.na(d.cali[[n]]),1,sum)==0 ] # remove na
      scv.vali[[n]] <- scv.vali[[n]][apply(is.na(d.vali[[n]]),1,sum)==0 ]
      
    }else{}
      
  }  
#################################################################################################################################
#################################################################################################################################
## initialise the list data

# --- model performance ---
modl.perf.nnet <- NULL
modl.perf.rsnns <- NULL
modl.perf.elm <- NULL
modl.perf.rf <- NULL
modl.perf.csrf <- NULL

modl.perf.bag <- NULL
modl.perf.egb <- NULL
modl.perf.gbm <- NULL

modl.perf.glm <-NULL
modl.perf.glmstep <- NULL
modl.perf.gam <- NULL

# --- model importance ---
varimp.nnet <- NULL
varimp.rsnns <- NULL
varimp.elm <- NULL
varimp.rf <- NULL
varimp.csrf <- NULL

varimp.bag <- NULL
varimp.egb <- NULL
varimp.gbm <- NULL

varimp.glm <- NULL
varimp.glmstep <- NULL
varimp.gam <- NULL

varimp.av <- NULL
mod.name <- NULL

count.while <- 0
pred.st.list <-list(NULL)

#####################################################################################################################################
#####################################################################################################################################      
## -------------------------   MODEL CALIBRATION ------------------------------------------
#####################################################################################################################################   
if (NNET == "NO") {}else{
  
  mod.tech <- "NNET"
  outpath <- paste(output,"/",mod.tech,sep="")
  creat.subDir (output,mod.tech)
  
  for (i in 1:n.rep)    
  {
    # creat x and y table
    #####################
    if (SVC=="YES")
    {
      y.in <- data.frame(d.cali[[i]][,names(d.cali[[i]])==var.y[tu]])
      y.in <- data.frame(y.in[scv.cali[[i]]==1,1])
      p.cal.svc <- as.numeric(row.names(d.cali[[i]]))[scv.cali[[i]]==1]
      names(y.in) <- var.y[tu]
      
      x.in <- d.cali[[i]][, na.omit( match(names(d.cali[[i]]), 
                                           var.list$var.name[var.list$var.type=="x"] ) ) ]
      x.in <- x.in[scv.cali[[i]]==1,]
      data.in <- cbind(y.in,x.in)
      
      x.vali <- d.vali[[i]][, na.omit( match(names(d.vali[[i]]), 
                                             var.list$var.name[var.list$var.type=="x"] ) ) ]
      x.vali <- x.vali[scv.vali[[i]]==1,]
      p.val.svc <- as.numeric(row.names(d.vali[[i]]))[scv.vali[[i]]==1]
      
      y.vali <- data.frame(d.vali[[i]][,names(d.vali[[i]])==var.y[tu]])
      y.vali <- data.frame(y.vali[scv.vali[[i]]==1,1])
      
    }else{ 
      y.in <- data.frame(d.cali[[i]][,names(d.cali[[i]])==var.y[tu]])
      names(y.in) <- var.y[tu]
      x.in <- d.cali[[i]][, na.omit( match(names(d.cali[[i]]), 
                                           var.list$var.name[var.list$var.type=="x"] ) ) ]
      data.in <- cbind(y.in,x.in)
      
      x.vali <- d.vali[[i]][, na.omit( match(names(d.vali[[i]]), 
                                             var.list$var.name[var.list$var.type=="x"] ) ) ]
      y.vali <- data.frame(d.vali[[i]][,names(d.vali[[i]])==var.y[tu]])
    }
    ######################
    ### NETWORK ANALYSIS
    #####################                       
    ## OPTIMIZED PARAMETER PERFORMANCE OF NNET
    ##########################################
    
    if (i ==1) {  cat("> NNET Parameter one layer ...", "\n",append = FALSE)
      
      ##(from the FAQ for a commercial neural network software company)  
      ## Number of inputs + outputs) * (2/3) -- Heaton 2008
      nb.hid1 <- round (ncol(data.in) *2/3)
      
      nb.hid2 <- round ( nb.hid1*2/3 )
      
      nb.hidden =c(nb.hid1,nb.hid2)
      
      WT.opt <- deca
      
    } else {}
    
    cat("> NNET MODEL run ",i,"sampling  ...", "\n",append = FALSE)
    
    net.tr <- nnet(x.in,y.in,data=data.in, decay = WT.opt, size = nb.hid1, 
                  linout = TRUE, maxit = it.max, Hess = TRUE) 
    
    # Save Model - R file
    file.name<-paste(outpath,"/",mod.tech,"_model_",i,sep='')
    save(list=c('net.tr'),file=file.name)
    
    ######################
    ## Variable importance
    ######################
    ## weight algorithm (Garson 1991)
    ##-------------------------------
    g  <- garson(net.tr, colnames(y.in), bar_plot = TRUE, x_lab = NULL,
                 y_lab = "Rel.importance", wts_only = FALSE)
    
    # Re-order data by row names
    order <- as.numeric(row.names(g$data))
    sort <- sort(order, index.return=TRUE)
    g$data[,1] <- g$data[sort$ix,1] 
    g$data[,2] <- g$data[sort$ix,2] 
    
    varimp.nnet <- rbind(varimp.nnet,g$data[,1])
    
    colnames(varimp.nnet) <- g$data[,2]
    
    ## MODEL PERFORMANCE
    ####################
    ## PREDICT THE OBS
    pred.y <- predict(net.tr, newdata=x.in,type="raw",na.rm = TRUE)
    # validate the prediction
    pred.y.vali <-predict(net.tr, newdata=x.vali,type="raw",na.rm = TRUE)
    
    ## recombin SVC classification and ML prediction
    ################################################
    if (SVC=="YES"){ 
      pred.y <- rbind(pred.y, matrix( rep(0,nrow(d.cali[[i]])-nrow(pred.y)) )) 
      pred.y.vali <- rbind(pred.y.vali, matrix( rep(0,nrow(d.vali[[i]])-nrow(pred.y.vali)) ) )
      
      y.in <- c(y.in[,1], d.cali[[i]][ scv.cali[[i]]==0,names(d.cali[[i]])==var.y[tu] ] )
      y.vali <- c(y.vali[,1], d.vali[[i]][ scv.vali[[i]]==0,names(d.vali[[i]])==var.y[tu] ] )
      
      p.cali <- as.numeric(row.names(x.in))
      p.vali <- as.numeric(row.names(x.vali))
      
      if (nrow(d.cali[[i]])==length(p.cal.svc)){}else{ 
        p.cali <- c(p.cali,
                    as.numeric(row.names(d.cali[[i]]))[!as.numeric(row.names(d.cali[[i]]))%in%p.cali]) }
      
      if (nrow(d.vali[[i]])==length(p.val.svc)){}else{ 
        p.vali <- c(p.vali,
                    as.numeric(row.names(d.vali[[i]]))[!as.numeric(row.names(d.vali[[i]]))%in%p.vali])}

    } else{
      y.in <- y.in[,1]
      y.vali <- y.vali[,1]
      
      p.cali <- as.numeric(row.names(x.in))
      p.vali <- as.numeric(row.names(x.vali))
    }
    
    # PERFORMANCE OF THE MODEL
    #	lin.corr <- lm(y.in ~ pred.y) ;summary(lin.corr)
    M.perf <- asses.lm.(pred.y, y.in)
    V.perf <- asses.lm.(pred.y.vali, y.vali)
    
    names(M.perf) <- paste(rep("M_",5),names(M.perf),sep="")
    names(V.perf) <- paste(rep("V_",5),names(V.perf),sep="")
    
    ##keep data and rescale it
    if (RES=="YES"){
      m.temp <- RES.data[,ncol(RES.data)][1]
      sd.temp <- RES.data[,ncol(RES.data)][2]
      
      pred.y <- (pred.y*sd.temp) + m.temp
      y.in <- (y.in*sd.temp) + m.temp
      pred.y.vali <- (pred.y.vali*sd.temp) + m.temp
      y.vali <- (y.vali*sd.temp) + m.temp
      
    } else{}
    
    ## keep data and retransformed it...
    if (TRANSF=="YES"){
      f <- as.function(inv.transf[[var.list$var.transf[var.list$var.name==var.y[tu]]]])
      cali <- data.frame(cbind(p.cali,rep("CALI",length(pred.y)), f(x=y.in), f(x=pred.y) ))
      vali <- data.frame(cbind(p.vali,rep("VALI",length(pred.y.vali)), f(x=y.vali), f(x=pred.y.vali) ))
    }else{
      cali <- data.frame(cbind(p.cali,rep("CALI",length(pred.y)), y.in, pred.y ))
      vali <- data.frame(cbind(p.vali,rep("VALI",length(pred.y.vali)), y.vali, pred.y.vali ))
    }
    
    names(cali) <-c("id","Type","obs.y","pred.y")
    names(vali) <- c("id","Type","obs.y","pred.y") 
    
    cali.vali <- rbind(cali,vali)
    
    write.table(cali.vali, file=paste(outpath,"/",mod.tech,"_cali.vali_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=FALSE,col.names=TRUE, quote=FALSE)
    
    ## CROSS - Validate model
    ########################
    if (CR.VAL=="YES"){
    CV.data  <- CV.NNET.(data.in, nnet.unit= nb.hid1, nnet.WT= WT.opt, it.max= it.max, 
                         K=fold.cv, cv.lim = 10, name.sp = var.y[tu] )
    
    ## OBS = fct(PRED) - CV set -- follow Pineiro 2008
    ##------------------------------------------------
    CV.perf <- asses.lm.(CV.data$predictions,CV.data$obs)
    
    # re-rescale the prediction
    if (RES=="YES"){
      m.temp <- RES.data[,ncol(RES.data)][1]
      sd.temp <- RES.data[,ncol(RES.data)][2]
      
      CV.data$predictions <- (CV.data$predictions*sd.temp) + m.temp
      CV.data$obs <- (CV.data$obs *sd.temp ) + m.temp
    }else{}
    
    # re-transform CV data..
    if (TRANSF=="YES"){
      f <- as.function(inv.transf[[var.list$var.transf[var.list$var.name==var.y[tu]]]])
      CV.data$predictions <- f(x=CV.data$predictions)
      CV.data$obs <- f(x=CV.data$obs)
    }else{}
    
    write.table(CV.data, file=paste(outpath,"/",mod.tech,"_CV_data_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    names(CV.perf) <- paste(rep("CV_",5),names(CV.perf),sep="")
    
    if (i==1) {
      png(paste(outpath,"/",mod.tech,"_OBSvsPRED_rep_",i,".png",sep=""),width = 30, height = 15,units="cm",res=150)
      
      ## Plot observed VS predicted for CV
      plot (CV.data$predictions, CV.data$obs, xlab="Predict_CV", ylab="Observed_CV")
      
      dev.off()
      
    } else {}
    }else{}
    
    ## AIC and AICC of the Best model
    ################################
    if (SVC=="YES"){ 
      RSS <- sum((y.in - (as.numeric(scv.cali[[i]]) *pred.y ))^2) 
    } else{
      RSS <- sum((y.in - pred.y)^2) 
    }
    
    aic.temp <- 2*sum(net.tr$wts!=0) - length(y.in)*log(RSS/length(y.in)) # AIC
    
    aicc.temp <- aic.temp + (2*sum(net.tr$wts!=0)+(sum(net.tr$wts!=0)+1))/
      (length(y.in) -  sum(net.tr$wts!=0)-1) #AICc
    
    if (CR.VAL=="YES"){
      modl.perf.nnet <- rbind(modl.perf.nnet, c(nb.perf=nrow(x.in), M.perf,
                                              nb.vali=nrow(x.vali),V.perf,CV.perf, 
                                              RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
       }else{
         modl.perf.nnet <- rbind(modl.perf.nnet, c(nb.perf=nrow(x.in), M.perf,
                                                   nb.vali=nrow(x.vali),V.perf, 
                                                   RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
      
    }
  } # end of loop repetition
  
  # write all model perf run
  write.table(modl.perf.nnet, file=paste(outpath,"/",mod.tech,"_model_perf_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  
  #write table for all run
  write.table(varimp.nnet, file=paste(outpath,"/",mod.tech,"_var_imp_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  
  # calculate avearge values for model perf and variable importance
  modl.perf.av <- rbind( modl.perf.av, apply(modl.perf.nnet,2,mean) )
  varimp.av <- rbind( varimp.av, apply(varimp.nnet,2,mean) ) 
  
  # name technique update
  mod.name.perf  <- rbind(mod.name.perf, c(species= var.y[tu],var.sel= var.sele[zu],technic=mod.tech) )
  mod.name <- rbind(mod.name, c(species= var.y[tu],var.sel= var.sele[zu],technic=mod.tech) )
  
} # end loop modeling technique
  
##########################################################################################################################################################
if (RSNNS == "NO") {}else{
  
  mod.tech <- "RSNNS"
  outpath <- paste(output,"/",mod.tech,sep="")
  creat.subDir (output,mod.tech)
  
  for (i in 1:n.rep)    
  { 
    # creat x and y table
    ######################
    if (SVC=="YES"){
      y.in <- data.frame(d.cali[[i]][,names(d.cali[[i]])==var.y[tu]])
      y.in <- data.frame(y.in[scv.cali[[i]]==1,1])
      names(y.in) <- var.y[tu]
      
      x.in <- d.cali[[i]][, na.omit( match(names(d.cali[[i]]), 
                                           var.list$var.name[var.list$var.type=="x"] ) ) ]
      x.in <- x.in[scv.cali[[i]]==1,]
      data.in <- cbind(y.in,x.in)
      
      x.vali <- d.vali[[i]][, na.omit( match(names(d.vali[[i]]), 
                                             var.list$var.name[var.list$var.type=="x"] ) ) ]
      x.vali <- x.vali[scv.vali[[i]]==1,]
      y.vali <- data.frame(d.vali[[i]][,names(d.vali[[i]])==var.y[tu]])
      y.vali <- data.frame(y.vali[scv.vali[[i]]==1,1])
      
      p.cal.svc <- as.numeric(row.names(d.cali[[i]]))[scv.cali[[i]]==1]
      p.val.svc <- as.numeric(row.names(d.vali[[i]]))[scv.vali[[i]]==1]
      
    }else{
    
    y.in <- data.frame(d.cali[[i]][,names(d.cali[[i]])==var.y[tu]])
    names(y.in) <- var.y[tu]
    x.in <- d.cali[[i]][, na.omit( match(names(d.cali[[i]]), 
                                         var.list$var.name[var.list$var.type=="x"] ) ) ]
    data.in <- cbind(y.in,x.in)
    
    x.vali <- d.vali[[i]][, na.omit( match(names(d.vali[[i]]), 
                                           var.list$var.name[var.list$var.type=="x"] ) ) ]
    y.vali <- data.frame(d.vali[[i]][,names(d.vali[[i]])==var.y[tu]])
    }
    
    ## Neural Networks (RSNNS) -- RSNNS package
    ##############################################
    ## OPTIMIZED PARAMETER PERFORMANCE OF RSNNS
    ##########################################
    
    if (i ==1) {  cat("> RSNNS Parameter two layer ...", "\n",append = FALSE)
      
      ##(from the FAQ for a commercial neural network software company)  
      ## Number of inputs + outputs) * (2/3) -- Heaton 2008
      nb.hid1 <- round (ncol(data.in) *2/3)
      
      nb.hid2 <- round (nb.hid1*2/3)
      
      nb.hidden =c(nb.hid1,nb.hid2)
      
    } else {}
    
    cat("> RSNNS MODEL run ",i,"sampling  ...", "\n",append = FALSE)
      
      rsnns.tr <- mlp(x.in,y.in, size= nb.hidden, maxit= it.max, learnFunc = "Std_Backpropagation")
      
    # Save Model - R file
    file.name <- paste(outpath,"/",mod.tech,"_model_",i,sep='')
    save(list=c('rsnns.tr'),file=file.name)
    
    ######################
    ## Variable importance
    ######################
    ## weight algorithm (Garson 1991)
    ##-------------------------------
    # compute a simple one layer neural network (var.importance not possible on two layer)
    rsnns.g <- mlp(x.in,y.in, size= nb.hidden[1], maxit= it.max, learnFunc = "Std_Backpropagation")
    
    g  <- garson(rsnns.g)
    
    # Re-order data by row names
    order <- as.numeric(row.names(g$data))
    sort <- sort(order, index.return=TRUE)
    g$data[,1] <- g$data[sort$ix,1] 
    g$data[,2] <- g$data[sort$ix,2] 
    
    varimp.rsnns <- rbind(varimp.rsnns,g$data[,1])
    
    colnames(varimp.rsnns) <- g$data[,2]
    
    ## MODEL PERFORMANCE
    ####################
    # PREDICT THE OBS
    pred.y <- predict(rsnns.tr, x.in) 
    # validate the prediction
    pred.y.vali <- predict(rsnns.tr, x.vali)
    
    ## recombin SVC classification and ML prediction
    ################################################
    if (SVC=="YES"){ 
      pred.y <- rbind(pred.y, matrix( rep(0,nrow(d.cali[[i]])-nrow(pred.y)) )) 
      pred.y.vali <- rbind(pred.y.vali, matrix( rep(0,nrow(d.vali[[i]])-nrow(pred.y.vali)) ) )
      
      y.in <- c(y.in[,1], d.cali[[i]][ scv.cali[[i]]==0,names(d.cali[[i]])==var.y[tu] ] )
      y.vali <- c(y.vali[,1], d.vali[[i]][ scv.vali[[i]]==0,names(d.vali[[i]])==var.y[tu] ] )
      
      p.cali <- as.numeric(row.names(x.in))
      p.vali <- as.numeric(row.names(x.vali))
      
      if (nrow(d.cali[[i]])==length(p.cal.svc)){}else{ 
        p.cali <- c(p.cali,
                    as.numeric(row.names(d.cali[[i]]))[!as.numeric(row.names(d.cali[[i]]))%in%p.cali]) }
      
      if (nrow(d.vali[[i]])==length(p.val.svc)){}else{ 
        p.vali <- c(p.vali,
                    as.numeric(row.names(d.vali[[i]]))[!as.numeric(row.names(d.vali[[i]]))%in%p.vali])}
      
    } else{
      y.in <- y.in[,1]
      y.vali <- y.vali[,1]
      
      p.cali <- as.numeric(row.names(x.in))
      p.vali <- as.numeric(row.names(x.vali))
    }
    
    # PERFORMANCE OF THE MODEL
    #	lin.corr <- lm(y.in ~ pred.y) ;summary(lin.corr)
    M.perf <- asses.lm.(pred.y, y.in)
    V.perf <- asses.lm.(pred.y.vali, y.vali)
    
    names(M.perf) <- paste(rep("M_",5),names(M.perf),sep="")
    names(V.perf) <- paste(rep("V_",5),names(V.perf),sep="")
    
    ##keep data and rescale it
    if (RES=="YES"){
      m.temp <- RES.data[,ncol(RES.data)][1]
      sd.temp <- RES.data[,ncol(RES.data)][2]
      
      pred.y <- (pred.y*sd.temp) + m.temp
      y.in <- (y.in*sd.temp) + m.temp
      pred.y.vali <- (pred.y.vali*sd.temp) + m.temp
      y.vali <- (y.vali*sd.temp) + m.temp
      
    } else{}
    
    ## keep data and retransformed it...
    if (TRANSF=="YES"){
        f <- as.function(inv.transf[[var.list$var.transf[var.list$var.name==var.y[tu]]]])
        cali <- data.frame(cbind(p.cali,rep("CALI",length(pred.y)), f(x=y.in), f(x=pred.y) ))
        vali <- data.frame(cbind(p.vali,rep("VALI",length(pred.y.vali)), f(x=y.vali), f(x=pred.y.vali) ))
      }else{
        cali <- data.frame(cbind(p.cali,rep("CALI",length(pred.y)), y.in, pred.y ))
        vali <- data.frame(cbind(p.vali,rep("VALI",length(pred.y.vali)), y.vali, pred.y.vali ))
      }
    
    names(cali) <-c("id","Type","obs.y","pred.y")
    names(vali) <- c("id","Type","obs.y","pred.y")
    
    cali.vali <- rbind(cali,vali)
    
    write.table(cali.vali, file=paste(outpath,"/",mod.tech,"_cali.vali_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=FALSE,col.names=TRUE, quote=FALSE)
    
    ## CROSS - Validate model
    ##########################
    if (CR.VAL=="YES"){
    CV.data  <- CV.RSNNS.(data.in, nral.hidden= nb.hidden , it.max = it.max, 
                           K=fold.cv, cv.lim = 10, name.sp = var.y[tu])
    
    ## OBS = fct(PRED) - CV set -- follow Pineiro 2008
    ##------------------------
    CV.perf <- asses.lm.(CV.data$predictions,CV.data$obs)
    
    # re-rescale the prediction
    if (RES=="YES"){
      m.temp <- RES.data[,ncol(RES.data)][1]
      sd.temp <- RES.data[,ncol(RES.data)][2]
      
      CV.data$predictions <- (CV.data$predictions*sd.temp) + m.temp
      CV.data$obs <- (CV.data$obs *sd.temp ) + m.temp
    }else{}
    
    # re-transform CV data..
    if (TRANSF=="YES"){
      CV.data$predictions <- f(x=CV.data$predictions)
      CV.data$obs <- f(x=CV.data$obs)
    }else{}
    
    write.table(CV.data, file=paste(outpath,"/",mod.tech,"_CV_data_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    names(CV.perf) <- paste(rep("CV_",5),names(CV.perf),sep="")
    
    if (i ==1) {
      png(paste(outpath,"/",mod.tech,"_OBSvsPRED_rep_",i,".png",sep=""),width = 30, height = 15,units="cm",res=150)
      
        ## Plot observed VS predicted for CV
        plot (CV.data$predictions, CV.data$obs, xlab="Predict_CV", ylab="Observed_CV")
      
      dev.off()
    }else{}
    
    }else{}
    
    ## AIC and AICC of the Best model
    ################################
    if (SVC=="YES"){ 
          RSS <- sum((y.in - (as.numeric(scv.cali[[i]]) *pred.y ) )^2) 
        } else{
          RSS <- sum((y.in - pred.y)^2) 
        }
    
    aic.temp <- 2*sum( weightMatrix(rsnns.tr)!=0 ) - length(y.in)*log(RSS/length(y.in)) # AIC
    
    aicc.temp <- aic.temp + (2*sum( weightMatrix(rsnns.tr)!=0 )+(sum( weightMatrix(rsnns.tr)!=0 )+1))/
      (length(y.in) -  sum( weightMatrix(rsnns.tr)!=0 )-1) #AICc
    
    if (CR.VAL=="YES"){
      modl.perf.rsnns <- rbind(modl.perf.rsnns, 
                             c(nb.perf=nrow(x.in), M.perf,
                               nb.vali=nrow(x.vali),V.perf,CV.perf, 
                               RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
    }else{
      modl.perf.rsnns <- rbind(modl.perf.rsnns, 
                               c(nb.perf=nrow(x.in), M.perf,
                                 nb.vali=nrow(x.vali),V.perf,
                                 RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
      
    }
    
  } # end of loop repetition
  
  # write all model perf run
  write.table(modl.perf.rsnns, file=paste(outpath,"/",mod.tech,"_model_perf_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  
  #write table for all run
  write.table(varimp.rsnns, file=paste(outpath,"/",mod.tech,"_var_imp_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  
  # calculate average values for model perf and variable importance
  modl.perf.av <- rbind( modl.perf.av, apply(modl.perf.rsnns,2,mean) )
  varimp.av <- rbind( varimp.av, apply(varimp.rsnns,2,mean) )
  
  # name technique update
  mod.name.perf  <- rbind(mod.name.perf, c(species= var.y[tu],var.sel= var.sele[zu],technic=mod.tech) )
  mod.name <- rbind(mod.name, c(species= var.y[tu],var.sel= var.sele[zu],technic=mod.tech) )
  
} # end loop modeling technique 
  
##########################################################################################################################################################  
if (ELM == "NO") {}else{
  
    mod.tech <- "ELM"
    outpath <- paste(output,"/",mod.tech,sep="")
    creat.subDir (output,mod.tech)
    
    for (i in 1:n.rep)    
    {
      # creat x and y table
      #####################
      if (SVC=="YES")
      {
        y.in <- data.frame(d.cali[[i]][,names(d.cali[[i]])==var.y[tu]])
        y.in <- data.frame(y.in[scv.cali[[i]]==1,1])
        names(y.in) <- var.y[tu]
        
        x.in <- d.cali[[i]][, na.omit( match(names(d.cali[[i]]), 
                                             var.list$var.name[var.list$var.type=="x"] ) ) ]
        x.in <- x.in[scv.cali[[i]]==1,]
        data.in <- cbind(y.in,x.in)
        
        x.vali <- d.vali[[i]][, na.omit( match(names(d.vali[[i]]), 
                                               var.list$var.name[var.list$var.type=="x"] ) ) ]
        x.vali <- x.vali[scv.vali[[i]]==1,]
        y.vali <- data.frame(d.vali[[i]][,names(d.vali[[i]])==var.y[tu]])
        y.vali <- data.frame(y.vali[scv.vali[[i]]==1,1])
        
        data.vali <- cbind(y.vali,x.vali)
        
        p.cal.svc <- as.numeric(row.names(d.cali[[i]]))[scv.cali[[i]]==1]
        p.val.svc <- as.numeric(row.names(d.vali[[i]]))[scv.vali[[i]]==1]
      
      }else{
        y.in <- data.frame(d.cali[[i]][,names(d.cali[[i]])==var.y[tu]])
        names(y.in) <- var.y[tu]
        
        x.in <- d.cali[[i]][, na.omit( match(names(d.cali[[i]]), 
                                             var.list$var.name[var.list$var.type=="x"] ) ) ]
        data.in <- cbind(y.in,x.in)
        
        x.vali <- d.vali[[i]][, na.omit( match(names(d.vali[[i]]), 
                                               var.list$var.name[var.list$var.type=="x"] ) ) ]
        
        y.vali <- data.frame(d.vali[[i]][,names(d.vali[[i]])==var.y[tu]])
        
        data.vali <- cbind(y.vali,x.vali) 
      }
    ######################
    ### NETWORK ANALYSIS
    #####################                       
    ## OPTIMIZED PARAMETER PERFORMANCE OF ELM
    ##########################################
    
    if (i ==1) {  cat("> ELM Parameter one layer ...", "\n",append = FALSE)
      
      ##(from the FAQ for a commercial neural network software company)  
      ## Number of inputs + outputs) * (2/3) -- Heaton 2008
      nb.hid1 <- round (ncol(data.in) *2/3)
      
      nb.hid2 <- round (( nb.hid1+ncol(y.in))*2/3 )
      
      nb.hidden =c(nb.hid1,nb.hid2)
      
      WT.opt <- 0.01
      
    } else {}
    
    cat("> ELM MODEL run ",i,"sampling  ...", "\n",append = FALSE)
    
    if (i ==1) {  cat("> ELM choose algorithm ...", "\n",append = FALSE)
    
              algo <- c("sig","sin","radbas","hardlim","hardlims","satlins",
                      "tansig","tribas","relu","purelin")
    
              for (s in 1:length(algo)) # test the best algorythm
                  {
                    # train the model
                    elm.tr <- elm_train(x=as.matrix(x.in), y=as.matrix(y.in),nhid=nb.hid1, actfun=algo[s] )
                    
                    # PREDICT THE OBS
                    pred.y <- elm_predict(elm.tr, newdata=as.matrix(data.in)[,2:ncol(data.in)])
                    
                    # PERFORMANCE OF THE MODEL
                    lin.corr <- asses.lm.(pred.y,y.in[,1]) 
                    
                    if (s==1){rse <-lin.corr[6] } else { rse <- c(rse,lin.corr[6]) }
                    
                  }
            # min rse = best algorithm
            b.algo <-  which(rse == min(rse), arr.ind = TRUE) 
            b.algo <-  b.algo[1]
              
            cat("> ", algo[b.algo], " is selected ...", "\n",append = FALSE)
            } else{}
   
    # train elm
    elm.tr <- elm_train(x=as.matrix(x.in), y=as.matrix(y.in),
                        nhid=nb.hid1, actfun=algo[b.algo] )
    
    # Save Model - R file
    file.name<-paste(outpath,"/",mod.tech,"_model_",i,sep='')
    save(list=c('elm.tr'),file=file.name)
    
    ######################
    ## Variable importance
    ######################
    #varimp <- var.imp.rf(elm.tr,x.in,var.list,nperm=n.perm)
    varimp <- var.imp.elm(elm.tr,x.in,names(x.in),nperm=100)
    varimp.elm <- rbind(varimp.elm,(varimp/sum(varimp)) )
    
    ## MODEL PERFORMANCE
    ####################
    ####################
    # PREDICT THE OBS
    pred.y <- elm_predict( elm.tr, newdata=as.matrix(data.in)[,2:ncol(data.in)] ) 
    pred.y.vali <- elm_predict( elm.tr, newdata=as.matrix(data.vali)[,2:ncol(data.vali)] ) 
    
    ## recombin SVC classification and ML prediction
    ################################################
    if (SVC=="YES"){ 
      pred.y <- rbind(pred.y, matrix( rep(0,nrow(d.cali[[i]])-nrow(pred.y)) )) 
      pred.y.vali <- rbind(pred.y.vali, matrix( rep(0,nrow(d.vali[[i]])-nrow(pred.y.vali)) ) ) 
      
      y.in <- c(y.in[,1], d.cali[[i]][ scv.cali[[i]]==0,names(d.cali[[i]])==var.y[tu] ] )
      y.vali <- c(y.vali[,1], d.vali[[i]][ scv.vali[[i]]==0,names(d.vali[[i]])==var.y[tu] ] )
      
      p.cali <- as.numeric(row.names(x.in))
      p.vali <- as.numeric(row.names(x.vali))
      
      if (nrow(d.cali[[i]])==length(p.cal.svc)){}else{ 
        p.cali <- c(p.cali,
                    as.numeric(row.names(d.cali[[i]]))[!as.numeric(row.names(d.cali[[i]]))%in%p.cali]) }
      
      if (nrow(d.vali[[i]])==length(p.val.svc)){}else{ 
        p.vali <- c(p.vali,
                    as.numeric(row.names(d.vali[[i]]))[!as.numeric(row.names(d.vali[[i]]))%in%p.vali])}
      
    } else{
      y.in <- y.in[,1]
      y.vali <- y.vali[,1]
      
      p.cali <- as.numeric(row.names(x.in))
      p.vali <- as.numeric(row.names(x.vali))
    }
    
    # PERFORMANCE OF THE MODEL
    #	lin.corr <- lm(y.in ~ pred.y) ;summary(lin.corr)
    M.perf <- asses.lm.(pred.y, y.in )
    V.perf <- asses.lm.(pred.y.vali, y.vali )
    
    names(M.perf) <- paste(rep("M_",5),names(M.perf),sep="")
    names(V.perf) <- paste(rep("V_",5),names(V.perf),sep="")
    
    ##keep data and rescale it
    if (RES=="YES"){
      m.temp <- RES.data[,ncol(RES.data)][1]
      sd.temp <- RES.data[,ncol(RES.data)][2]
      
      pred.y <- (pred.y*sd.temp) + m.temp
      y.in <- (y.in*sd.temp) + m.temp
      pred.y.vali <- (pred.y.vali*sd.temp) + m.temp
      y.vali <- (y.vali*sd.temp) + m.temp
      
    } else{}
    
    ## keep data and retransformed it...
    if (TRANSF=="YES"){
        f <- as.function(inv.transf[[var.list$var.transf[var.list$var.name==var.y[tu]]]])
        cali <- data.frame(cbind(p.cali,rep("CALI",length(pred.y)), f(x=y.in), f(x=pred.y) ))
        vali <- data.frame(cbind(p.vali,rep("VALI",length(pred.y.vali)), f(x=y.vali), f(x=pred.y.vali) ))
      }else{
        cali <- data.frame(cbind(p.cali,rep("CALI",length(pred.y)), y.in, pred.y ))
        vali <- data.frame(cbind(p.vali,rep("VALI",length(pred.y.vali)), y.vali, pred.y.vali ))
      }
    
    names(cali) <-c("id","Type","obs.y","pred.y")
    names(vali) <- c("id","Type","obs.y","pred.y")  	
    
    cali.vali <- rbind(cali,vali)
    
    write.table(cali.vali, file=paste(outpath,"/",mod.tech,"_cali.vali_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    ## CROSS - Validate model
    #########################
    if (CR.VAL=="YES"){
    CV.data <- CV.ELM. (data.in, elm.tr, K=fold.cv, cv.lim = 10, name.sp = var.y[tu] )
    
    ## OBS = fct(PRED) - CV set -- follow Pineiro 2008
    ##------------------------
    CV.perf <- asses.lm.(CV.data$predictions,CV.data$obs)
    
    # re-rescale the prediction
    if (RES=="YES"){
      m.temp <- RES.data[,ncol(RES.data)][1]
      sd.temp <- RES.data[,ncol(RES.data)][2]
      
      CV.data$predictions <- (CV.data$predictions*sd.temp) + m.temp
      CV.data$obs <- (CV.data$obs *sd.temp ) + m.temp
    }else{}
    
    # re-transform CV data..
    if (TRANSF=="YES"){
      CV.data$predictions <- f(x=CV.data$predictions)
      CV.data$obs <- f(x=CV.data$obs)
    }else{}
    
    write.table(CV.data, file=paste(outpath,"/",mod.tech,"_CV_data_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    names(CV.perf) <- paste(rep("CV_",5),names(CV.perf),sep="")
    
    if (i==1) {
      png(paste(outpath,"/",mod.tech,"_OBSvsPRED_rep_",i,".png",sep=""),width = 30, 
          height = 15,units="cm",res=150)
      
        ## Plot observed VS predicted for CV
        plot (CV.data$predictions, CV.data$obs, xlab="Predict_CV", ylab="Observed_CV")
      
      dev.off()
      
    } else {}
    }else{}
    
    ## AIC and AICC of the Best model
    ################################
    if (SVC=="YES"){ 
          RSS <- sum((y.in - (as.numeric(scv.cali[[i]]) *pred.y ) )^2) 
        } else{
          RSS <- sum( (y.in - pred.y)^2) 
        }
    
    aic.temp <- 2*sum(elm.tr$outweight!=0) - length(y.in)*log(RSS/length(y.in)) # AIC
    
    aicc.temp <- aic.temp + (2*sum(elm.tr$outweight!=0)+(sum(elm.tr$outweight!=0)+1))/
      (length(y.in) -  sum(elm.tr$outweight!=0)-1) #AICc
    
    if (CR.VAL=="YES"){
      modl.perf.elm <- rbind(modl.perf.elm, 
                           c(nb.perf=nrow(x.in), M.perf,
                             nb.vali=nrow(x.vali),V.perf,CV.perf, 
                             RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
    }else{
      modl.perf.elm <- rbind(modl.perf.elm, 
                             c(nb.perf=nrow(x.in), M.perf,
                               nb.vali=nrow(x.vali),V.perf, 
                               RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
      
    }
    
  } # end of loop repetition

  # write all model perf run
  write.table(modl.perf.elm, file=paste(outpath,"/",mod.tech,"_model_perf_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)

  #write table for all run
  write.table(varimp.elm, file=paste(outpath,"/",mod.tech,"_var_imp_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
  # calculate avearge values for model perf and variable importance
  modl.perf.av <- rbind( modl.perf.av, apply(modl.perf.elm,2,mean) )
  varimp.av <- rbind( varimp.av, apply(varimp.elm,2,mean) )
  
  # name technique update
  mod.name.perf  <- rbind(mod.name.perf, c(species= var.y[tu],var.sel= var.sele[zu],technic=mod.tech) )
  mod.name <- rbind(mod.name, c(species= var.y[tu],var.sel= var.sele[zu],technic=mod.tech) )
  
} # end loop modeling technique
  
#####################################################################################################################################
if (RF == "NO") {}else{
  
  mod.tech <- "RF"
  outpath <- paste(output,"/",mod.tech,sep="")
  creat.subDir (output,mod.tech) 
  
  for (i in 1:n.rep)    
  {
    # creat x and y table
    #####################
    if (SVC=="YES")
    { 
      y.in <- data.frame(d.cali[[i]][,names(d.cali[[i]])==var.y[tu]])
      y.in <- data.frame(y.in[scv.cali[[i]]==1,1])
      names(y.in) <- var.y[tu]
      
      x.in <- d.cali[[i]][, na.omit( match(names(d.cali[[i]]), 
                                           var.list$var.name[var.list$var.type=="x"] ) ) ]
      x.in <- x.in[scv.cali[[i]]==1,]
      data.in <- cbind(y.in,x.in)
      
      x.vali <- d.vali[[i]][, na.omit( match(names(d.vali[[i]]), 
                                             var.list$var.name[var.list$var.type=="x"] ) ) ]
      x.vali <- x.vali[scv.vali[[i]]==1,]
      
      y.vali <- data.frame(d.vali[[i]][,names(d.vali[[i]])==var.y[tu]])
      y.vali <- data.frame(y.vali[scv.vali[[i]]==1,1])
      
      p.cal.svc <- as.numeric(row.names(d.cali[[i]]))[scv.cali[[i]]==1]
      p.val.svc <- as.numeric(row.names(d.vali[[i]]))[scv.vali[[i]]==1]
      
      }else{
      y.in <- data.frame(d.cali[[i]][,names(d.cali[[i]])==var.y[tu]])
      names(y.in) <- var.y[tu]
      x.in <- d.cali[[i]][, na.omit( match(names(d.cali[[i]]), 
                                           var.list$var.name[var.list$var.type=="x"] ) ) ]
      data.in <- cbind(y.in,x.in)
      
      x.vali <- d.vali[[i]][, na.omit( match(names(d.vali[[i]]), 
                                             var.list$var.name[var.list$var.type=="x"] ) ) ]
      y.vali <- data.frame(d.vali[[i]][,names(d.vali[[i]])==var.y[tu]])
    }
    #####################
    ## Random Forest
    #####################
    cat("> RANDOM FOREST MODEL run ",i,"sampling  ...", "\n",append = FALSE)
    
    # First we tune the parameter
    if (i==1) { rfpa <- tryCatch({tuneRF (x = x.in,y = t(y.in),data = data.in, 
                                          ntree = nb.tree, na.action=na.omit)})
          }else{}
    
    if (overfit=="YES") {nb.hid1 <- round (ncol(data.in) *2/3)
                  rf <- tryCatch({ randomForest(eval(parse(text= paste(colnames(y.in),"~.",sep=""))), 
                                  data =data.in,ntree = nb.tree, mtry = rfpa,importance=TRUE, 
                                  na.action=na.omit, maxnodes=nb.hid1) })
        }else{ 
                  rf <- tryCatch({ randomForest(eval(parse(text= paste(colnames(y.in),"~.",sep=""))), data =data.in, 
                      ntree = nb.tree, mtry = rfpa,importance=TRUE, na.action=na.omit ) })      
        }
    # Save Model - R file
    file.name <- paste(outpath,"/",mod.tech,"_model_",i,sep='')
    save(list=c('rf'),file=file.name) 
    
    ######################
    ## Variable importance type 2 mean decrease in node impurity
    ######################
    varimp <- rf$importance[,2]
    varimp.rf <- rbind(varimp.rf,(varimp/sum(varimp)) ) 
    
    ## MODEL PERFORMANCE
    ####################
    # PREDICT THE OBS
    pred.y <- predict( rf, newdata=x.in ) 
    pred.y.vali <- predict(rf, newdata=x.vali )
    
    ## recombin SVC classification and ML prediction
    ################################################
    if (SVC=="YES"){  
      pred.y <- c(pred.y,  rep(0,nrow(d.cali[[i]])-length(pred.y) ) )
      pred.y.vali <- c(pred.y.vali,  rep(0,nrow(d.vali[[i]])-length(pred.y.vali))  ) 
      
      y.in <- c(y.in[,1], d.cali[[i]][ scv.cali[[i]]==0,names(d.cali[[i]])==var.y[tu] ] )
      y.vali <- c(y.vali[,1], d.vali[[i]][ scv.vali[[i]]==0,names(d.vali[[i]])==var.y[tu] ] )
      
      p.cali <- as.numeric(row.names(x.in))
      p.vali <- as.numeric(row.names(x.vali))
      
      if (nrow(d.cali[[i]])==length(p.cal.svc)){}else{ 
        p.cali <- c(p.cali,
                    as.numeric(row.names(d.cali[[i]]))[!as.numeric(row.names(d.cali[[i]]))%in%p.cali]) }
      
      if (nrow(d.vali[[i]])==length(p.val.svc)){}else{ 
        p.vali <- c(p.vali,
                    as.numeric(row.names(d.vali[[i]]))[!as.numeric(row.names(d.vali[[i]]))%in%p.vali])}
      
    } else{
      y.in <- y.in[,1]
      y.vali <- y.vali[,1]
      
      p.cali <- as.numeric(row.names(x.in))
      p.vali <- as.numeric(row.names(x.vali))
    }
  
    # PERFORMANCE OF THE MODEL
    #	lin.corr <- lm(y.in ~ pred.y) ;summary(lin.corr)
    M.perf <- asses.lm.(pred.y, y.in )
    V.perf <- asses.lm.(pred.y.vali, y.vali )
    
    names(M.perf) <- paste(rep("M_",5),names(M.perf),sep="")
    names(V.perf) <- paste(rep("V_",5),names(V.perf),sep="")
    
    ##keep data and rescale it
    if (RES=="YES"){
      m.temp <- RES.data[,ncol(RES.data)][1]
      sd.temp <- RES.data[,ncol(RES.data)][2]
      
      pred.y <- (pred.y*sd.temp) + m.temp
      y.in <- (y.in*sd.temp) + m.temp
      pred.y.vali <- (pred.y.vali*sd.temp) + m.temp
      y.vali <- (y.vali*sd.temp) + m.temp
      
    } else{}
    
    ## keep data and retransformed it...
    if (TRANSF=="YES"){
        f <- as.function(inv.transf[[var.list$var.transf[var.list$var.name==var.y[tu]]]])
        cali <- data.frame(cbind(p.cali,rep("CALI",length(pred.y)), f(x=y.in), f(x=pred.y) ))
        vali <- data.frame(cbind(p.vali, rep("VALI",length(pred.y.vali)), f(x=y.vali), f(x=pred.y.vali) ))
      }else{
        cali <- data.frame(cbind(p.cali,rep("CALI",length(pred.y)), y.in, pred.y ))
        vali <- data.frame(cbind(p.vali, rep("VALI",length(pred.y.vali)), y.vali, pred.y.vali ))
      }
    
    names(cali) <-c("id","Type","obs.y","pred.y")
    names(vali) <- c("id","Type","obs.y","pred.y") 
    
    cali.vali <- rbind(cali,vali)
    
    write.table(cali.vali, file=paste(outpath,"/",mod.tech,"_cali.vali_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    ###CROSS - Validate model
    ##########################
    if (CR.VAL=="YES"){
    CV.data  <- tryCatch({ CV.RF.(rf, data.in,  K=fold.cv, cv.lim = 10, name.sp = var.y[tu])  })
    
    ## OBS = fct(PRED) - CV set -- follow Pineiro 2008
    ##------------------------
    CV.perf <- asses.lm.(CV.data$predictions,CV.data$obs)
    
    # re-rescale the prediction
    if (RES=="YES"){
      m.temp <- RES.data[,ncol(RES.data)][1]
      sd.temp <- RES.data[,ncol(RES.data)][2]
      
      CV.data$predictions <- (CV.data$predictions*sd.temp) + m.temp
      CV.data$obs <- (CV.data$obs *sd.temp ) + m.temp
    }else{}
    
    # re-transform CV data..
    if (TRANSF=="YES"){
      CV.data$predictions <- f(x=CV.data$predictions)
      CV.data$obs <- f(x=CV.data$obs)
    }else{}
    
    write.table(CV.data, file=paste(outpath,"/",mod.tech,"_CV_data_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    names(CV.perf) <- paste(rep("CV_",5),names(CV.perf),sep="")
    
    if (i ==1) {
      png(paste(outpath,"/",mod.tech,"_OBSvsPRED_rep_",i,".png",sep=""),width = 30, 
          height = 15,units="cm",res=150)
      
        ## Plot observed VS predicted for CV
        plot (CV.data$predictions, CV.data$obs, xlab="Predict_CV", ylab="Observed_CV")
      
      dev.off()
    }else{} 
    }else{}
    
    ## AIC and AICC of the Best model
    ################################
    if (SVC=="YES"){ 
          RSS <- sum((y.in - (as.numeric(scv.cali[[i]]) *pred.y ))^2) 
        } else{
          RSS <- sum((y.in - pred.y)^2) 
        }
     
    aic.temp <- 2*mean(treesize(rf)) - length(y.in)*log(RSS/length(y.in)) # AIC
    
    aicc.temp <- aic.temp + (2*mean(treesize(rf))+(mean(treesize(rf))+1))/
      (length(y.in) -  mean(treesize(rf))-1) #AICc
    
    if (CR.VAL=="YES"){
    modl.perf.rf <- rbind(modl.perf.rf,
                          c(nb.perf=nrow(x.in), M.perf,
                            nb.vali=nrow(x.vali),V.perf,CV.perf, 
                            RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
    }else{
      modl.perf.rf <- rbind(modl.perf.rf,
                            c(nb.perf=nrow(x.in), M.perf,
                              nb.vali=nrow(x.vali),V.perf, 
                              RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
      
    }
    
  } # end of loop repetition
  
  # write all model perf run
  write.table(modl.perf.rf, file=paste(outpath,"/",mod.tech,"_model_perf_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  
  #write table for all run
  write.table(varimp.rf, file=paste(outpath,"/",mod.tech,"_var_imp_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
      
  # calculate avearge values for model perf and variable importance
  modl.perf.av <- rbind( modl.perf.av, apply(modl.perf.rf,2,mean) )
  varimp.av <- rbind( varimp.av, apply(varimp.rf,2,mean) )
  
  # name technique update
  mod.name.perf  <- rbind(mod.name.perf, c(species= var.y[tu],var.sel= var.sele[zu],technic=mod.tech) )
  mod.name <- rbind(mod.name, c(species= var.y[tu],var.sel= var.sele[zu],technic=mod.tech) )
  
} # end loop modeling technique

#####################################################################################################################################
if (CSRF == "NO") {}else{
  mod.tech <- "CSRF"
  outpath <- paste(output,"/",mod.tech,sep="")
  creat.subDir (output,mod.tech)  
  
  for (i in 1:n.rep)    
  {
    # creat x and y table
    #####################
    if (SVC=="YES")
    { 
      y.in <- data.frame(d.cali[[i]][,names(d.cali[[i]])==var.y[tu]])
      y.in <- data.frame(y.in[scv.cali[[i]]==1,1])
      names(y.in) <- var.y[tu]
      
      x.in <- d.cali[[i]][, na.omit( match(names(d.cali[[i]]), 
                                           var.list$var.name[var.list$var.type=="x"] ) ) ]
      x.in <- x.in[scv.cali[[i]]==1,]
      data.in <- cbind(x.in,y.in)
      
      x.vali <- d.vali[[i]][, na.omit( match(names(d.vali[[i]]), 
                                             var.list$var.name[var.list$var.type=="x"] ) ) ]
      x.vali <- x.vali[scv.vali[[i]]==1,]
      y.vali <- data.frame(d.vali[[i]][,names(d.vali[[i]])==var.y[tu]])
      y.vali <- data.frame(y.vali[scv.vali[[i]]==1,1])
      
      p.cal.svc <- as.numeric(row.names(d.cali[[i]]))[scv.cali[[i]]==1]
      p.val.svc <- as.numeric(row.names(d.vali[[i]]))[scv.vali[[i]]==1]
      
      }else{
      y.in <- data.frame(d.cali[[i]][,names(d.cali[[i]])==var.y[tu]])
      names(y.in) <- var.y[tu]
      x.in <- d.cali[[i]][, na.omit( match(names(d.cali[[i]]), 
                                         var.list$var.name[var.list$var.type=="x"] ) ) ]
      data.in <- cbind(x.in,y.in)
    
      x.vali <- d.vali[[i]][, na.omit( match(names(d.vali[[i]]), 
                                           var.list$var.name[var.list$var.type=="x"] ) ) ]
      y.vali <- data.frame(d.vali[[i]][,names(d.vali[[i]])==var.y[tu]])
    }
    ##############################
    ## Case specific Random forest
    ##############################
    cat("> Case specific RANDOM FOREST MODEL run ",i,"sampling  ...", "\n",append = FALSE)
    
    # First we tune the parameter
    if (i==1) {tryCatch( rfpa <- tuneRF (x = x.in,y = t(y.in),data = data.in, ntree = nb.tree, na.action=na.omit) )
                rfpa <- rfpa[match( max(rfpa[,2]), rfpa[,2] ),1]
          }else{}
    
    if (overfit=="YES") {
        nb.hid1 <- round (ncol(data.in) *2/3)
    
        csrf <- ranger(eval(parse(text= paste(colnames(y.in),"~.",sep=""))), data =as.data.frame(data.in), 
                   num.trees = nb.tree,mtry= rfpa, importance="impurity", 
                   min.node.size = nb.hid1, splitrule = "maxstat", alpha = alpha.lim)
      }else{
        csrf <- ranger(eval(parse(text= paste(colnames(y.in),"~.",sep=""))), data =as.data.frame(data.in), 
                     num.trees = nb.tree,mtry= rfpa, importance="impurity")
      }
    
    # Save Model - R file
    file.name <- paste(outpath,"/",mod.tech,"_model_",i,sep='')
    save(list=c('csrf'),file=file.name) 
    
    ######################
    ## Variable importance
    ######################
    varimp <- csrf$variable.importance
    varimp.csrf <- rbind(varimp.csrf,(varimp/sum(varimp)) ) 
    
    ## MODEL PERFORMANCE
    ####################
    # PREDICT THE OBS
    pred.y <- predict(csrf, dat=x.in) 
    pred.y.vali <- predict(csrf, dat=x.vali )
    
    ## recombin SVC classification and ML prediction
    ################################################
    if (SVC=="YES"){  
      pred.y <- c(pred.y$predictions,  rep(0,nrow(d.cali[[i]])-length(pred.y$predictions)) ) 
      pred.y.vali <- c(pred.y.vali$predictions, rep(0,nrow(d.vali[[i]])-length(pred.y.vali$predictions)) ) 
      
      y.in <- c(y.in[,1], d.cali[[i]][ scv.cali[[i]]==0,names(d.cali[[i]])==var.y[tu] ] )
      y.vali <- c(y.vali[,1], d.vali[[i]][ scv.vali[[i]]==0,names(d.vali[[i]])==var.y[tu] ] )
      
      p.cali <- as.numeric(row.names(x.in))
      p.vali <- as.numeric(row.names(x.vali))
      
      if (nrow(d.cali[[i]])==length(p.cal.svc)){}else{ 
        p.cali <- c(p.cali,
                    as.numeric(row.names(d.cali[[i]]))[!as.numeric(row.names(d.cali[[i]]))%in%p.cali]) }
      
      if (nrow(d.vali[[i]])==length(p.val.svc)){}else{ 
        p.vali <- c(p.vali,
                    as.numeric(row.names(d.vali[[i]]))[!as.numeric(row.names(d.vali[[i]]))%in%p.vali])}
      
    } else{
      y.in <- y.in[,1]
      y.vali <- y.vali[,1]
      
      p.cali <- as.numeric(row.names(x.in))
      p.vali <- as.numeric(row.names(x.vali))
    }
    
    # PERFORMANCE OF THE MODEL
    #	lin.corr <- lm(y.in ~ pred.y) ;summary(lin.corr)
    M.perf <- asses.lm.(pred.y, y.in )
    V.perf <- asses.lm.(pred.y.vali, y.vali )
    
    names(M.perf) <- paste(rep("M_",5),names(M.perf),sep="")
    names(V.perf) <- paste(rep("V_",5),names(V.perf),sep="")
    
    ##keep data and rescale it
    if (RES=="YES"){
      m.temp <- RES.data[,ncol(RES.data)][1]
      sd.temp <- RES.data[,ncol(RES.data)][2]
      
      pred.y <- (pred.y*sd.temp) + m.temp
      y.in <- (y.in*sd.temp) + m.temp
      pred.y.vali <- (pred.y.vali*sd.temp) + m.temp
      y.vali <- (y.vali*sd.temp) + m.temp
      
    } else{}
    
    ## keep data and retransformed it...
    if (TRANSF=="YES"){
         f <- as.function(inv.transf[[var.list$var.transf[var.list$var.name==var.y[tu]]]])
        cali <- data.frame(cbind(p.cali, rep("CALI",length(pred.y)), f(x=y.in), f(x=pred.y) ))
        vali <- data.frame(cbind(p.vali, rep("VALI",length(pred.y.vali)), f(x=y.vali), f(x=pred.y.vali) ))
      }else{
        cali <- data.frame(cbind(p.cali, rep("CALI",length(pred.y)), y.in, pred.y ))
        vali <- data.frame(cbind(p.vali, rep("VALI",length(pred.y.vali)), y.vali, pred.y.vali ))
      }
    
    names(cali) <-c("id","Type","obs.y","pred.y")
    names(vali) <- c("id","Type","obs.y","pred.y") 
    
    cali.vali <- rbind(cali,vali)
    
    write.table(cali.vali, file=paste(outpath,"/",mod.tech,"_cali.vali_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    ###CROSS - Validate model
    #########################
    if (CR.VAL=="YES"){
    CV.data  <- tryCatch({ CV.CSRF.(csrf, as.data.frame(data.in),  K=fold.cv, 
                                    cv.lim = 10, name.sp = var.y[tu]) }) 
    
    ## OBS = fct(PRED) - CV set -- follow Pineiro 2008
    ##------------------------
    CV.perf <- tryCatch({asses.lm.(CV.data$predictions,CV.data$obs) })
    
    # re-rescale the prediction
    if (RES=="YES"){
      m.temp <- RES.data[,ncol(RES.data)][1]
      sd.temp <- RES.data[,ncol(RES.data)][2]
      
      CV.data$predictions <- (CV.data$predictions*sd.temp) + m.temp
      CV.data$obs <- (CV.data$obs *sd.temp ) + m.temp
    }else{}
    
    # re-transform CV data..
    if (TRANSF=="YES"){
      CV.data$predictions <- f(x=CV.data$predictions)
      CV.data$obs <- f(x=CV.data$obs)
    }else{}
    
    write.table(CV.data, file=paste(outpath,"/",mod.tech,"_CV_data_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    names(CV.perf) <- paste(rep("CV_",5),names(CV.perf),sep="")
    
      if (i ==1) {
        png(paste(outpath,"/",mod.tech,"_OBSvsPRED_rep_",i,".png",sep=""),width = 30, 
            height = 15,units="cm",res=150)
        
        ## Plot observed VS predicted for CV
        plot (CV.data$predictions, CV.data$obs, xlab="Predict_CV", ylab="Observed_CV")
        
        dev.off()
      }else{} 
    }else{}
    
    ## AIC and AICC of the Best model
    ################################
    if (SVC=="YES"){ 
          RSS <- sum((y.in - (as.numeric(scv.cali[[i]]) * pred.y) )^2) 
        } else{
          RSS <- sum((y.in - pred.y )^2)  
        }
    
    aic.temp <- NA # AIC
    aicc.temp <- NA #AICc
    
    if (CR.VAL=="YES"){
      modl.perf.csrf <- rbind(modl.perf.csrf, 
                            c(nb.perf=nrow(x.in), M.perf,
                              nb.vali=nrow(x.vali),V.perf,CV.perf, 
                              RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
    }else{
      modl.perf.csrf <- rbind(modl.perf.csrf, 
                              c(nb.perf=nrow(x.in), M.perf,
                                nb.vali=nrow(x.vali),V.perf, 
                                RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
      
    }
    
  } # end of loop repetition
  
  # write all model perf run
  write.table(modl.perf.csrf, file=paste(outpath,"/",mod.tech,"_model_perf_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  
  #write table for all run
  write.table(varimp.csrf, file=paste(outpath,"/",mod.tech,"_var_imp_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  
  # calculate avearge values for model perf and variable importance
  modl.perf.av <- rbind( modl.perf.av, apply(modl.perf.csrf,2,mean) )
  varimp.av <- rbind( varimp.av, apply(varimp.csrf,2,mean) )
  
  # name technique update
  mod.name.perf  <- rbind(mod.name.perf, c(species= var.y[tu],var.sel= var.sele[zu],technic=mod.tech) )
  mod.name <- rbind(mod.name, c(species= var.y[tu],var.sel= var.sele[zu],technic=mod.tech) )
  
} # end loop modeling technique
  
#####################################################################################################################################
if (BAG == "NO") {}else{
  
  mod.tech <- "BAG"
  outpath <- paste(output,"/",mod.tech,sep="")
  creat.subDir (output,mod.tech)
    
  for (i in 1:n.rep)    
  {
    ## creat x and y table
    #####################
    if (SVC=="YES")
    {
      y.in <- data.frame(d.cali[[i]][,names(d.cali[[i]])==var.y[tu]])
      y.in <- data.frame(y.in[scv.cali[[i]]==1,1])
      names(y.in) <- var.y[tu]
      
      x.in <- d.cali[[i]][, na.omit( match(names(d.cali[[i]]), 
                                           var.list$var.name[var.list$var.type=="x"] ) ) ]
      x.in <- x.in[scv.cali[[i]]==1,]
      data.in <- cbind(x.in,y.in)
      
      x.vali <- d.vali[[i]][, na.omit( match(names(d.vali[[i]]), 
                                             var.list$var.name[var.list$var.type=="x"] ) ) ]
      x.vali <- x.vali[scv.vali[[i]]==1,]
      y.vali <- data.frame(d.vali[[i]][,names(d.vali[[i]])==var.y[tu]])
      y.vali <- data.frame(y.vali[scv.vali[[i]]==1,1])
      
      data.vali <- cbind(y.vali,x.vali)
      
      p.cal.svc <- as.numeric(row.names(d.cali[[i]]))[scv.cali[[i]]==1]
      p.val.svc <- as.numeric(row.names(d.vali[[i]]))[scv.vali[[i]]==1]
    }else{
    y.in <- data.frame(d.cali[[i]][,names(d.cali[[i]])==var.y[tu]])
    names(y.in) <- var.y[tu]
    x.in <- d.cali[[i]][, na.omit( match(names(d.cali[[i]]), 
                                         var.list$var.name[var.list$var.type=="x"] ) ) ]
    data.in <- cbind(x.in,y.in)
    
    x.vali <- d.vali[[i]][, na.omit( match(names(d.vali[[i]]), 
                                           var.list$var.name[var.list$var.type=="x"] ) ) ]
    y.vali <- data.frame(d.vali[[i]][,names(d.vali[[i]])==var.y[tu]])
    
    data.vali <- cbind(y.vali,x.vali)
    }
    ######################
    ### BAGGING ANALYSIS -- adabag
    #####################                       
    ## OPTIMIZED PARAMETER PERFORMANCE OF BAG
    ##########################################
    cat("> BAG MODEL run ",i,"sampling  ...", "\n",append = FALSE)
    
    bag.tr <- ipredbagg (y.in[,1], x.in, nbag=nb.tree)
    
    # Save Model - R file
    file.name<-paste(outpath,"/",mod.tech,"_model_",i,sep='')
    save(list=c('bag.tr'),file=file.name)
    
    ######################
    ## Variable importance --> node impurity (similar as permutation)
    ######################
    for (n in 1: nb.tree)
      { 
        varimp.t <- bag.tr$mtrees[[n]]$btree$variable.importance # extract the var. imp
        
        if (is.null(varimp.t)) { } else{ # check if not null
            sort <- sort(names(varimp.t), index.return=TRUE) # sort by alphabetic order
            varimp.t <- varimp.t[sort$ix]
            varimp.t <- varimp.t[names(varimp.t)!=var.y[tu]]
            # remove the var. dependent
            if (n==1){ varimp <- varimp.t} else { varimp <- rbind( varimp, varimp.t )}
          }              
        }
    
    varimp.bag <- rbind( varimp.bag, apply(varimp, 2, sum)/ sum(varimp) )
    
    ## MODEL PERFORMANCE
    ####################
    # PREDICT THE OBS
    pred.y <- predict(bag.tr, newdata=data.frame(data.in), type="raw",na.rm = TRUE)
    pred.y.vali <- predict(bag.tr, newdata=data.frame(data.vali), type="raw",na.rm = TRUE)
    
    ## recombin SVC classification and ML prediction
    ################################################
    if (SVC=="YES"){  
      pred.y <- c(pred.y,  rep(0,nrow(d.cali[[i]])-length(pred.y)) ) 
      pred.y.vali <- c(pred.y.vali,  rep(0,nrow(d.vali[[i]])-length(pred.y.vali)) ) 
      
      y.in <- c(y.in[,1], d.cali[[i]][ scv.cali[[i]]==0,names(d.cali[[i]])==var.y[tu] ] )
      y.vali <- c(y.vali[,1], d.vali[[i]][ scv.vali[[i]]==0,names(d.vali[[i]])==var.y[tu] ] )
      
      p.cali <- as.numeric(row.names(x.in))
      p.vali <- as.numeric(row.names(x.vali))
      
      if (nrow(d.cali[[i]])==length(p.cal.svc)){}else{ 
        p.cali <- c(p.cali,
                    as.numeric(row.names(d.cali[[i]]))[!as.numeric(row.names(d.cali[[i]]))%in%p.cali]) }
      
      if (nrow(d.vali[[i]])==length(p.val.svc)){}else{ 
        p.vali <- c(p.vali,
                    as.numeric(row.names(d.vali[[i]]))[!as.numeric(row.names(d.vali[[i]]))%in%p.vali])}
      
      
    } else{
      y.in <- y.in[,1]
      y.vali <- y.vali[,1]
      
      p.cali <- as.numeric(row.names(x.in))
      p.vali <- as.numeric(row.names(x.vali))
    }
    
    # PERFORMANCE OF THE MODEL
    #	lin.corr <- lm(y.in ~ pred.y) ;summary(lin.corr)
    M.perf <- asses.lm.(pred.y, y.in )
    V.perf <- asses.lm.(pred.y.vali, y.vali)
    
    names(M.perf) <- paste(rep("M_",5),names(M.perf),sep="")
    names(V.perf) <- paste(rep("V_",5),names(V.perf),sep="")
    
    ##keep data and rescale it
    if (RES=="YES"){
      m.temp <- RES.data[,ncol(RES.data)][1]
      sd.temp <- RES.data[,ncol(RES.data)][2]
      
      pred.y <- (pred.y*sd.temp) + m.temp
      y.in <- (y.in*sd.temp) + m.temp
      pred.y.vali <- (pred.y.vali*sd.temp) + m.temp
      y.vali <- (y.vali*sd.temp) + m.temp
      
    } else{}
    
    ## keep data and retransformed it...
    if (TRANSF=="YES"){
        f <- as.function(inv.transf[[var.list$var.transf[var.list$var.name==var.y[tu]]]])
        cali <- data.frame(cbind(p.cali,rep("CALI",length(pred.y)), f(x=y.in), f(x=pred.y) ))
        vali <- data.frame(cbind(p.vali,rep("VALI",length(pred.y.vali)), f(x=y.vali), f(x=pred.y.vali) ))
      }else{
        cali <- data.frame(cbind(p.cali, rep("CALI",length(pred.y)), y.in, pred.y ))
        vali <- data.frame(cbind(p.vali, rep("VALI",length(pred.y.vali)), y.vali, pred.y.vali ))
      }
    
    names(cali) <-c("id","Type","obs.y","pred.y")
    names(vali) <- c("id","Type","obs.y","pred.y")  
    
    cali.vali <- rbind(cali,vali)
    
    write.table(cali.vali, file=paste(outpath,"/",mod.tech,"_cali.vali_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    ## CROSS - Validate model
    #########################
    if (CR.VAL=="YES"){
    CV.data  <- CV.BAG. (data.in, K= fold.cv, cv.lim = 10, nb.tree = nb.tree, name.sp = var.y[tu])
    
    ## OBS = fct(PRED) - CV set -- follow Pineiro 2008
    ##------------------------
    CV.perf <- asses.lm.(CV.data$predictions,CV.data$obs)
    
    # re-rescale the prediction
    if (RES=="YES"){
      m.temp <- RES.data[,ncol(RES.data)][1]
      sd.temp <- RES.data[,ncol(RES.data)][2]
      
      CV.data$predictions <- (CV.data$predictions*sd.temp) + m.temp
      CV.data$obs <- (CV.data$obs *sd.temp ) + m.temp
    }else{}
    
    # re-transform CV data..
    if (TRANSF=="YES"){
      CV.data$predictions <- f(x=CV.data$predictions)
      CV.data$obs <- f(x=CV.data$obs)
    }else{}
    
    write.table(CV.data, file=paste(outpath,"/",mod.tech,"_CV_data_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    names(CV.perf) <- paste(rep("CV_",5),names(CV.perf),sep="")
    
    if (i==1) {
      png(paste(outpath,"/",mod.tech,"_OBSvsPRED_rep_",i,".png",sep=""),width = 30, height = 15,units="cm",res=150)
      
        ## Plot observed VS predicted for CV
        plot (CV.data$predictions, CV.data$obs, xlab="Predict_CV", ylab="Observed_CV")
      
      dev.off()
      
    } else {}
    
    }else{}
    
    ## AIC and AICC of the Best model
    ################################
    if (SVC=="YES"){ 
          RSS <- sum((y.in - (as.numeric(scv.cali[[i]]) * pred.y) )^2) 
        } else{
          RSS <- sum((y.in - pred.y)^2)  
        }
    
    for (n in 1: nb.tree)
        { 
          if (n==1){ nb.weight <- sum(bag.tr$mtrees[[n]]$btree$frame$wt!=0) }
            else{ nb.weight <- c( nb.weight, sum(bag.tr$mtrees[[n]]$btree$frame$wt!=0) )}
        }
    nb.weight <- mean(nb.weight)
    
    aic.temp <- 2*nb.weight - length(y.in)*log(RSS/length(y.in)) # AIC
    
    aicc.temp <- aic.temp + (2*nb.weight +(nb.weight+1) )/
      (length(y.in) -  nb.weight-1) #AICc
    
    if (CR.VAL=="YES"){
    modl.perf.bag <- rbind(modl.perf.bag, 
                           c(nb.perf=nrow(x.in), M.perf,
                             nb.vali=nrow(x.vali),V.perf,CV.perf, 
                             RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
    }else{
      modl.perf.bag <- rbind(modl.perf.bag, 
                            c(nb.perf=nrow(x.in), M.perf,
                              nb.vali=nrow(x.vali),V.perf,
                              RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
    }
    
  } # end of loop repetition
  
  # write all model perf run
  write.table(modl.perf.bag, file=paste(outpath,"/",mod.tech,"_model_perf_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  #write table for all run
  write.table(varimp.bag, file=paste(outpath,"/",mod.tech,"_var_imp_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  
  # calculate avearge values for model perf and variable importance
  modl.perf.av <- rbind( modl.perf.av, apply(modl.perf.bag,2,mean) )
  varimp.av <- rbind( varimp.av, apply(varimp.bag,2,mean) )
  
  # name technique update
  mod.name.perf  <- rbind(mod.name.perf, c(species= var.y[tu],var.sel= var.sele[zu],technic=mod.tech) )
  mod.name <- rbind(mod.name, c(species= var.y[tu],var.sel= var.sele[zu],technic=mod.tech) )
  
} # end loop modeling technique

#####################################################################################################################################
if (EGB == "NO") {}else{
  mod.tech <- "EGB"
  outpath <- paste(output,"/",mod.tech,sep="")
  creat.subDir (output,mod.tech)  
  
  for (i in 1:n.rep)    
  {
    # creat x and y table
    #####################
    if (SVC=="YES")
    {
      y.in <- data.frame(d.cali[[i]][,names(d.cali[[i]])==var.y[tu]])
      y.in <- data.frame(y.in[scv.cali[[i]]==1,1])
      names(y.in) <- var.y[tu]
      x.in <- d.cali[[i]][, na.omit( match(names(d.cali[[i]]), 
                                           var.list$var.name[var.list$var.type=="x"] ) ) ]
      x.in <- x.in[scv.cali[[i]]==1,]
      data.in <- cbind(x.in,y.in)
      
      x.vali <- d.vali[[i]][, na.omit( match(names(d.vali[[i]]), 
                                             var.list$var.name[var.list$var.type=="x"] ) ) ]
      x.vali <- x.vali[scv.vali[[i]]==1,]
      y.vali <- data.frame(d.vali[[i]][,names(d.vali[[i]])==var.y[tu]])
      y.vali <- data.frame(y.vali[scv.vali[[i]]==1,1])
      
      p.cal.svc <- as.numeric(row.names(d.cali[[i]]))[scv.cali[[i]]==1]
      p.val.svc <- as.numeric(row.names(d.vali[[i]]))[scv.vali[[i]]==1]
    }else{
      y.in <- data.frame(d.cali[[i]][,names(d.cali[[i]])==var.y[tu]])
      names(y.in) <- var.y[tu]
      x.in <- d.cali[[i]][, na.omit( match(names(d.cali[[i]]), 
                                           var.list$var.name[var.list$var.type=="x"] ) ) ]
      data.in <- cbind(x.in,y.in)
      
      x.vali <- d.vali[[i]][, na.omit( match(names(d.vali[[i]]), 
                                             var.list$var.name[var.list$var.type=="x"] ) ) ]
      y.vali <- data.frame(d.vali[[i]][,names(d.vali[[i]])==var.y[tu]]) 
      
    }
    #####################
    ## EGB -- START --> get info on http://xgboost.readthedocs.io/en/latest/get_started
    #####################
    cat("> Extreme Gradient Boosting run ",i,"sampling  ...", "\n",append = FALSE)
    
    # set the parameter of the booster  (default)
    para <- list(booster= "gbtree", max_depth = 6, eta = 0.3, silent = 1, nthread = 2,
                 objective = "reg:linear", eval_metric = "rmse")
    
    # run model
    egb.model <- xgboost(data = as.matrix(x.in), label = row.names(x.in), #t(as.vector(y.in)), missing = NA, weight = NULL,
                params = para , nrounds = it.max, verbose = 0 )
            
    # Save Model - R file
    file.name <- paste(outpath,"/",mod.tech,"_model_",i,sep='')
    xgb.save(egb.model, file.name) 
    
    ######################
    ## Variable importance
    ######################
    varimp <- xgb.importance(colnames(x.in), model = egb.model)
    
    # weird but need the plot to calculate the importance overwise don't do
    xgb.plot.importance (varimp, rel_to_first = TRUE, xlab = "Relative importance")
    
    #re-organized to have alphabetic order
    order <- varimp$Feature
    sort <- sort(order, index.return=TRUE)
    varimp <- varimp$Importance[sort$ix]
    
    varimp.egb <- rbind(varimp.egb, varimp ) 
    
    ## MODEL PERFORMANCE
    ####################
    # PREDICT THE OBS
    pred.y <- predict(egb.model, as.matrix(x.in)) 
    pred.y.vali <- predict(egb.model, as.matrix(x.vali))
    
    ## recombin SVC classification and ML prediction
    ################################################
    if (SVC=="YES"){ 
      pred.y <- c(pred.y, rep(0,nrow(d.cali[[i]])-length(pred.y)) )
      pred.y.vali <- c(pred.y.vali, rep(0,nrow(d.vali[[i]])-length(pred.y.vali)) ) 
      
      y.in <- c(y.in[,1], d.cali[[i]][ scv.cali[[i]]==0,names(d.cali[[i]])==var.y[tu] ] )
      y.vali <- c(y.vali[,1], d.vali[[i]][ scv.vali[[i]]==0,names(d.vali[[i]])==var.y[tu] ] )
      
      p.cali <- as.numeric(row.names(x.in))
      p.vali <- as.numeric(row.names(x.vali))
      
      if (nrow(d.cali[[i]])==length(p.cal.svc)){}else{ 
        p.cali <- c(p.cali,
                    as.numeric(row.names(d.cali[[i]]))[!as.numeric(row.names(d.cali[[i]]))%in%p.cali]) }
      
      if (nrow(d.vali[[i]])==length(p.val.svc)){}else{ 
        p.vali <- c(p.vali,
                    as.numeric(row.names(d.vali[[i]]))[!as.numeric(row.names(d.vali[[i]]))%in%p.vali])}
      
    } else{
      y.in <- y.in[,1]
      y.vali <- y.vali[,1]
      
      p.cali <- as.numeric(row.names(x.in))
      p.vali <- as.numeric(row.names(x.vali))
    }
    
    # PERFORMANCE OF THE MODEL
    #	lin.corr <- lm(y.in ~ pred.y) ;summary(lin.corr)
    M.perf <- asses.lm.(pred.y, y.in )
    V.perf <- asses.lm.(pred.y.vali, y.vali )
    
    names(M.perf) <- paste(rep("M_",5),names(M.perf),sep="")
    names(V.perf) <- paste(rep("V_",5),names(V.perf),sep="")
    
    ##keep data and rescale it
    if (RES=="YES"){
      m.temp <- RES.data[,ncol(RES.data)][1]
      sd.temp <- RES.data[,ncol(RES.data)][2]
      
      pred.y <- (pred.y*sd.temp) + m.temp
      y.in <- (y.in*sd.temp) + m.temp
      pred.y.vali <- (pred.y.vali*sd.temp) + m.temp
      y.vali <- (y.vali*sd.temp) + m.temp
      
    } else{}
    
    ## keep data and retransformed it...
    if (TRANSF=="YES"){
        f <- as.function(inv.transf[[var.list$var.transf[var.list$var.name==var.y[tu]]]])
        cali <- data.frame(cbind(p.cali, rep("CALI",length(pred.y)), f(x=y.in), f(x=pred.y) ))
        vali <- data.frame(cbind(p.vali, rep("VALI",length(pred.y.vali)), f(x=y.vali), f(x=pred.y.vali) ))
      }else{
        cali <- data.frame(cbind(p.cali, rep("CALI",length(pred.y)), y.in, pred.y ))
        vali <- data.frame(cbind(p.vali, rep("VALI",length(pred.y.vali)), y.vali, pred.y.vali ))
      }
    
    names(cali) <-c("id","Type","obs.y","pred.y")
    names(vali) <- c("id","Type","obs.y","pred.y")
    
    cali.vali <- rbind(cali,vali)
    
    write.table(cali.vali, file=paste(outpath,"/",mod.tech,"_cali.vali_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    ###CROSS - Validate model
    #########################
    if (CR.VAL=="YES"){
    CV.data  <- CV.EGB.(data.in,  K=fold.cv, cv.lim = 10, it.max= it.max, name.sp = var.y[tu])  
    
    ## OBS = fct(PRED) - CV set -- follow Pineiro 2008
    ##------------------------
    CV.perf <- asses.lm.(CV.data$predictions,CV.data$obs)
    
    # re-rescale the prediction
    if (RES=="YES"){
      m.temp <- RES.data[,ncol(RES.data)][1]
      sd.temp <- RES.data[,ncol(RES.data)][2]
      
      CV.data$predictions <- (CV.data$predictions*sd.temp) + m.temp
      CV.data$obs <- (CV.data$obs *sd.temp ) + m.temp
    }else{}
    
    # re-transform CV data..
    if (TRANSF=="YES"){
      CV.data$predictions <- f(x=CV.data$predictions)
      CV.data$obs <- f(x=CV.data$obs)
    }else{}
    
    write.table(CV.data, file=paste(outpath,"/",mod.tech,"_CV_data_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    names(CV.perf) <- paste(rep("CV_",5),names(CV.perf),sep="")
    
    if (i ==1) {
      png(paste(outpath,"/",mod.tech,"_OBSvsPRED_rep_",i,".png",sep=""),width = 30, 
          height = 15,units="cm",res=150)
      
      ## Plot observed VS predicted for CV
      plot (CV.data$predictions, CV.data$obs, xlab="Predict_CV", ylab="Observed_CV")
      
      dev.off()
    }else{} 
    } else{}
    
    ## AIC and AICC of the Best model
    ################################
    if (SVC=="YES"){ 
          RSS <- sum((y.in - (as.numeric(scv.cali[[i]]) *pred.y ))^2) 
        } else{
          RSS <- sum((y.in - pred.y)^2)  
        }
    
    aic.temp <- NA # AIC
    aicc.temp <- NA #AICc
    
    if (CR.VAL=="YES"){
    modl.perf.egb <- rbind(modl.perf.egb, 
                           c(nb.perf=nrow(x.in), M.perf,
                             nb.vali=nrow(x.vali),V.perf,CV.perf, 
                             RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
    }else{
      modl.perf.egb <- rbind(modl.perf.egb, 
                             c(nb.perf=nrow(x.in), M.perf,
                               nb.vali=nrow(x.vali),V.perf,
                               RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
      
    }
    
  } # end of loop repetition
  
  # write all model perf run
  write.table(modl.perf.egb, file=paste(outpath,"/",mod.tech,"_model_perf_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  
  #write table for all run
  write.table(varimp.egb, file=paste(outpath,"/",mod.tech,"_var_imp_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  
  # calculate avearge values for model perf and variable importance
  modl.perf.av <- rbind( modl.perf.av, apply(modl.perf.egb,2,mean) )
  varimp.av <- rbind( varimp.av, apply(varimp.egb,2,mean) )
  
  # name technique update
  mod.name.perf  <- rbind(mod.name.perf, c(species= var.y[tu],var.sel= var.sele[zu],technic=mod.tech) )
  mod.name <- rbind(mod.name, c(species= var.y[tu],var.sel= var.sele[zu],technic=mod.tech) )
  
} # end loop modeling technique

#####################################################################################################################################
if (GBM == "NO") {}else{
  
    mod.tech <- "GBM"
    outpath <- paste(output,"/",mod.tech,sep="")
    creat.subDir (output,mod.tech)
    
  for (i in 1:n.rep)    
  {
    # creat x and y table
    #####################
    if (SVC=="YES")
    { 
        y.in <- data.frame(d.cali[[i]][,names(d.cali[[i]])==var.y[tu]])
        y.in <- data.frame(y.in[scv.cali[[i]]==1,1])
        names(y.in) <- var.y[tu]
        x.in <- d.cali[[i]][, na.omit( match(names(d.cali[[i]]), 
                                             var.list$var.name[var.list$var.type=="x"] ) ) ]
        x.in <- x.in[scv.cali[[i]]==1,]
        data.in <- cbind(x.in,y.in)
        
        x.vali <- d.vali[[i]][, na.omit( match(names(d.vali[[i]]), 
                                               var.list$var.name[var.list$var.type=="x"] ) ) ]
        x.vali <- x.vali[scv.vali[[i]]==1,]
        y.vali <- data.frame(d.vali[[i]][,names(d.vali[[i]])==var.y[tu]])
        y.vali <- data.frame(y.vali[scv.vali[[i]]==1,1])
        
        p.cal.svc <- as.numeric(row.names(d.cali[[i]]))[scv.cali[[i]]==1]
        p.val.svc <- as.numeric(row.names(d.vali[[i]]))[scv.vali[[i]]==1]
      }else{
        y.in <- data.frame(d.cali[[i]][,names(d.cali[[i]])==var.y[tu]])
        names(y.in) <- var.y[tu]
        x.in <- d.cali[[i]][, na.omit( match(names(d.cali[[i]]), 
                                             var.list$var.name[var.list$var.type=="x"] ) ) ]
        data.in <- cbind(x.in,y.in)
        
        x.vali <- d.vali[[i]][, na.omit( match(names(d.vali[[i]]), 
                                               var.list$var.name[var.list$var.type=="x"] ) ) ]
        y.vali <- data.frame(d.vali[[i]][,names(d.vali[[i]])==var.y[tu]])
      }
    ## OPTIMIZED PARAMETER PERFORMANCE OF GBM
    ##########################################
    pred <- var.list$var.name[var.list$var.type=="x"]
    form <- as.formula(paste(colnames(y.in), "~", paste(pred[!pred %in% colnames(y.in)],
                                                        collapse = " + ")))
    
    cat("> GBM MODEL run ",i,"sampling  ...", "\n",append = FALSE)
    
    gbm.tr <- gbm(form, data=as.data.frame(data.in), n.trees=nb.tree, distribution = "gaussian" )
    
    # Save Model - R file
    file.name<-paste(outpath,"/",mod.tech,"_model_",i,sep='')
    save(list=c('gbm.tr'),file=file.name)
    
    ######################
    ## Variable importance -- permutation method
    ######################
    gbm.imp <- relative.influence(gbm.tr, n.trees= nb.tree, scale.=FALSE, sort.=FALSE)

    varimp.gbm <- rbind(varimp.gbm,(gbm.imp /sum(gbm.imp)) ) # relative importance to 1
    
    ## MODEL PERFORMANCE
    ####################
    # PREDICT THE OBS
    pred.y <- predict.gbm(gbm.tr, newdata=as.data.frame(x.in), 
                        n.trees=nb.tree, type="response",na.rm = TRUE)
    pred.y.vali <- predict.gbm(gbm.tr, newdata=as.data.frame(x.vali), 
                               n.trees=nb.tree, type="response",na.rm = TRUE)
    ## recombin SVC classification and ML prediction
    ################################################
    if (SVC=="YES"){ 
      pred.y <- c(pred.y, rep(0,nrow(d.cali[[i]])-length(pred.y)) )
      pred.y.vali <- c(pred.y.vali,  rep(0,nrow(d.vali[[i]])-length(pred.y.vali)) ) 
      
      y.in <- c(y.in[,1], d.cali[[i]][ scv.cali[[i]]==0,names(d.cali[[i]])==var.y[tu] ] )
      y.vali <- c(y.vali[,1], d.vali[[i]][ scv.vali[[i]]==0,names(d.vali[[i]])==var.y[tu] ] )
      
      p.cali <- as.numeric(row.names(x.in))
      p.vali <- as.numeric(row.names(x.vali))
      
      if (nrow(d.cali[[i]])==length(p.cal.svc)){}else{ 
        p.cali <- c(p.cali,
                    as.numeric(row.names(d.cali[[i]]))[!as.numeric(row.names(d.cali[[i]]))%in%p.cali]) }
      
      if (nrow(d.vali[[i]])==length(p.val.svc)){}else{ 
        p.vali <- c(p.vali,
                    as.numeric(row.names(d.vali[[i]]))[!as.numeric(row.names(d.vali[[i]]))%in%p.vali])}
      
      
    } else{
      y.in <- y.in[,1]
      y.vali <- y.vali[,1]
      
      p.cali <- as.numeric(row.names(x.in))
      p.vali <- as.numeric(row.names(x.vali))
    }
    
    # PERFORMANCE OF THE MODEL
    #	lin.corr <- lm(y.in ~ pred.y) ;summary(lin.corr)
    M.perf <- asses.lm.(pred.y, y.in )
    V.perf <- asses.lm.(pred.y.vali, y.vali )
    
    names(M.perf) <- paste(rep("M_",5),names(M.perf),sep="")
    names(V.perf) <- paste(rep("V_",5),names(V.perf),sep="")
    
    ##keep data and rescale it
    if (RES=="YES"){
      m.temp <- RES.data[,ncol(RES.data)][1]
      sd.temp <- RES.data[,ncol(RES.data)][2]
      
      pred.y <- (pred.y*sd.temp) + m.temp
      y.in <- (y.in*sd.temp) + m.temp
      pred.y.vali <- (pred.y.vali*sd.temp) + m.temp
      y.vali <- (y.vali*sd.temp) + m.temp
      
    } else{}
    
    ## keep data and retransformed it...
    if (TRANSF=="YES"){
        f <- as.function(inv.transf[[var.list$var.transf[var.list$var.name==var.y[tu]]]])
        cali <- data.frame(cbind(p.cali, rep("CALI",length(pred.y)), f(x=y.in), f(x=pred.y) ))
        vali <- data.frame(cbind(p.vali, rep("VALI",length(pred.y.vali)), f(x=y.vali), f(x=pred.y.vali) ))
      }else{
        cali <- data.frame(cbind(p.cali, rep("CALI",length(pred.y)), y.in, pred.y ))
        vali <- data.frame(cbind(p.vali, rep("VALI",length(pred.y.vali)), y.vali, pred.y.vali ))
      }
    
    names(cali) <-c("id","Type","obs.y","pred.y")
    names(vali) <- c("id", "Type","obs.y","pred.y")  	
    
    cali.vali <- rbind(cali,vali)
    
    write.table(cali.vali, file=paste(outpath,"/",mod.tech,"_cali.vali_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    ## CROSS - Validate model
    #########################
    if (CR.VAL=="YES"){
    CV.data  <- CV.GBM.(gbm.tr, data.in, nb.tree= nb.tree, K=fold.cv, cv.lim = 10, 
                        name.sp = var.y[tu])
    
    ## OBS = fct(PRED) - CV set -- follow Pineiro 2008
    ##------------------------
    CV.perf <- asses.lm.(CV.data$predictions,CV.data$obs)
    
    # re-rescale the prediction
    if (RES=="YES"){
      m.temp <- RES.data[,ncol(RES.data)][1]
      sd.temp <- RES.data[,ncol(RES.data)][2]
      
      CV.data$predictions <- (CV.data$predictions*sd.temp) + m.temp
      CV.data$obs <- (CV.data$obs *sd.temp ) + m.temp
    }else{}
    
    # re-transform CV data..
    if (TRANSF=="YES"){
      CV.data$predictions <- f(x=CV.data$predictions)
      CV.data$obs <- f(x=CV.data$obs)
    }else{}
    
    write.table(CV.data, file=paste(outpath,"/",mod.tech,"_CV_data_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    names(CV.perf) <- paste(rep("CV_",5),names(CV.perf),sep="")
    
    if (i==1) {
      png(paste(outpath,"/",mod.tech,"_OBSvsPRED_rep_",i,".png",sep=""),width = 30, height = 15,units="cm",res=150)
      
      ## Plot observed VS predicted for CV
      plot (CV.data$predictions, CV.data$obs, xlab="Predict_CV", ylab="Observed_CV")
      
      dev.off()
      
    } else {}
    
    }else{}
    
    ## AIC and AICC of the Best model
    ################################
    if (SVC=="YES"){ 
        RSS <- sum((y.in - (as.numeric(scv.cali[[i]]) *pred.y ))^2) 
      } else{
        RSS <- sum((y.in - pred.y)^2) 
      }
     
    for (n in 1: nb.tree)
      { 
        if (n==1){ nb.weight <- sum(unlist(gbm.tr$trees[[n]])!=0) }
        else{ nb.weight <- c( nb.weight, sum(unlist(gbm.tr$trees[[n]])!=0) )}
      }
    nb.weight <- mean(nb.weight)
   
    aic.temp <- 2*nb.weight - length(y.in)*log(RSS/length(y.in)) # AIC
    
    aicc.temp <- aic.temp + (2*nb.weight +(nb.weight+1) )/
      (length(y.in) -  nb.weight-1) #AICc
    
    if (CR.VAL=="YES"){
    modl.perf.gbm <- rbind(modl.perf.gbm, 
                           c(nb.perf=nrow(x.in), M.perf,
                             nb.vali=nrow(x.vali),V.perf,CV.perf, 
                             RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
    }else{
      modl.perf.gbm <- rbind(modl.perf.gbm, 
                             c(nb.perf=nrow(x.in), M.perf,
                               nb.vali=nrow(x.vali),V.perf,
                               RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
      
    }
  } # end of loop repetition

  # write all model perf run
  write.table(modl.perf.gbm, file=paste(outpath,"/",mod.tech,"_model_perf_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  
  #write table for all run
  write.table(varimp.gbm, file=paste(outpath,"/",mod.tech,"_var_imp_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  
  # calculate avearge values for model perf and variable importance
  modl.perf.av <- rbind( modl.perf.av, apply(modl.perf.gbm,2,mean) )
  varimp.av <- rbind( varimp.av, apply(varimp.gbm,2,mean) )
  
  # name technique update
  mod.name.perf  <- rbind(mod.name.perf, c(species= var.y[tu],var.sel= var.sele[zu],technic=mod.tech) )
  mod.name <- rbind(mod.name, c(species= var.y[tu],var.sel= var.sele[zu],technic=mod.tech) )
  
} # end loop modeling technique
  
##################################################################################################################################### 
if (GLM == "NO") {}else{
    mod.tech <- "GLM"
    outpath <- paste(output,"/",mod.tech,sep="")
    creat.subDir (output,mod.tech) 
  
    for (i in 1:n.rep)    
    {
      # creat x and y table
      #####################
      if (SVC=="YES")
      { 
        y.in <- data.frame(d.cali[[i]][,names(d.cali[[i]])==var.y[tu]])
        y.in <- data.frame(y.in[scv.cali[[i]]==1,1])
        names(y.in) <- var.y[tu]
        
        x.in <- d.cali[[i]][, na.omit( match(names(d.cali[[i]]), 
                                             var.list$var.name[var.list$var.type=="x"] ) ) ]
        x.in <- x.in[scv.cali[[i]]==1,]
        data.in <- cbind(x.in,y.in)
        
        x.vali <- d.vali[[i]][, na.omit( match(names(d.vali[[i]]), 
                                               var.list$var.name[var.list$var.type=="x"] ) ) ]
        x.vali <- x.vali[scv.vali[[i]]==1,]
        y.vali <- data.frame(d.vali[[i]][,names(d.vali[[i]])==var.y[tu]])
        y.vali <- data.frame(y.vali[scv.vali[[i]]==1,1])
        
        data.vali <- cbind(y.vali,x.vali)
        
        p.cal.svc <- as.numeric(row.names(d.cali[[i]]))[scv.cali[[i]]==1]
        p.val.svc <- as.numeric(row.names(d.vali[[i]]))[scv.vali[[i]]==1]
        
        }else{
        
        y.in <- data.frame(d.cali[[i]][,names(d.cali[[i]])==var.y[tu]])
        names(y.in) <- var.y[tu]
        x.in <- d.cali[[i]][, na.omit( match(names(d.cali[[i]]), 
                                             var.list$var.name[var.list$var.type=="x"] ) ) ]
        data.in <- cbind(x.in,y.in)
        
        x.vali <- d.vali[[i]][, na.omit( match(names(d.vali[[i]]), 
                                               var.list$var.name[var.list$var.type=="x"] ) ) ]
        y.vali <- data.frame(d.vali[[i]][,names(d.vali[[i]])==var.y[tu]])
        
        data.vali <- cbind(y.vali,x.vali) 
      }
      
    #####################
    ## GLM - Multivariate
    #####################
	  cat("> GLM MODEL run ",i,"sampling  ...", "\n",append = FALSE)
	  
    data.cal <- data.frame(data.in)    
    
    names.pred <- names(x.in)
    
    name.sp <- colnames(y.in)
    
    for (p in 1:length(names.pred)) 
        { 
          if (p == 1) 
          {
            fmula.pred.glm <- as.vector(paste("poly(",names.pred[p],",2)",sep=""),mode="character")  
          }
          else  
          {
            fmula.pred.glm <- paste(fmula.pred.glm," + poly(",names.pred[p],",2)",sep="") 
          }
        }		 	
    
    ###########
    # GLM FIT #
    ###########
    glm.tmp <- glm(eval(parse(text = paste(paste(name.sp), "~",fmula.pred.glm, collapse = ""))),
                             data=data.cal,family=gaussian,maxit = 100)    
  
    # Save Model - R file
    file.name <- paste(outpath,"/",mod.tech,"_model_",i,sep='')
    save(list=c('glm.tmp'),file=file.name) 
    
    ######################
    ## Variable importance
    ######################
    varimp <- var.imp.glm.step(glm.tmp,data.cal,names.pred,nperm = n.perm)
    
    varimp.glm <- rbind(varimp.glm, (varimp /sum(varimp)) )
  
    ## MODEL PERFORMANCE
    ####################
    # PREDICT THE OBS
    pred.y <- predict(glm.tmp, newdata = data.cal, type="response",na.rm = TRUE)
    pred.y.vali <- predict(glm.tmp, newdata = data.frame(data.vali), 
                           type="response",na.rm = TRUE)
    ## recombin SVC classification and ML prediction
    ################################################
    if (SVC=="YES"){ 
      pred.y <- c(pred.y, rep(0,nrow(d.cali[[i]])-length(pred.y)) ) 
      pred.y.vali <- c(pred.y.vali, rep(0,nrow(d.vali[[i]])-length(pred.y.vali)) ) 
      
      y.in <- c(y.in[,1], d.cali[[i]][ scv.cali[[i]]==0,names(d.cali[[i]])==var.y[tu] ] )
      y.vali <- c(y.vali[,1], d.vali[[i]][ scv.vali[[i]]==0,names(d.vali[[i]])==var.y[tu] ] )
      
      p.cali <- as.numeric(row.names(x.in))
      p.vali <- as.numeric(row.names(x.vali))
      
      if (nrow(d.cali[[i]])==length(p.cal.svc)){}else{ 
        p.cali <- c(p.cali,
                    as.numeric(row.names(d.cali[[i]]))[!as.numeric(row.names(d.cali[[i]]))%in%p.cali]) }
      
      if (nrow(d.vali[[i]])==length(p.val.svc)){}else{ 
        p.vali <- c(p.vali,
                    as.numeric(row.names(d.vali[[i]]))[!as.numeric(row.names(d.vali[[i]]))%in%p.vali])}
      
      
    } else{
      y.in <- y.in[,1]
      y.vali <- y.vali[,1]
      
      p.cali <- as.numeric(row.names(x.in))
      p.vali <- as.numeric(row.names(x.vali))
    }
    
    # PERFORMANCE OF THE MODEL
    #	lin.corr <- lm(y.in ~ pred.y) ;summary(lin.corr)
    M.perf <- asses.lm.(pred.y, y.in )
    V.perf <- asses.lm.(pred.y.vali, y.vali )
    
    names(M.perf) <- paste(rep("M_",5),names(M.perf),sep="")
    names(V.perf) <- paste(rep("V_",5),names(V.perf),sep="")
    
    ##keep data and rescale it
    if (RES=="YES"){
      m.temp <- RES.data[,ncol(RES.data)][1]
      sd.temp <- RES.data[,ncol(RES.data)][2]
      
      pred.y <- (pred.y*sd.temp) + m.temp
      y.in <- (y.in*sd.temp) + m.temp
      pred.y.vali <- (pred.y.vali*sd.temp) + m.temp
      y.vali <- (y.vali*sd.temp) + m.temp
      
    } else{}
    
    ## keep data and retransformed it...
    if (TRANSF=="YES"){
        f <- as.function(inv.transf[[var.list$var.transf[var.list$var.name==var.y[tu]]]])
        cali <- data.frame(cbind(p.cali, rep("CALI",length(pred.y)), f(x=y.in), f(x=pred.y) ))
        vali <- data.frame(cbind(p.vali, rep("VALI",length(pred.y.vali)), f(x=y.vali), f(x=pred.y.vali) ))
      }else{
        cali <- data.frame(cbind(p.cali, rep("CALI",length(pred.y)), y.in, pred.y ))
        vali <- data.frame(cbind(p.vali, rep("VALI",length(pred.y.vali)), y.vali, pred.y.vali ))
      }
    
    names(cali) <-c("id", "Type","obs.y","pred.y")
    names(vali) <- c("id","Type","obs.y","pred.y")
    
    cali.vali <- rbind(cali,vali)
    
    write.table(cali.vali, file=paste(outpath,"/",mod.tech,"_cali.vali_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    ## CROSS - Validate model
    #########################
    if (CR.VAL=="YES"){
    CV.data  <- CV.glm.(glm.tmp, K=fold.cv, cv.lim = 10, name.sp = var.y[tu])  
    
    ## OBS = fct(PRED) - CV set -- follow Pineiro 2008
    ##------------------------
    CV.perf <- asses.lm.(CV.data$predictions,CV.data$obs)
    
    # re-rescale the prediction
    if (RES=="YES"){
      m.temp <- RES.data[,ncol(RES.data)][1]
      sd.temp <- RES.data[,ncol(RES.data)][2]
      
      CV.data$predictions <- (CV.data$predictions*sd.temp) + m.temp
      CV.data$obs <- (CV.data$obs *sd.temp ) + m.temp
    }else{}
    
    # re-transform CV data..
    if (TRANSF=="YES"){
      CV.data$predictions <- f(x=CV.data$predictions)
      CV.data$obs <- f(x=CV.data$obs)
    }else{}
    
    write.table(CV.data, file=paste(outpath,"/",mod.tech,"_CV_data_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    names(CV.perf) <- paste(rep("CV_",5),names(CV.perf),sep="")
    
    if (i ==1) {
      png(paste(outpath,"/",mod.tech,"_OBSvsPRED_rep_",i,".png",sep=""),width = 30, 
          height = 15,units="cm",res=150)
      
        ## Plot observed VS predicted for CV
        plot (CV.data$predictions, CV.data$obs, xlab="Predict_CV", ylab="Observed_CV")
      
      dev.off()
    }else{}
    }else{}
    
    ## AIC and AICC of the Best model
    ################################
    if (SVC=="YES"){ 
          RSS <- sum((y.in - (as.numeric(scv.cali[[i]]) *pred.y ))^2) 
        } else{
          RSS <- sum((y.in - pred.y)^2) 
        }
    
    aic.temp <- AIC(glm.tmp) # AIC
    aicc.temp <- AICc(glm.tmp) #AICc
    
    if (CR.VAL=="YES"){
    modl.perf.glm <- rbind(modl.perf.glm, 
                           c(nb.perf=nrow(x.in), M.perf,
                             nb.vali=nrow(x.vali),V.perf,CV.perf, 
                             RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
    }else{
      modl.perf.glm <- rbind(modl.perf.glm, 
                             c(nb.perf=nrow(x.in), M.perf,
                               nb.vali=nrow(x.vali),V.perf, 
                               RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
      
    }
    
    } # end of loop repetition
  
  # write all model perf run
  write.table(modl.perf.glm, file=paste(outpath,"/",mod.tech,"_model_perf_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  
  #write table for all run
  write.table(varimp.glm, file=paste(outpath,"/",mod.tech,"_var_imp_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  
  # calculate avearge values for model perf and variable importance
  modl.perf.av <- rbind( modl.perf.av, apply(modl.perf.glm,2,mean) )
  varimp.av <- rbind( varimp.av, apply(varimp.glm,2,mean) )
  
  # name technique update
  mod.name.perf  <- rbind(mod.name.perf, c(species= var.y[tu],var.sel= var.sele[zu],technic=mod.tech) )
  mod.name <- rbind(mod.name, c(species= var.y[tu],var.sel= var.sele[zu],technic=mod.tech) )
  
} # end loop modeling technique  
  
##################################################################################################################################### 
if (GLM_sw == "NO") {}else{
  
    mod.tech <- "GLM_sw"
    outpath <- paste(output,"/",mod.tech,sep="")
    creat.subDir (output,mod.tech) 
    
    for (i in 1:n.rep)    
    {
      # creat x and y table
      #####################
      if (SVC=="YES")
      {
        y.in <- data.frame(d.cali[[i]][,names(d.cali[[i]])==var.y[tu]])
        y.in <- data.frame(y.in[scv.cali[[i]]==1,1])
        names(y.in) <- var.y[tu]
        
        x.in <- d.cali[[i]][, na.omit( match(names(d.cali[[i]]), 
                                             var.list$var.name[var.list$var.type=="x"] ) ) ]
        x.in <- x.in[scv.cali[[i]]==1,]
        data.in <- cbind(y.in,x.in)
        
        x.vali <- d.vali[[i]][, na.omit( match(names(d.vali[[i]]), 
                                               var.list$var.name[var.list$var.type=="x"] ) ) ]
        x.vali <- x.vali[scv.vali[[i]]==1,]
        y.vali <- data.frame(d.vali[[i]][,names(d.vali[[i]])==var.y[tu]])
        y.vali <- data.frame(y.vali[scv.vali[[i]]==1,1])
        
        data.vali <- cbind(y.vali,x.vali)
        
        p.cal.svc <- as.numeric(row.names(d.cali[[i]]))[scv.cali[[i]]==1]
        p.val.svc <- as.numeric(row.names(d.vali[[i]]))[scv.vali[[i]]==1]
        
      } else{
        y.in <- data.frame(d.cali[[i]][,names(d.cali[[i]])==var.y[tu]])
        names(y.in) <- var.y[tu]
        x.in <- d.cali[[i]][, na.omit( match(names(d.cali[[i]]), 
                                             var.list$var.name[var.list$var.type=="x"] ) ) ]
        data.in <- cbind(y.in,x.in)
        
        x.vali <- d.vali[[i]][, na.omit( match(names(d.vali[[i]]), 
                                               var.list$var.name[var.list$var.type=="x"] ) ) ]
        y.vali <- data.frame(d.vali[[i]][,names(d.vali[[i]])==var.y[tu]])
        
        data.vali <- cbind(y.vali,x.vali)
        
      }
    
    ###############
    ## GLM -STEP ##
    ###############
    data.cal <- data.frame(data.in)    
    
    names.pred <- names(x.in)
    
    name.sp <- colnames(y.in)
    
    glm.tmp.step <- step(glm(eval(parse(text = paste(paste(name.sp), "~",fmula.pred.glm, collapse = ""))),
                             data=data.cal,family=gaussian,maxit = 100),trace = F,direction="both")
    # Save Model - R file
    file.name <- paste(outpath,"/",mod.tech,"_model_",i,sep='')
    save(list=c('glm.tmp.step'),file=file.name) 
    
    ######################
    ## Variable importance
    ######################
    varimp <- var.imp.glm.step(glm.tmp.step,data.cal,names.pred,nperm = n.perm)
    
    varimp.glmstep <- rbind(varimp.glmstep, (varimp/sum(varimp)) )
    
    ####################
    ## MODEL PERFORMANCE
    ####################
    # PREDICT THE OBS
    pred.y <- predict(glm.tmp.step, newdata = data.cal, type="response",na.rm = TRUE)
    pred.y.vali <- predict(glm.tmp.step, newdata = data.frame(data.vali), 
                           type="response",na.rm = TRUE)
    
    ## recombin SVC classification and ML prediction
    ################################################
    if (SVC=="YES"){ 
      pred.y <- c(pred.y, rep(0,nrow(d.cali[[i]])-length(pred.y)) )
      pred.y.vali <- c(pred.y.vali, rep(0,nrow(d.vali[[i]])-length(pred.y.vali)) )
      
      y.in <- c(y.in[,1], d.cali[[i]][ scv.cali[[i]]==0,names(d.cali[[i]])==var.y[tu] ] )
      y.vali <- c(y.vali[,1], d.vali[[i]][ scv.vali[[i]]==0,names(d.vali[[i]])==var.y[tu] ] )
      
      p.cali <- as.numeric(row.names(x.in))
      p.vali <- as.numeric(row.names(x.vali))
      
      if (nrow(d.cali[[i]])==length(p.cal.svc)){}else{ 
        p.cali <- c(p.cali,
                    as.numeric(row.names(d.cali[[i]]))[!as.numeric(row.names(d.cali[[i]]))%in%p.cali]) }
      
      if (nrow(d.vali[[i]])==length(p.val.svc)){}else{ 
        p.vali <- c(p.vali,
                    as.numeric(row.names(d.vali[[i]]))[!as.numeric(row.names(d.vali[[i]]))%in%p.vali])}
      
      
    } else{
      y.in <- y.in[,1]
      y.vali <- y.vali[,1]
      
      p.cali <- as.numeric(row.names(x.in))
      p.vali <- as.numeric(row.names(x.vali))
    }
    
    # PERFORMANCE OF THE MODEL
    #	l#######################
    M.perf <- asses.lm.(pred.y, y.in )
    V.perf <- asses.lm.(pred.y.vali, y.vali )
    
    names(M.perf) <- paste(rep("M_",5),names(M.perf),sep="")
    names(V.perf) <- paste(rep("V_",5),names(V.perf),sep="")
    
    ##keep data and rescale it
    if (RES=="YES"){
      m.temp <- RES.data[,ncol(RES.data)][1]
      sd.temp <- RES.data[,ncol(RES.data)][2]
      
      pred.y <- (pred.y*sd.temp) + m.temp
      y.in <- (y.in*sd.temp) + m.temp
      pred.y.vali <- (pred.y.vali*sd.temp) + m.temp
      y.vali <- (y.vali*sd.temp) + m.temp
      
    } else{}
    
    ## keep data and retransformed it...
    if (TRANSF=="YES"){
        f <- as.function(inv.transf[[var.list$var.transf[var.list$var.name==var.y[tu]]]])
        cali <- data.frame(cbind(p.cali, rep("CALI",length(pred.y)), f(x=y.in), f(x=pred.y) ))
        vali <- data.frame(cbind(p.vali, rep("VALI",length(pred.y.vali)), f(x=y.vali), f(x=pred.y.vali) ))
      }else{
        cali <- data.frame(cbind(p.cali, rep("CALI",length(pred.y)), y.in, pred.y ))
        vali <- data.frame(cbind(p.vali, rep("VALI",length(pred.y.vali)), y.vali, pred.y.vali ))
      }
    
    names(cali) <-c("id", "Type","obs.y","pred.y")
    names(vali) <- c("id","Type","obs.y","pred.y") 
    
    cali.vali <- rbind(cali,vali)
    
    write.table(cali.vali, file=paste(outpath,"/",mod.tech,"_cali.vali_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    ## CROSS - Validate model
    #########################
    if (CR.VAL=="YES"){
    CV.data  <- CV.glm.step(glm.tmp.step, K=fold.cv, cv.lim = 10, name.sp = var.y[tu]) 
    
    ## OBS = fct(PRED) - CV set -- follow Pineiro 2008
    ##------------------------
    CV.perf <- asses.lm.(CV.data$predictions,CV.data$obs)
    
    # re-rescale the prediction
    if (RES=="YES"){
      m.temp <- RES.data[,ncol(RES.data)][1]
      sd.temp <- RES.data[,ncol(RES.data)][2]
      
      CV.data$predictions <- (CV.data$predictions*sd.temp) + m.temp
      CV.data$obs <- (CV.data$obs *sd.temp ) + m.temp
    }else{}
    
    # re-transform CV data..
    if (TRANSF=="YES"){
      CV.data$predictions <- f(x=CV.data$predictions)
      CV.data$obs <- f(x=CV.data$obs)
    }else{}
    
    write.table(CV.data, file=paste(outpath,"/",mod.tech,"_CV_data_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    names(CV.perf) <- paste(rep("CV_",5),names(CV.perf),sep="")
    
    if (i ==1) {
      png(paste(outpath,"/",mod.tech,"_OBSvsPRED_rep_",i,".png",sep=""),width = 30, height = 15,units="cm",res=150)
      
          ## Plot observed VS predicted for CV
          plot (CV.data$predictions, CV.data$obs, xlab="Predict_CV", ylab="Observed_CV")
      
      dev.off()
    }else{}  
    
    }else{}
    
    ## AIC and AICC of the Best model
    ################################
    if (SVC=="YES"){ 
          RSS <- sum((y.in - (as.numeric(scv.cali[[i]]) *pred.y ))^2) 
        } else{
          RSS <- sum((y.in - pred.y)^2) 
        }
    
    aic.temp <- AIC(glm.tmp.step) # AIC
    aicc.temp <- AICc(glm.tmp.step) #AICc
    
    if (CR.VAL=="YES"){
      modl.perf.glmstep <- rbind(modl.perf.glmstep, 
                               c(nb.perf=nrow(x.in), M.perf,
                                 nb.vali=nrow(x.vali),V.perf,CV.perf, 
                                 RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
    }else{
      modl.perf.glmstep <- rbind(modl.perf.glmstep, 
                                 c(nb.perf=nrow(x.in), M.perf,
                                   nb.vali=nrow(x.vali),V.perf, 
                                   RSS=RSS, AIC=aic.temp, AICc=aicc.temp))  
      
    }
    
  } # end of loop repetition

# write all model perf run
write.table(modl.perf.glmstep, file=paste(outpath,"/",mod.tech,"_model_perf_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)

#write table for all run
write.table(varimp.glmstep, file=paste(outpath,"/",mod.tech,"_var_imp_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)

# calculate avearge values for model perf and variable importance
modl.perf.av <- rbind( modl.perf.av, apply(modl.perf.glmstep,2,mean) )
varimp.av <- rbind( varimp.av, apply(varimp.glmstep,2,mean) )

# name technique update
mod.name.perf  <- rbind(mod.name.perf, c(species= var.y[tu],var.sel= var.sele[zu],technic=mod.tech) )
mod.name <- rbind(mod.name, c(species= var.y[tu],var.sel= var.sele[zu],technic=mod.tech) )

} # end loop modeling technique  
  
##################################################################################################################################### 
if (GAM == "NO") {}else{
    mod.tech <- "GAM"
    
    outpath <- paste(output,"/",mod.tech,sep="")
    
    creat.subDir (output,mod.tech) 
    
  for (i in 1:n.rep)    
  {
    # creat x and y table
    #####################
    if (SVC=="YES")
    { 
      y.in <- data.frame(d.cali[[i]][,names(d.cali[[i]])==var.y[tu]])
      y.in <- data.frame(y.in[scv.cali[[i]]==1,1])
      names(y.in) <- var.y[tu]
      
      x.in <- d.cali[[i]][, na.omit( match(names(d.cali[[i]]), 
                                           var.list$var.name[var.list$var.type=="x"] ) ) ]
      x.in <- x.in[scv.cali[[i]]==1,]
      data.in <- cbind(x.in,y.in)
      
      x.vali <- d.vali[[i]][, na.omit( match(names(d.vali[[i]]), 
                                             var.list$var.name[var.list$var.type=="x"] ) ) ]
      x.vali <- x.vali[scv.vali[[i]]==1,]
      y.vali <- data.frame(d.vali[[i]][,names(d.vali[[i]])==var.y[tu]])
      y.vali <- data.frame(y.vali[scv.vali[[i]]==1,1])
      
      data.vali <- cbind(y.vali,x.vali)
      
      p.cal.svc <- as.numeric(row.names(d.cali[[i]]))[scv.cali[[i]]==1]
      p.val.svc <- as.numeric(row.names(d.vali[[i]]))[scv.vali[[i]]==1]
      
      }else{
        
      y.in <- data.frame(d.cali[[i]][,names(d.cali[[i]])==var.y[tu]])
      names(y.in) <- var.y[tu]
      x.in <- d.cali[[i]][, na.omit( match(names(d.cali[[i]]), 
                                         var.list$var.name[var.list$var.type=="x"] ) ) ]
      data.in <- cbind(x.in,y.in)
    
      x.vali <- d.vali[[i]][, na.omit( match(names(d.vali[[i]]), 
                                           var.list$var.name[var.list$var.type=="x"] ) ) ]
      y.vali <- data.frame(d.vali[[i]][,names(d.vali[[i]])==var.y[tu]])
    
      data.vali <- cbind(y.vali,x.vali)
      }
    #####################
    ## GAM - Multivariate
    #####################
    cat("> GAM MODEL run ",i,"sampling  ...", "\n",append = FALSE)
    
    data.cal <- data.frame(data.in)    
    
    names.pred <- names(x.in)
    
    name.sp <- colnames(y.in)
    
    for (p in 1:length(names.pred)) 
      { 
        if (p == 1) 
        {
          fmula.pred.glm <- as.vector(paste("poly(",names.pred[p],",2)",sep=""),mode="character")  
        }
        else  
        {
          fmula.pred.glm <- paste(fmula.pred.glm," + poly(",names.pred[p],",2)",sep="") 
        }
      }		 	
    
    ###########
    # GAM FIT #
    ###########
    gam.tmp <- gam(eval(parse(text = paste(paste(name.sp), "~",fmula.pred.glm, collapse = ""))),
                   data=data.cal,family=gaussian,maxit = 100)    
    
    # Save Model - R file
    file.name <- paste(outpath,"/",mod.tech,"_model_",i,sep='')
    save(list=c('gam.tmp'),file=file.name) 
    
    ######################
    ## Variable importance
    ######################
    varimp <- var.imp.glm.step(gam.tmp,data.cal,names.pred,nperm = n.perm)
      
    varimp.gam <- rbind(varimp.gam, (varimp/sum(varimp)) )
    
    ## MODEL PERFORMANCE
    ####################
    ####################
    # PREDICT THE OBS
    pred.y <- predict(gam.tmp, newdata = data.cal, type="response",na.rm = TRUE)
    
    pred.y.vali <- predict(gam.tmp, newdata = data.frame(data.vali), 
                           type="response",na.rm = TRUE)
    ## recombin SVC classification and ML prediction
    ################################################
    if (SVC=="YES"){ 
      pred.y <- c(pred.y,  rep(0,nrow(d.cali[[i]])-length(pred.y)) ) 
      pred.y.vali <- c(pred.y.vali, rep(0,nrow(d.vali[[i]])-length(pred.y.vali)) ) 
      
      y.in <- c(y.in[,1], d.cali[[i]][ scv.cali[[i]]==0,names(d.cali[[i]])==var.y[tu] ] )
      y.vali <- c(y.vali[,1], d.vali[[i]][ scv.vali[[i]]==0,names(d.vali[[i]])==var.y[tu] ] )
      
      p.cali <- as.numeric(row.names(x.in))
      p.vali <- as.numeric(row.names(x.vali))
      
      if (all(p.cali==p.cal.svc) ==TRUE){}else{ p.cali <- p.cali[!p.cali%in%p.cal.svc]}
      if (all(p.vali==p.val.svc) ==TRUE){}else{ p.vali <- p.vali[!p.vali%in%p.val.svc]}
      
    } else{
      y.in <- y.in[,1]
      y.vali <- y.vali[,1]
      
      p.cali <- as.numeric(row.names(x.in))
      p.vali <- as.numeric(row.names(x.vali))
    }
    
    # PERFORMANCE OF THE MODEL
    ##########################
    M.perf <- asses.lm.(pred.y, y.in )
    V.perf <- asses.lm.(pred.y.vali, y.vali )
    
    names(M.perf) <- paste(rep("M_",5),names(M.perf),sep="")
    names(V.perf) <- paste(rep("V_",5),names(V.perf),sep="")
    
    ##keep data and rescale it
    if (RES=="YES"){
      m.temp <- RES.data[,ncol(RES.data)][1]
      sd.temp <- RES.data[,ncol(RES.data)][2]
      
      pred.y <- (pred.y*sd.temp) + m.temp
      y.in <- (y.in*sd.temp) + m.temp
      pred.y.vali <- (pred.y.vali*sd.temp) + m.temp
      y.vali <- (y.vali*sd.temp) + m.temp
      
    } else{}
    
    ## keep data and retransformed it...
    if (TRANSF=="YES"){
        f <- as.function(inv.transf[[var.list$var.transf[var.list$var.name==var.y[tu]]]])
        cali <- data.frame(cbind(p.cali, rep("CALI",length(pred.y)), f(x=y.in), f(x=pred.y) ))
        vali <- data.frame(cbind(p.vali, rep("VALI",length(pred.y.vali)), f(x=y.vali), f(x=pred.y.vali) ))
      }else{
        cali <- data.frame(cbind(p.cali, rep("CALI",length(pred.y)), y.in, pred.y ))
        vali <- data.frame(cbind(p.vali, rep("VALI",length(pred.y.vali)), y.vali, pred.y.vali ))
      }
    
    names(cali) <-c("id", "Type","obs.y","pred.y")
    names(vali) <- c("id", "Type","obs.y","pred.y")  
   
    cali.vali <- rbind(cali,vali)
    
    write.table(cali.vali, file=paste(outpath,"/",mod.tech,"_cali.vali_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    ## CROSS - Validate model
    #########################
    if (CR.VAL=="YES"){
    CV.data  <- CV.glm.(gam.tmp, K=fold.cv, cv.lim = 10, name.sp = var.y[tu])  
    
    ## OBS = fct(PRED) - CV set -- follow Pineiro 2008
    ##------------------------
    CV.perf <- asses.lm.(CV.data$predictions,CV.data$obs)
    
    # re-rescale the prediction
    if (RES=="YES"){
      m.temp <- RES.data[,ncol(RES.data)][1]
      sd.temp <- RES.data[,ncol(RES.data)][2]
      
      CV.data$predictions <- (CV.data$predictions*sd.temp) + m.temp
      CV.data$obs <- (CV.data$obs *sd.temp ) + m.temp
    }else{}
    
    # re-transform CV data..
    if (TRANSF=="YES"){
      f <- as.function(inv.transf[[var.list$var.transf[var.list$var.name==var.y[tu]]]])
      CV.data$predictions <- f(x=CV.data$predictions)
      CV.data$obs <- f(x=CV.data$obs)
    }else{}
    
    write.table(CV.data, file=paste(outpath,"/",mod.tech,"_CV_data_rep_",i,".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
    names(CV.perf) <- paste(rep("CV_",5),names(CV.perf),sep="")
    
    if (i ==1) {
      png(paste(outpath,"/",mod.tech,"_OBSvsPRED_rep_",i,".png",sep=""),width = 30, 
          height = 15,units="cm",res=150)
      
        ## Plot observed VS predicted for CV
        plot (CV.data$predictions, CV.data$obs, xlab="Predict_CV", ylab="Observed_CV")
      
      dev.off()
    }else{}
    
    }else{}
    
    ## AIC and AICC of the Best model
    ################################
    if (SVC=="YES"){ 
        RSS <- sum((y.in - (as.numeric(scv.cali[[i]]) *pred.y ))^2) 
      } else{
        RSS <- sum((y.in - pred.y)^2) 
      }
    
    aic.temp <- AIC(gam.tmp) # AIC
    aicc.temp <- AICc(gam.tmp) #AICc
    
    if (CR.VAL =="YES"){
    modl.perf.gam <- rbind(modl.perf.gam, 
                           c(nb.perf=nrow(x.in), M.perf,
                             nb.vali=nrow(x.vali),V.perf,CV.perf, 
                             RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
    }else{
      modl.perf.gam <- rbind(modl.perf.gam, 
                             c(nb.perf=nrow(x.in), M.perf,
                               nb.vali=nrow(x.vali),V.perf, 
                               RSS=RSS, AIC=aic.temp, AICc=aicc.temp))
      
    }
    
  } # end of loop repetition

# write all model perf run
write.table(modl.perf.gam, file=paste(outpath,"/",mod.tech,"_model_perf_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)

#write table for all run
write.table(varimp.gam, file=paste(outpath,"/",mod.tech,"_var_imp_x-run.txt",sep=""),sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)

# calculate avearge values for model perf and variable importance
modl.perf.av <- rbind( modl.perf.av, apply(modl.perf.gam,2,mean) )
varimp.av <- rbind( varimp.av, apply(varimp.gam,2,mean) )

# name technique update
mod.name.perf  <- rbind(mod.name.perf, c(species= var.y[tu],var.sel= var.sele[zu],technic=mod.tech) )
mod.name <- rbind(mod.name, c(species= var.y[tu],var.sel= var.sele[zu],technic=mod.tech) )

} # end loop modeling technique  

##################################################################################################################################### 
## Write average values
# define folder
outpath <- paste(outpath.0,"/",var.y[tu],"/",var.sele[zu],sep="")

# calculate avearge values for variable importance
varimp.av <- data.frame(cbind(mod.name,varimp.av)) 

write.table(varimp.av, file=paste(outpath,"/var_imp_mean.txt",sep=""),sep="\t", 
            append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
} # end loop too much zero in data (RA.O parameter)
} # end loop diff x variable list
  } # end loop does the species exist
    
} # end of loop species
  
#####################################################################################################################################
##################################################################################################################################### 
 ## Write average values
# define folder
outpath <- paste(workdir,"/output/",folder,sep="")

# calculate avearge values for model perf 
modl.perf.av <- data.frame(cbind(mod.name.perf, modl.perf.av)) 
 
#write table for average perf
write.table(modl.perf.av, file=paste(outpath.0 ,"/model_perf_mean.txt",sep=""),sep="\t", 
            append=FALSE, row.names=FALSE,col.names=TRUE, quote=FALSE)

#write table for average perf of SVC
write.table(sv.perf.av, file=paste(outpath.0 ,"/modelSVC_perf_mean.txt",sep=""),sep="\t", 
            append=FALSE, row.names=FALSE,col.names=TRUE, quote=FALSE)

# endCluster() # END OF MULTICORE CALCULATION

# proc.time() - ptm # check time
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#### END SCRIPT ---- COFFE TIMES ----
#################################################################################################################################
#################################################################################################################################




