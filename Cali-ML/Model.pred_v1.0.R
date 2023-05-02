####################################################################################################################################
####################################################################################################################################
###					                     # ##
###  MODEL PREDICTION ONLY		    ###   
###                              ###					
####################################################################################################################################
####################################################################################################################################
## Historic
## --------
## v1.0 DECEMBER 2022 - Boris DROZ -UCC

# PRODOM project - School of Biological, Earth and Environmental Sciences
# Environmental Research Institute (ERI)
# University College Cork
####################################################################################################################################
## DESCRIPTION
###############
## RUN under R-4.2.0 
##
## Performed several prediction machine technics defined by previously 
## calibrate model using script 3.Model.proj_v14.4.R
## 
# Input: - data.table : input predictive variable table
#=======  - folder with previously calibrate model
#
# Output: - model average prediction of several DBPs species
#=======  
# 
####################################################################################################################################
## Load library
################
## Set library path
#.libPaths("C:/Users/Public/Documents/R/win-library/3.1")
# Cheak and download or reload package
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


check.lib(lib.list)
####################################################################################################################################
## SCRIPT PARAMETER  --> NEED TO GO TROUH AND MAKE THE APPROPRIATE MODIFICATION !!!
#####################
# Set folder 
## =========
workdir <- "C:/Users/bodro/Documents/Actual/EPA 2019-W-MS43 - PRODOM/Data_analysis/ML_pred"

model.folder <- "/output/2022-11-09_mod.proj_100run"

f.dat <- "input/input_data_MLv3.csv"# File name of data points
f.var <-"input/var.list.txt"# File with the variable list

## choose the name of the variable selection you will use 
## -- var.list column name with 
## 0 for not use and 1 for use
var.sele <- c("ML_opt3") # one option at the time

## -- option her should be similar as calibrated model
#######################################################
TRANSF <-"YES" # transform the data according to the var.list file ## YES or NO
RES <- "YES" # rescalling to center the data ## YES or NO
n.rep <- 100 ## number of replicate
nb.tree <- 100 # number of tree
####################################################################################################################################
####################################################################################################################################
## ADD IN FUNCTION
#########################################################
#######################################################################################################
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

####################################################################################################################################
####################################################################################################################################
## SCRIPT START HERE
####################################################################################################################################

################################################################################################################################
# -- LOADING DATA SET -----
##############################################
# Open the data set
setwd(workdir)
data.0 <- read.table(f.dat,header = TRUE ,
                      na.strings = "NA", sep=",")

var.list <- read.table(f.var, header= TRUE)

###############################################################################
###############################################################################
## ML - DBP function....

## set model path
inpath <- paste(workdir,model.folder,sep="")
setwd(inpath) # set environment

# interogate how mainy folder
list.species <- list.dirs(path = inpath, full.names = FALSE, recursive = FALSE)

## load subset of data.in
##################
data.in <- data.0

pos <-match( var.list$var.name[var.list[,names(var.list)==var.sele]==1 & 
                                 var.list$var.type=="x"], names(data.in))
data.in <- data.in[,pos]

out.plot <- NULL

for (q in 1:length(list.species)) # species loop
{
  cat("> MODEL SPECIES ",list.species[q],"is running  ...", "\n",append = FALSE)
  cat("#######################################", "\n",append = FALSE)
  
  list.var <- list.dirs(path = paste(inpath,"/",list.species[q],sep=""),
                        full.names = FALSE, recursive = FALSE)
  j<- match(var.sele, list.var)  # var. selection
  
  list.tech <- list.dirs(path = paste(inpath,"/",list.species[q],"/",list.var[j],sep=""),
                         full.names = FALSE, recursive = FALSE)
  
  list.tech <- list.tech[!list.tech=="SVC"] # only quantitative model
  list.tech <-list.tech[!list.tech=="sens.analyse"] # not sens.analysis
  
  # check if one model SVC was build
  output <- paste(inpath,"/",list.species[q],"/",list.var[j],sep="")
  
  model.list <- paste( rep(paste(output,"/","SVC/SVC_model_",sep=''),n.rep),
                       seq(from=1, to=n.rep, by=1), sep="")
  
  if ( sum( file.exists(model.list) )>=1)
  {
    ## svm
    svm.f <- NULL
    pos <- 1
    for (m in 1:n.rep)  
    {
      if ( file.exists(paste(output,"/","SVC/SVC_model_",m,sep='')) ==TRUE)
      {
        mod.tech <- "SVC"
        load(paste(output,"/",mod.tech,"/",mod.tech,"_model_",m,sep=''))
        
        pred <- predict(fit.svm, newdata=data.in, 
                        decision.values = FALSE,
                        probability = FALSE, na.action = na.omit)
        
        pred <- as.numeric(as.character(pred))
        
        if (pos == 1) 
        {svm.f <- pred
        }else{ svm.f <- cbind(svm.f,pred)}
        
        pos <- pos+1
        
      }else{ }
    }# end loop svm rep 
    svm.f <- apply(svm.f,1,mean, na.rm=TRUE)
    svm.f <- round(svm.f,0)
  }# end loop svm
  
  if (identical(list.tech, character(0)) ) # check if quantitative model exist
    {
    ##qualitative model applied - 0/1
    t.model <- rep("quali",length(svm.f))
    
    if (q == 1) 
        {out.plot <- data.frame(type= t.model,mean= svm.f) 
        names(out.plot)[1] <- list.species[q]
    }else{ df <- data.frame(type= t.model,mean= svm.f)
        names(df)[1] <- list.species[q]
        out.plot <- cbind(out.plot,df) }
    
    }else{ 
       ## quantitative model applied
      ## applied modified function 
      for (m in 1:ncol(data.in) )
      {
        f <- as.function(transf[[var.list$var.transf
                                 [var.list$var.name==names(data.in[m]) ]  ]]) # define function from the list
        data.in[,m] <- f(x=data.in[,m])# applied modification
      }
      
      # applied rescaling function (Z-score)  
      par.resc <- read.table(paste(output,"/RESCALING_values.txt", sep=""),header=TRUE)
      
      for (m in 1:ncol(data.in) )
          {
            ## RESCALE THE VARIABLE with normal rescaling
            #f.res <- function(x) {(x-mean(x, na.rm=TRUE))/sd(x,na.rm=TRUE)}
            data.in[,m] <- (data.in[,m] - par.resc[1,m])/par.resc[2,m]
          }
      
      ### check if Na or inf data --> just remove the line
      data.in <- na.omit(data.in)
      data.in <- data.in[!is.infinite(rowSums(data.in)),]
      
      ##########################################################################################################
      ## (3) RUN ML MODEL
      ##########################################################################################################
      t.model <- NULL
      ###############################################################################################################
      if (any("NNET" == list.tech)) {
        
        for (i in 1:n.rep)    
        {
          cat("> NNET MODEL run ",i,"predict  ...", "\n",append = FALSE)
          
          mod.tech <- "NNET"
          
          load(paste(output,"/",mod.tech,"/",mod.tech,"_model_",i,sep=''))
          
          pred <- predict(net.tr, newdata=data.in, type='raw', na.rm = TRUE)
          
          # unscale the data
          #if (res.q== "NO") {
          #}else{
          pred.mean <- par.resc[1,11]
          pred.sd <- par.resc[2,11]
          
          pred <- (pred * pred.sd) + pred.mean
          #}
          
          # back transformed
          f <- as.function(inv.transf[[var.list$var.transf[var.list$var.name==list.species[q]]]])
          pred <- f(pred)
          
          if (i == 1) 
          {t.model <- pred
          }else{ t.model <- cbind(t.model,pred) }
          
        }# end loop model replicate
        
      } # end of if model 
      ###############################################################################################################
      if (any("BAG" == list.tech)) {
        
        for (i in 1:n.rep)    
        {
          cat("> BAG MODEL run ",i,"predict  ...", "\n",append = FALSE)
          
          mod.tech <- "BAG"
          
          load(paste(output,"/",mod.tech,"/",mod.tech,"_model_",i,sep=''))
          
          pred <- predict(bag.tr, newdata=data.frame(data.in), type="raw",na.rm = TRUE ) 
          
          #if (res.q== "NO") {
          #   }else{
          # unscale the data
          pred.mean <- par.resc[1,11]
          pred.sd <- par.resc[2,11]
          
          pred <- (pred * pred.sd) + pred.mean 
          #}
          
          # back transformed
          f <- as.function(inv.transf[[var.list$var.transf[var.list$var.name==list.species[q]]]])
          pred <- f(pred)
          
          #if (i == 1) 
          # {t.model <- pred
          # }else{ 
          t.model <- cbind(t.model,pred)#}
          
        }# end loop model replicate
        
      } # end of if model 
      ############################################################################################################### 
      if (any("GBM" == list.tech)) {
        
        for (i in 1:n.rep)    
        {
          cat("> GBM MODEL run ",i,"predict  ...", "\n",append = FALSE)
          
          mod.tech <- "GBM"
          
          load(paste(output,"/",mod.tech,"/",mod.tech,"_model_",i,sep=''))
          
          pred <- predict( gbm.tr , newdata=as.data.frame(data.in),
                           n.trees=nb.tree, type="response",na.rm = TRUE) 
          
          # if (res.q== "NO") {
          #    }else{
          # unscale the data
          pred.mean <- par.resc[1,11]
          pred.sd <- par.resc[2,11]
          
          pred <- (pred * pred.sd) + pred.mean 
          #  }
          
          # back transformed
          f <- as.function(inv.transf[[var.list$var.transf[var.list$var.name==list.species[q]]]])
          pred <- f(pred)
          
          # if (i == 1) 
          # {t.model <- pred
          #}else{ 
          t.model <- cbind(t.model,pred) #}
          
        }# end loop model replicate
        
      } # end of if model 
      ############################################################################################################### 
      ## creat average and sd
      #########################
       ## - check if svm
      if ( sum( file.exists(model.list) )>=1)
      {
        t.model <-  t.model *svm.f
      } else{}
      
      m.model <-apply(t.model,1,mean)
      sd.model <- apply(t.model,1,sd)
      t.model <- rep("quanti",length(t.model))
      
      if (q == 1) 
      {out.plot <- data.frame(type= t.model,mean= m.model,sd= sd.model) 
      names(out.plot)[1] <- list.species[q]
      }else{ df <- data.frame(type= t.model,mean= m.model,sd= sd.model)
              names(df)[1] <- list.species[q]
        out.plot <- cbind(out.plot,df) }
}# end if quantitative model exist
  
} # end loop for diff species

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#### END SCRIPT ---- COFFE TIMES ----
#################################################################################################################################
#################################################################################################################################




