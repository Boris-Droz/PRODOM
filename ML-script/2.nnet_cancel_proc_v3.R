####################################################################################################################################
####################################################################################################################################
###					                ###
###  Cancelation procedure 
###   to help defining opt. pred variable	    ###   
###                         ###					
####################################################################################################################################
####################################################################################################################################
## Historic
## --------
## v2.0 - June 2015 -  Boris Droz &Gerrad Jones ETHZ & EAWAG --> function develloped
## v3.0 - October 2022 - B. Droz --> used for PROMOD proj.

# PRODOM project - School of Biological, Earth and Environmental Sciences
# Environmental Research Institute (ERI)
# University College Cork
####################################################################################################################################
## DESCRIPTION
###############
##
## Performed euronal Network input cancellation FUNCTION 
##  as a procedure recommaneded by H?ctor F. Satiz?bal M. (2007)

##        Machine learning model
##        ----------------------
##                  - Neural Networks (nnet -- nnet package v.7.3-17)
##
# Input: - data.table : with cali and vali dataset
#=======  
#
# Output: - result of the cancelation procedure
#=======  
# 
####################################################################################################################################
## Load library
##################

## Set library path
# work under R-4.2.0 [64bits]
#.libPaths("C:/Users/Public/Documents/R/win-library/4.2")
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
lib.list <- c("nnet","MASS","qpcR", "e1071")

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
var.sele <- c("Test_Cancel") 

## SAMPLING of each proportion of data at each replic
#====================================================
sampl <- 80 ## percent kept for analysis --usual 90

# transform the data according to the var.list file
TRANSF <-"YES" ## YES or NO

# rescalling to center the data
RES <- "YES" ## YES or NO

# suport vector classification
SVC <- "NO"
k.nel <- "radial" # kernel type

# ratio occurrence -- 
## if proportion in cali data set bigger stop it
RA.O <- 0.7

## -- CHOOSE MODELLING TECHN:
## YES or NO
#############################
# machine learning method
NNET <- "YES"

## NNET PARAMETER
##################
## Weight decay Folow recomendation of B. D. Ripley: "Pattern Recognition and Neural Networks", Cambridge, 1996.
## between 0.1 and 0.01
deca <- 0.1
UNIT.Max <- 15 # NUMBER OF HIDDEN UNITS --> between the number of input nodes and number of output nodes 
it.max <- 10000 # Need to be 10'000 to be MONTE-CARLO Permutation

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

#########################################################
########################################################
## Neuronal Network input cancellation FUNCTION
## ============================================
##
## v2.0 - June 2015 -  Boris Droz &Gerrad Jones ETHZ & EAWAG
## v3.0 - October 2022
## 
## Procedure recommaneded by H?ctor F. Satiz?bal M. (2007)
###########################################################
## Parameter
## =========
## n : number of hidden nods of the nnet
## w : number of weight decay

nnet.cancel <- function (x.train,y.train,n,w,it.max= 10000)
{
  
  require(nnet);
  require(MASS);
  require(qpcR);
  
  # INITIATE A NULL TABLE
  rmse.table <- NULL;
  
  x.tr <- x.train # temporary x data set for training nnet
  
  for (n in 1:(ncol(x.train)-1))
  {
    
    cat(paste("var.",names(x.train)[n],"is running..."),append=TRUE, sep="\n")
    
    # TRAIN NEURAL NETS
    data.in <- cbind(y.train,x.tr)
    
    net <- nnet(x.tr,y.train, data=data.in, size = n, 
                linout = TRUE, maxit = it.max, decay = w,
                trace = FALSE, MaxNWts = 5000)
    
    # PREDICT THE OBS on the train data set 
    pred.y <- predict(net, newdata=x.tr,type="raw",na.rm = TRUE)
    
    # PERFORMANCE OF THE MODEL
    lin.corr <- lm(y.train ~ pred.y) #Matrix ."slope","Int","R2","RMSE"
    
    if (nrow(summary(lin.corr)$coefficients)==1) # Fuck up check
    { 
      # TRAIN NEURAL NETS
      net <- nnet(x.tr,y.train, data=data.in, size = n, 
                  linout = TRUE, maxit = it.max, decay = w,
                  trace = FALSE)
      
      # PREDICT THE OBS
      pred.y <- predict(net, newdata=x.vali,type="raw",na.rm = TRUE)
      
      # PERFORMANCE OF THE MODEL
      lin.corr <- lm(y.train ~ pred.y)  
      
    }else{ }
    
    slope <- summary(lin.corr)$coefficients[2,1] #slope
    std.slope <- summary(lin.corr)$coefficients[2,2]
    int. <- summary(lin.corr)$coefficients[1,1] #Intercept
    std.int <- summary(lin.corr)$coefficients[1,2]
    rsquare <- summary(lin.corr)$r.squared #r squared
    RMSE <- sqrt(mean((y.train - pred.y)^2,na.rm = TRUE)) #RMSE
    
    perf.model <- c(slope=slope, std.slope=std.slope, int.=int., std.int=std.int, rsquare=rsquare, RMSE=RMSE)
    
    train.rmse <- perf.model[6]
    
    # NAME OF VARIABLE USED
    var.sel <- paste((dimnames(x.tr)[[2]]),"/",sep="",collapse="")
    
    # KEEP RESULT TO SEE EVOLUTION OF THE PERFOMANCE
    if (n == 1) { rmse.table <- c(Var.sel = var.sel, perf.model) }else{  
      rmse.table <- rbind( rmse.table, c(Var.sel = var.sel, perf.model))}
    
    diff.table <- NULL
    
    for (m in 1:ncol(x.tr))
    {
      # REPLACE ONE VARIABLE BY THAN AVERAGE
      x.temp <- x.tr
      x.temp [,m] <- mean(x.tr[,m])
      
      # CALCULATE RMSE FOR NEW SET 
      new.rmse <- sqrt(mean( (y.train - predict(net, newdata=x.temp,type="raw" ))^2 ))
      
      diff <- (train.rmse - new.rmse)
      
      # APPEND EACH diff TO A VECTOR
      if (m == 1)  diff.table <- diff else  diff.table <- rbind( diff.table, diff);
      
    }
    rownames(diff.table) <- colnames(x.tr) # RE-NAME NEW TABLE
    
    ## FIND THE LESS RELEVANT INPUT
    l.r <- order(diff.table[,1],decreasing = TRUE)
    
    ## DISCARD THE LESS RELEVANT VARIABLE
    x.tr <- x.tr[,-l.r[1]]
  }
  
  rownames(rmse.table) <- rev(seq (from=1, to= (ncol(x.train)-1), by=1))   
  rmse.table <- as.data.frame(rmse.table)
  
  rmse.table <- rmse.table[order(rmse.table$RMSE,decreasing = FALSE),]
  
  return(rmse.table)
}

####################################################################################################################################
####################################################################################################################################
## SCRIPT START HERE
####################################################################################################################################
# Set folder --> 
date <- Sys.Date()
inpath <- paste(workdir,"/input/", sep="")
folder <- paste("/",date,"_cancel.proc", sep="")
outpath <- creat.subDir(paste(workdir,"/output",sep=""), folder)

## PARRALLEL CORE SETUP# 
##--------------------
# ptm <- proc.time()# ignite timer

# beginCluster(detectCores()) # ACTIVATE THIS MULTI CORE CALCULATION 
###############################################################################
### write the option of the run
################################
f.info <- paste(outpath,"/AA_INFO_PARA.txt",sep="")

cat( paste("*** CANCALATION PROCEDURE RUNNING OUT ---", Sys.Date()), file= f.info, sep="\n")
cat("###########################",file= f.info,append=TRUE, sep="\n")
cat("R-script nnet_cancel_proc_v1 - Droz 2022",file= f.info,append=TRUE, sep="\n")
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

cat("###########################",file= f.info,append=TRUE, sep="\n")
cat( "FUNCTION PARAMETER" ,file= f.info,append=TRUE, sep="\n")
cat("###########################",file= f.info,append=TRUE, sep="\n")            
cat(paste("NNET para - weight decay: ", deca),file= f.info,append=TRUE, sep="\n")
cat(paste("NNET para - number max of hidden units: ", UNIT.Max),file= f.info,append=TRUE, sep="\n")
cat(paste("NNET para - iteration max: ", it.max),file= f.info,append=TRUE, sep="\n")

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

# temporary value to compare data at the end
# --- AVERAGE TABLE ---
mod.name.perf <- NULL
modl.perf.av <- NULL
sv.perf.av <- NULL

for (tu in 1:length(var.y))# selected variable y to pred
{
  cat(paste("SELECTION FOR",var.y[tu],"is running..."),append=TRUE, sep="\n")
  cat("########################################",append=TRUE, sep="\n")
  
  for (zu in 1:length(var.sele)) ## list of selected variable
  {
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
  i=1 # list of sampling data
    n.s <- sample(1:nrow(data.in),nrow(data.in)*sampl/100, replace=FALSE)
    d.cali[[i]] <- data.in[n.s,]
    d.vali[[i]] <- data.in[-n.s,]
    d.cali[[i]] <- na.omit(d.cali[[i]])
    d.vali[[i]] <- na.omit(d.vali[[i]])
    d.cali[[i]] <- d.cali[[i]][is.finite(rowSums(d.cali[[i]]) ),]
    d.vali[[i]] <- d.vali[[i]][is.finite(rowSums(d.vali[[i]]) ),]
  
  ### APPLIED SVC - Support vector classification
  ###############################################
  if (SVC=="YES"){
    
    scv.cali <- list(NULL)
    scv.vali <- list(NULL)
    acc.svc <- NULL
    
      y <- d.cali[[i]][,ncol(d.cali[[i]])] 
      y[y>0] <- 1 # classifier
      
      y.vali <-d.vali[[i]][,ncol(d.vali[[i]])] # data real
      y.vali[y.vali>0] <- 1 # classifier
      
      if (length(y)==sum(y)|length(y.vali)==sum(y.vali)) # if no zero
      {
        scv.cali[[i]] <- y
        scv.vali[[i]] <- y.vali
        
        acc.svc <-rbind(acc.svc,c(n.run=i,acc.cali="no zero", acc.vali="no zero"))
        
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
      
        acc.svc <- rbind(acc.svc,c(n.run=i,acc.cali=acc.cali, acc.vali=acc.vali))
      }
    
      write.table(acc.svc, file=paste(outpath,"/",mod.tech,"_",var.y[tu] ,"_accuracy_x_run.txt",sep=""),
                  sep="\t", append=FALSE, row.names=FALSE,col.names=TRUE, quote=FALSE)
    
    }else{}
  
  #########################
  ### APPLIED DATA MODIFICATION
  ##########################
  if (TRANSF=="YES"){
    for (n in 1:length(d.cali))
        {
        for (i in 1:(ncol(d.cali[[n]])-1) )
            {
              f <- as.function(transf[[var.list$var.transf[i]]]) # define function from the list
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
    d.mean <- apply(d.all,2,FUN=mean,na.rm=TRUE)
    d.sd <- apply(d.all,2,FUN=sd, na.rm=TRUE)
    
    RES.data <- rbind(mean.var=d.mean,sd.var=d.sd)
    
    # write RESCALE table
    write.table(RES.data, file=paste(output,"/RESCALING_values.txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
    
   
        ## RESCALE THE VARIABLE with normal rescaling
        #f.res <- function(x) {(x-mean(x, na.rm=TRUE))/sd(x,na.rm=TRUE)}
        d.cali[[1]] <- (d.cali[[1]]-d.mean)/d.sd
        d.vali[[1]] <- (d.vali[[1]]-d.mean)/d.sd
        }
   else{}
  
    # remove NA/NAN/Inf and produced clean data set of cali vali
      d.cali[[1]] <-na.omit(d.cali[[1]]) # remove na
      d.vali[[1]] <-na.omit(d.vali[[1]])
      
      d.cali[[1]] <-d.cali[[1]][is.finite(rowSums(d.cali[[1]])),] # remove inf
      d.vali[[1]] <-d.vali[[1]][is.finite(rowSums(d.vali[[1]])),] 
      if (SVC=="YES"){
        scv.cali[[1]] <- scv.cali[[1]][apply(is.na(d.cali[[1]]),1,sum)==0 ] # remove na
        scv.vali[[1]] <- scv.vali[[1]][apply(is.na(d.vali[[1]]),1,sum)==0 ]
        
        scv.cali[[1]] <- scv.cali[[1]][is.finite(rowSums(d.cali[[1]]))]# remove inf
        scv.vali[[1]] <- scv.vali[[1]][is.finite(rowSums(d.vali[[1]]))]
      }else{}
      
#################################################################################################################################
#################################################################################################################################
## initialise the list data
# --- model performance ---
modl.perf.nnet <- NULL

count.while <- 0
pred.st.list <-list(NULL)

#####################################################################################################################################
#####################################################################################################################################      
## -------------------------   MODEL CALIBRATION ------------------------------------------
#####################################################################################################################################   
  
i=1 ### just one run #####
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
  
}else{ 
  y.in <- data.frame(d.cali[[i]][,names(d.cali[[i]])==var.y[tu]])
  names(y.in) <- var.y[tu]
  x.in <- d.cali[[i]][, na.omit( match(names(d.cali[[i]]), 
                                       var.list$var.name[var.list$var.type=="x"] ) ) ]
  data.in <- cbind(y.in,x.in)
  
}
  y.in <- y.in[[1]]
  
  if (length(y.in)<=4) {} else{ # need of minimal 4 data to performed lin.reg.
    ######################
    ### NETWORK ANALYSIS
    #####################                       
    ## OPTIMIZED PARAMETER PERFORMANCE OF NNET
    ##########################################
      ##(from the FAQ for a commercial neural network software company)  
      ## Number of inputs + outputs) * (2/3) -- Heaton 2008
      nb.hid1 <- round (ncol(data.in) *2/3)
      nb.hid2 <- round ( nb.hid1*2/3 )
      nb.hidden =c(nb.hid1,nb.hid2)
      
      ## cancelation fucntion
      t.out <- nnet.cancel (x.in,y.in,n=nb.hidden,w=deca,it.max= 10000)
    
    write.table(t.out, file=paste(outpath,"/cancel.proc_",var.y[tu],".txt",sep=""),
                sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)
  }   
 }# end loop diff x variable list
  
} # end of loop species

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#### END SCRIPT ---- COFFE TIMES ----
#################################################################################################################################
#################################################################################################################################




