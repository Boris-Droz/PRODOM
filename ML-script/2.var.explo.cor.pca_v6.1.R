####################################################################################################################################
####################################################################################################################################
###					
###  EXPLONATORY ANALYSIS
### ===========================
###-->  VARIABLE SELECTION
###					
####################################################################################################################################
####################################################################################################################################
## HISTORIC:
## --------
##          v1.0, Decembre 2015 - Boris DROZ & Gerrad Jones, ETHZ & EAWAG
##                Projet Se around the world
##          v2-4 modified December 2016 - B.Droz copper project!
##		      v5.5 July 2017 - - choose analysis and divers para.
##          v6 July 2022 - adapt for non spatial data

# PRODOM project - School of Biological, Earth and Environmental Sciences
# Environmental Research Institute (ERI)
# University College Cork
####################################################################################################################################
## DESCRIPTION
###############
# 
# Transform all continuous var. --> Normalized !!!
#
# Performed divers analysis in order to select the best var. for further model
#           iii)   cor plot
#           iv)  cluster analysis
#           v)   PCA 

# Input: 
# =====
# - data set
# - var.list --> defined all variable ( header: var.type (x or y)	var.name	var.transf)
# - remove list --> defined during the PARAFAC analysis -- just copie it on the appropriate folder

# Output: 
#=======
# Boxplot of all data between resample and remain
# Correlation predict between all variable
# PCA analyis for the data set
# NNET:
#   - Model performance: i) plot OBS vs PRED --> slope, Intercept, R2 
#                            and Root of mean square error between the outputs and targets (RMSE)
#                         ii) MOdel --> AIC, AICc and ratio between nb iteration and nb iteration max (ratio_CONV)
#
####################################################################################################################################
## Load library
##################

## Set library path
#.libPaths("C:/Users/Public/Documents/R/win-library/3.1")

## Cheak and download or releaod package
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

lib.list <- c("ade4","boot","MASS","pracma","qpcR","pvclust","nFactors","ape",
              "foreign","ggplot2","gplots","plyr","Hmisc","corrplot","stringr",
              "rms","gam","nnet","neuralnet","NeuralNetTools","glmnet",
              "randomForest", "AICcmodavg","snowfall", "parallel")

check.lib(lib.list)

####################################################################################################################################
## SCRIPT PARAMETER
#####################
####################################################################################################################################
## SCRIPT PARAMETER  
#####################
## FIRST creat an input and output folder on the workdir
workdir <- "C:/Users/bodro/Documents/Actual/EPA 2019-W-MS43 - PRODOM/Data_analysis/ML_pred"

f.dat <- "input_data_MLv3.csv"# File name of data points
f.var <-"var.list.txt"# File with the variable list
f.rem <- "exclude.list.txt" # File with the list of the sample to remove

## choose the name of the variable selection you will use -- var.list column name with 
## 0 for not use and 1 for use
var.sele <- "Test_Cancel"

## Choose which analysis to performed
#######################################
CORR <- "YES"
CLUST <- "NO"
PCA <- "YES"

## NUMBER OF replication
# ========================
n.rep <- 5
## SAMPLING of each proportion of data at each replic
#====================================================
resampling <- 100 ## percent kept for analysis

# transform the data according to the var.list file
TRANSF <-"NO" ## YES or NO

# rescalling to center the data
RES <- "NO" ## YES or NO

# cluster group for pca analysis
n.g <- 3
####################################################################################################################################
####################################################################################################################################
####################################################################################################################################
##########################
##########################
##########################
##                        ##
##    ADD-IN FUNCTIONS    ##
##                        ##
##########################
########################################################################################################################################
#########################################################
########################################################
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

#################################################################################################################################################
# FUNCTION check and produced subDir folder
###########################################
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
#################################################################################################################################################

####################################################################################################################################
####################################################################################################################################
## SCRIPT STAR HERE
####################################################################################################################################
# Set folder --> 
date <- Sys.Date()
inpath <- paste(workdir,"/input/", sep="")
folder <- paste("/",date,"_variable.explo_",var.sele, paste="", sep="")
outpath <- creat.subDir(paste(workdir,"/output",sep=""), folder)
#outpath <-paste(workdir,"/output",folder,"/",sep="")
#creat.subDir(paste(workdir,"/output",sep=""), folder)

## set your work space
setwd(inpath)

## PARRALLEL CORE SETUP
##--------------------
# beginCluster(detectCores()) # ACTIVATE THIS MULTI CORE CALCULATION 

##############################################
# -- LOADING AND ORGANIZED DATA SET -----
##############################################
# Open the data set
data.in <- read.table(paste(inpath,f.dat,sep=""),header = TRUE ,
                      na.strings = "NA", sep=",")

var.list <- read.table(paste(inpath,f.var,sep=""), header= TRUE)
rem.list <- read.table(paste(inpath,f.rem,sep=""), header= FALSE)

########################################
## Creat Organized data set for the test
## --> remove selected sampl ( remove list)
## --> select only variable (x and y) who as one in the selected list 
## --> transformed to normalized it
var.list <-na.omit(var.list)
data.in <- data.in[!(data.in$sample %in% rem.list$V1),]
p. <- names(var.list)== var.sele
pos <-match(var.list$var.name[var.list[,p.]==1], names(data.in))
data.in <- data.in[,pos]
var.list <- var.list[var.list[,p.]==1,]

##########################
### APPLIED DATA MODIFICATION 
if (TRANSF=="YES"){
for (i in 1:ncol(data.in))
  {
  f <- as.function(transf[[var.list$var.transf[i]]]) # define function from the list
  data.in[,i] <- f(x=data.in[,i]) # applied modification
  #print(paste(names(data.in[i]), length(data.in[,i][is.finite(data.in[,i]) ]) ) ) 
  }
    }else{}

if (RES=="YES"){
## RESCALE THE VARIABLE with normal rescaling
f.res <- function(x) {(x-mean(x, na.rm=TRUE))/sd(x,na.rm=TRUE)}
data.in <- apply(data.in,2,f.res)} else{}
  
data.in <- na.omit(data.in) # delet row with NA value

x.clu <- apply(data.in, 2,function(x) {is.finite(x)})
data.in <-data.in[apply(x.clu, 1,all), ] # remove if infinite value

d.samp <- data.in # keep local tble

### write the option of the run
################################
f.info <- paste(outpath,"/A_INFO_PARA.txt",sep="")
cat( paste("*** RUN VAR SELECTION ---", Sys.Date()), file= f.info, sep="\n")
cat("###########################",file= f.info,append=TRUE, sep="\n")
cat(paste("RUN for the data set:", f.dat),file= f.info,append=TRUE, sep="\n")
cat(paste("Are variable transform?", TRANSF),file= f.info,append=TRUE, sep="\n")
cat(paste("Do you have rescaling?", RES),file= f.info,append=TRUE, sep="\n")
cat(paste("name of variable list do you use?", var.sele),file= f.info,append=TRUE, sep="\n")
cat(paste("remove sample list:", f.rem),file= f.info,append=TRUE, sep="\n")
cat("###########################",file= f.info,append=TRUE, sep="\n")
cat( "VARIABLE SELECTION TEST" ,file= f.info,append=TRUE, sep="\n")
cat("###########################",file= f.info,append=TRUE, sep="\n")
cat(paste("correlation plot?", CORR),file= f.info,append=TRUE, sep="\n")
cat(paste("cluster group?", CLUST),file= f.info,append=TRUE, sep="\n")
cat(paste("if cluster group, how much group", n.g),file= f.info,append=TRUE, sep="\n")
cat(paste("principale component analysis?", PCA),file= f.info,append=TRUE, sep="\n")
cat("###########################",file= f.info,append=TRUE, sep="\n")
cat(paste("number of replication?", n.rep),file= f.info,append=TRUE, sep="\n")
cat(paste("percent of sampling:", resampling),file= f.info,append=TRUE, sep="\n")
cat(paste("Total number of data (cali+test):", nrow(data.in)),file= f.info,append=TRUE, sep="\n")
 
####################################################################
################################################################################################################################
## VARIABLE SELECTION
##=====================
  # Prepare output matrix - variable
  mat.var.pca <- matrix(0,n.rep,(length(var.list)-2),dimnames=list(c(1:n.rep),
                paste("Var.Axe_",seq(from=1,to=(length(var.list)-2),by=1),sep="")))
  
  mat.var.pca.l <- list (mat.var.pca, mat.var.pca ,mat.var.pca, mat.var.pca)
  
  # PCA LIST
  # Prepare output matrix - all 
  mat.all.pca <- matrix(0,n.rep,(length(var.list)-2),dimnames=list(c(1:n.rep),
                paste("All.Axe_",seq(from=1,to=(length(var.list)-2),by=1),sep="")))
  
  mat.all.pca.l <- list (mat.all.pca,mat.all.pca,mat.all.pca,mat.all.pca)
  
  # CORPLOT LIST
  df.cor.pred.var.l <- list(list(),list(),list(),list())
  
  df.cor.pred.all.l <- list(list(),list(),list(),list())
  
  ## --> end ignitia temporary list
  #################################################################################################################################
  ## START THE VARIABLE SELECTION
  ## ============================
  for (i in 1:n.rep)
      {
      data.in <- d.samp
    
      # resample data 
      pos <- sample(seq(from=1,to=nrow(data.in), by=1), 
                    size = round(nrow(data.in)*resampling/100) )
      data.in <- data.in[pos,]
    
      # creat x and y table
      y.in <- data.in[,match( var.list$var.name[var.list$var.type=="y"],
                            names(data.in) ) ]
    
      x.in <- data.in[,match( var.list$var.name[var.list$var.type=="x"],
                            names(data.in) ) ]
    
        ## SAMPLING
        cat("> Replication ",i," is starting now...", "\n",append = FALSE)
    	  cat("#############################################", "\n",append = FALSE)
        
        t <- 1 # used before with several sampling technique   
        
  #################################################################################################################################     
  #################################################################################################################################
  ## ii) CORR PLOT BETWEEN PREDICTOR
  ## ===============================
  
  if (CORR=="NO"){ }else{
  
  df.cor.pred.var.l[[t]][[i]] <-data.frame(round(cor(x.in),3))
  
  # Save correlations in your workspace
  write.table(df.cor.pred.var.l[[t]][[i]],file=paste(outpath,"/Cor_Pred_var_rep_",i,".txt",sep=""),
              sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
  
  df.cor.pred.all.l[[t]][[i]] <-data.frame(round(cor(data.in),3))
  
  # Save correlations in your workspace
  write.table(df.cor.pred.all.l[[t]][[i]],file=paste(outpath,"/Cor_Pred_all_rep_",i,".txt",sep=""),
              sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
  
  } # end technique
  #################################################################################################################################
  ## iii) CLUSTER ANALYSIS of the pred.variable
  ## ======================
  # ref: Suzuki, R. and Shimodaira, H. (2006)
  
  if (CLUST=="NO"){ }else{
  
  cat("-> Cluster analysis...", "\n",append = FALSE)
  
          # Ward Hierarchical Clustering with Bootstrapped p values
          fit <- pvclust(x.in, method.hclust="ward.D2",
                         method.dist="euclidean",nboot=1000)
  
          jpeg(paste(outpath,"/Pred.var_cluster_rep_",i,".jpg",sep=""),
               width = 30, height = 15,units="cm",res=150)
          
              plot(fit) # dendogram with p values
              # add rectangles around groups highly supported by the data
              pvrect(fit, alpha=.95) 
          
          dev.off()
  
  } # end technique
  #################################################################################################################################
  ## iV) PCA
  ## ========
  if (PCA=="NO"){ }else{
  
        ## PCA - variable ##
        ######################
        pca.pred.var <- dudi.pca(x.in,scannf = FALSE, nf= ncol(x.in) )
        
        # Variance explained by all axe 
        mat.var.pca.l[[t]][i,] <- pca.pred.var$eig[1:(length(var.list)-2)]/sum(pca.pred.var$eig)
        
        ## PCA - all ##
        ######################
        pca.pred.all <- dudi.pca(data.in,scannf = FALSE, nf= ncol(data.in) )
        
        # Variance explained by all axe in ratio
        mat.all.pca.l[[t]][i,] <- pca.pred.all$eig[1:(length(var.list)-2)]/sum(pca.pred.all$eig)
        
        # Plot PCA correlation circle #
        ###############################
        #if (i == 1) { # write a plot only if first run
          jpeg(paste(outpath,"/PCA_CorCircle_rep_",i,".jpeg",sep=""),width = 30, height = 15,units="cm",res=150)
            par(mfrow=c(1,2))
          
              s.corcircle(pca.pred.var$co,clabel=.7, cgrid = 2, full = FALSE, sub = "var.indep", 
                          csub = 2.5, possub = "bottomleft", box = TRUE)
              
              s.corcircle(pca.pred.all$co,clabel=.7, cgrid = 2, full = FALSE, sub = "all.var", csub = 2.5, 
                          possub = "bottomleft", box = TRUE)
          
          dev.off()
        
        
         pdf(paste(outpath,"/PCA_CorCircle_rep_",i,".pdf",sep=""),width = 11, height = 5.5)
             par(mfrow=c(1,2))
              
              s.corcircle(pca.pred.var$co,clabel=.7, cgrid = 2, full = FALSE, sub = "var.indep", 
                          csub = 2.5, possub = "bottomleft", box = TRUE)
              
              s.corcircle(pca.pred.all$co,clabel=.7, cgrid = 2, full = FALSE, sub = "all.var", csub = 2.5, 
                          possub = "bottomleft", box = TRUE)
          
          dev.off()
        #}else{}
        
        ## Compute root mean square loading factor --> pca var. importance
        ##########################################
        mydata <- x.in
        
        # Determine Number of Factors to Extract based on Kaiser (1960) equation
        #library(nFactors)
        ev <- eigen(cor(mydata)) # get eigenvalues
        ap <- parallel(subject=nrow(mydata),var=ncol(mydata),
                       rep=100,cent=.05)
        nf <- nScree(x=ev$values, aparallel=ap$eigen$qevpea) #plotnScree(nf) 
        nf <- nf$Components$nkaiser # extract Kaiser method of selecting PCA var number
        
        load.f <- pca.pred.var$co[,1:nf] # loading factor
        
        f <- function(x){(sum(x^2))^0.5}
        
        rms.load <- apply(load.f,1,f) # root mean square of important var.
        
        if (i==1) { mat.rms.load <- rms.load  
              }else{ mat.rms.load <- rbind(mat.rms.load,rms.load) 
                    row.names(mat.rms.load) <- seq(from=1,to=i, by=1) }
        
        # Save correlations in your workspace
        write.table(mat.rms.load,file=paste(outpath,"/pca_var_imp_",i,".txt",sep=""),
                    sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
  
  ## compute euclidien distance (on PCA1 and PCA2)
  ##########################################
  
  fun <- function(x){sqrt(sum(x[1]^2,x[2]^2))} # distance of vector two
  
  euc.dist <- apply(pca.pred.var$co[,1:2],1, fun)
  
  if (i==1) { mat.euc.dist<- euc.dist  
              }else{ mat.euc.dist <- rbind(mat.euc.dist,euc.dist) 
  			row.names(mat.euc.dist) <- seq(from=1,to=i, by=1) }
  
        ## Plot PCA point classes 
        ###########################
  
        ## classified data point per cluster analysis
        ##############################################
        ## Zonal analysis 
        ############################################
        # Ward Hierarchical Clustering without Bootstrapped 
        # simple analysis by 5 classes
        d <- dist(mydata, method = "euclidean") # distance matrix
        fit <- hclust(d, method="ward.D")
  
        groups <- cutree(fit, k= n.g) # cut tree into 5 clusters
        
        jpeg(paste(outpath,"/PCA_clustering_rep_",i,".jpg",sep=""),width = 30, height = 15,units="cm",res=150)
            par(mfrow=c(1,2))  
            
            plot(fit)
            # draw dendogram with red borders around the 5 clusters
            rect.hclust(fit, k= n.g, border="red") 
  
            # data point
            s.class(pca.pred.var$li, cgrid = 2, fac= as.factor(groups),cstar = 0, 
                cellipse = 0, col=c("#a5cf84","#f052b1","#3d82eb","#e86519","#545454"))
  
        dev.off()
  
  } # end technique
  #################################################################################################################################
  
  }# End of loop n rep
  
  ################################################################################################################################
  #######################   
  # Export output means #
  ######################
  creat.subDir(mainDir= outpath, subDir = "average")
  
    #for ( t in 1:length(s.tech))
         # {
      t <- 1
          ##########################
          ## III) PCA and Cor. PLot
          ##########################
          
          var <- cbind(colMeans(mat.var.pca.l[[t]], na.rm = TRUE), apply(mat.var.pca.l[[t]],2,sd))
          all <- cbind(colMeans(mat.all.pca.l[[t]], na.rm = TRUE), apply(mat.all.pca.l[[t]],2,sd))
          
          dimnames(var) <- list(paste("Var._",seq(from=1,to=(length(var.list)-2),by=1)),c("Average","SD"))
    
          dimnames(all) <- list(paste("Var._",seq(from=1,to=(length(var.list)-2),by=1)),c("Average","SD"))
    
          write.table(var, file=paste(outpath,"/average/PCA_var.txt",sep=""),
                      sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
    
        	write.table(all, file=paste(outpath,"/average/PCA_all.txt",sep=""),
                      sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
                
          ## COR.Predict.table
          ##-------------------
          df.cor.pred.var <- df.cor.pred.var.l[[t]]
          df.cor.pred.all <- df.cor.pred.all.l[[t]]
          
          ## extract data and mean for cor. pred
          p <- nrow(df.cor.pred.all[[1]])
          names <- dimnames(df.cor.pred.all[[1]])
          
          std.cor.var <- matrix(0, p, p, dimnames = names)
          std.cor.all <- matrix(0, p, p, dimnames = names)
          
          for (q in 1:p)
              {
                for (r in 1:p)
                    {
                      temp.pred.var <- c()
                      temp.all <- c()
                      for (d in 1:n.rep)
                          {
                            temp.pred.var <- c(temp.pred.var, df.cor.pred.var[[d]][q,r])
                            temp.all <- c(temp.all,df.cor.pred.all[[d]][q,r])
                            
                          }
                          std.cor.var[q,r] <- sd(temp.pred.var)
                          std.cor.all[q,r] <- sd(temp.all)    
                    }  
              }
          
          mean.cor.var <- matrix(0, p, p, dimnames = names)
          mean.cor.all <- matrix(0, p, p, dimnames = names)
          
          for (n in 1:n.rep)
              {
                mean.cor.var <- df.cor.pred.var[[n]] + mean.cor.var
                mean.cor.all <- df.cor.pred.all[[n]] + mean.cor.all 
              }
          
          mean.cor.var <- mean.cor.var/n.rep
          mean.cor.all <- mean.cor.all/n.rep
          
          ## write mean and std cor predict
          write.table(std.cor.var,file=paste(outpath,"/average/cor.var.std.cor.txt",sep=""),
                      sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
          write.table(std.cor.all,file=paste(outpath,"/average/cor.all.std.txt",sep=""),
                      sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
          
          write.table(mean.cor.var,file=paste(outpath,"/average/cor.var.mean.txt",sep=""),
                      sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
          write.table(mean.cor.all,file=paste(outpath,"/average/cor.all.mean.txt",sep=""),
                      sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
          
          ## root mean square load of PCA
          write.table(mat.rms.load,file=paste(outpath,"/PCA_rms_loading_x-run.txt",sep=""),
                      sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
          
          # compute average
          av.rms.load <- apply(mat.rms.load,2,mean)
          
          write.table(av.rms.load,file=paste(outpath,"/average/PCA_mean_rms_loading.txt",sep=""),
                      sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
  
  	#compute euclidienne distance
  	write.table(mat.euc.dist,file=paste(outpath,"/PCA_euc.dist_x-run.txt",sep=""),
                      sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
          
          # compute average
          av.euc.dist <- apply(mat.euc.dist,2,mean)
          
          write.table(av.euc.dist,file=paste(outpath,"/average/PCA_mean_euc.dist.txt",sep=""),
                      sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
  
          
          # std of the avearge
          sd.rms.load <- apply(mat.rms.load,2,sd)
          
          write.table(sd.rms.load,file=paste(outpath,"/average/PCA_std_rms_loading.txt",sep=""),
                      sep="\t",append=FALSE,row.names=TRUE,col.names=TRUE,quote=FALSE)
          
          ## COR.Plot
          #----------
          ## need to define as matrix
          var <- as.matrix(mean.cor.var)
          all <- as.matrix(mean.cor.all)
            
          ## Plot option  
          plot.method<-"color"     # options= "circle", "square", "ellipse", "number", "shade", "color", "pie"
          plot.type<-"lower"         # options= "full", "lower", "upper"
          plot.order<-"original"     # options= "original", "AOE", "FPC", "hclust", "alphabet"
                
          col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                                     "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))  
          
          pdf(paste(outpath,"/average/corr_plot.pdf",sep=""),width = 11, height = 5.5)
              par(mfrow=c(1,2))
              
              corrplot(var,method=plot.method,type=plot.type, order=plot.order,outline=T,tl.col="black",
                       tl.offset=0.5,cl.length=10,sig.level=0.05,insig="pch",pch=1,pch.cex=23,col=col2(20))
              
              corrplot(all,method=plot.method,type=plot.type, order=plot.order,outline=T,tl.col="black",
                       tl.offset=0.5,cl.length=10,sig.level=0.05,insig="pch",pch=1,pch.cex=23,col=col2(20))
          
          dev.off()
        
# endCluster() # END OF MULTICORE CALCULATION

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#### COFFE TIMES
#################################################################################################################################
#################################################################################################################################




