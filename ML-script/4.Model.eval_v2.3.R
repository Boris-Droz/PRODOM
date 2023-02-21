# ENSEMBLE Model evaluation
############################
####################################################################################################################################
## Historic
## --------
# v2.0 - B.Droz November 2022 

# PRODOM project - School of Biological, Earth and Environmental Sciences
# Environmental Research Institute (ERI)
# University College Cork
################################################################################
## DESCRIPTION
###############
## RUN under R-4.2.0 
## parameter
#
workdir <- "C:/Users/bodro/Documents/Actual/EPA 2019-W-MS43 - PRODOM/Data_analysis/ML_pred"

model.folder <- "/output/2022-11-09_mod.proj_100run"

list.var.sele <- "ML_UV" # one option at the time

# Load packages
library(data.table)
library(tidyr)
library(dplyr)
library(wesanderson)

############################################################################
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
  RMSE <- sqrt(mean((x - y)^2,na.rm = TRUE)) #root mean square error
  RMSPE <- sqrt(mean(((x - y)/x)^2,na.rm = TRUE)) #root mean square percent error
  
  lm.vector <- c(slope=slope, std.slope=std.slope, int.=int., std.int=std.int, 
                 rsquare=rsquare, MAE=MAE, MSE=MSE, RMSE=RMSE, RMSPE=RMSPE)
  
  return(lm.vector) 
}  

### Zhang 2018 npj model eval
##############################

bias2 <- function(m.pred.vali,obs.vali){
  
      bias2 <- 1/length(obs.vali)*sum((m.pred.vali-obs.vali)^2) 
      return(bias2)}

vari2 <- function(m.pred.vali,m.pred.cali,obs.vali ){ # cali vali for the same id
  
  vari2 <- 1/length(obs.vali)*sum( 1/length(m.pred.cali) * sum( (m.pred.cali-m.pred.vali)^2) )
  return(vari2)}
##############################################################
##############################################################
##
##  -- SCRIPT START HERE ------
##
###############################################################

inpath <- paste(workdir,model.folder,sep="")

setwd(inpath) # set environment

# interogate how mainy folder
list.species <- list.dirs(path = inpath, full.names = FALSE, recursive = FALSE)

modl.perf.av <- NULL
mod.name.perf <- NULL

for (i in 1:length(list.species)) # species loop
    {
  
  list.var <- list.dirs(path = paste(inpath,"/",list.species[i],sep=""),
                                         full.names = FALSE, recursive = FALSE)
  j<- match(list.var.sele, list.var)  # var. selection
  
 # for (j in 1:length(list.var)) # variable selection loop
  #{
    list.tech <- list.dirs(path = paste(inpath,"/",list.species[i],"/",list.var[j],sep=""),
                          full.names = FALSE, recursive = FALSE)
    
    list.tech <- list.tech[!list.tech=="SVC"] # not quantitative model
    list.tech <-list.tech[!list.tech=="sens.analyse"] # not sens.analysis
    
    if (identical(list.tech, character(0)) ) {}else{
      
      mean.model <- NULL
      for (k in 1:length(list.tech)) # list ML technic loop
      { 
        outpath <- paste(inpath,"/",list.species[i],"/",list.var[j],
                         "/",list.tech[k],sep="")
        
        list.run <- list.files(path = outpath ,
                               pattern=paste(list.tech[k],"_cali.vali_rep_",sep=""))
        mean.tech <- NULL
        
        for ( l in 1:length(list.run) ) # variable selection loop
          {  
            dat <- read.table(paste(inpath,"/",list.species[i],"/",list.var[j],
                                    "/",list.tech[k],"/",list.run[l],sep=""),
                              header=TRUE, sep="")
            
            mean.tech <-  rbind(mean.tech,dat) # merge all data together per tech
          }
        
        mean.model <-rbind(mean.model,mean.tech) ## merge all data together
            
          # make mean and sd of the data per id and CALI/VALI
                avg.tech <-   
          mean.tech %>%   # First step in the next string of statements
          group_by(id,Type) %>%   # Groups the summary file by Plot number
          summarize(           # Coding for how we want our CWMs summarized
            m.obs.y = mean(obs.y),
            sd.obs.y = sd(obs.y),
            m.pred.y = mean(pred.y),
            sd.pred.y = sd(pred.y)   )
      
        # evalated the average model per techn
        fit.cali <- asses.lm.(avg.tech$m.pred.y[avg.tech$Type=="CALI"],
                              avg.tech$m.obs.y[avg.tech$Type=="CALI"] )
        
        fit.vali <- asses.lm.(avg.tech$m.pred.y[avg.tech$Type=="VALI"],
                              avg.tech$m.obs.y[avg.tech$Type=="VALI"] )
        
        # evaluated the prediction interval of the model for 95% confidence interval
        si.cali = qt(0.975,df=length(avg.tech$m.obs.y[avg.tech$Type=="CALI"])-1)*
                                sigma(lm(avg.tech$m.obs.y[avg.tech$Type=="CALI"] ~
                                  avg.tech$m.pred.y[avg.tech$Type=="CALI"]))
        
        si.vali = qt(0.975,df=length(avg.tech$m.obs.y[avg.tech$Type=="VALI"])-1)*
                            sigma(lm(avg.tech$m.obs.y[avg.tech$Type=="VALI"] ~
                             avg.tech$m.pred.y[avg.tech$Type=="VALI"])) 
          
        # precision and accuracy of the model
        prec.cali <- avg.tech$sd.pred.y[avg.tech$Type=="CALI"]/
                      avg.tech$m.pred.y[avg.tech$Type=="CALI"]*100
        m.prec.cali <-mean(prec.cali, na.rm=TRUE)
        sd.prec.cali <-sd(prec.cali, na.rm=TRUE)
        
        prec.vali <- avg.tech$sd.pred.y[avg.tech$Type=="VALI"]/
                    avg.tech$m.pred.y[avg.tech$Type=="VALI"]*100
        m.prec.vali <-mean(prec.vali, na.rm=TRUE)
        sd.prec.vali <-sd(prec.vali, na.rm=TRUE)
          
        resi <- (avg.tech$m.pred.y[avg.tech$Type=="CALI"] - avg.tech$m.obs.y[avg.tech$Type=="CALI"])/
                  avg.tech$m.obs.y[avg.tech$Type=="CALI"]*100
        resi <-resi[is.finite(resi)]
        m.resi.cali <- mean(resi,na.rm=TRUE)
        sd.resi.cali <- sd(resi,na.rm=TRUE)
        
        resi <- (avg.tech$m.pred.y[avg.tech$Type=="VALI"] - avg.tech$m.obs.y[avg.tech$Type=="VALI"])/
                  avg.tech$m.obs.y[avg.tech$Type=="VALI"]*100
        resi <-resi[is.finite(resi)]
        m.resi.vali <- mean(resi,na.rm=TRUE)
        sd.resi.vali <- sd(resi,na.rm=TRUE)
        
        # bias and variance of the model
        bi <- bias2(avg.tech$m.pred.y[avg.tech$Type=="VALI"],
              avg.tech$m.obs.y[avg.tech$Type=="VALI"])
        
        # set up the number of cali-vali dataset
        p.cali <- match(avg.tech$id[avg.tech$Type=="VALI"],
              avg.tech$id[avg.tech$Type=="CALI"] )
        p.cali <-na.omit(p.cali)
        p.vali <- match(avg.tech$id[avg.tech$Type=="CALI"],
                        avg.tech$id[avg.tech$Type=="VALI"] )
        p.vali <- na.omit(p.vali)
        
        var <- vari2(avg.tech$m.pred.y[avg.tech$Type=="VALI"][p.vali],
                     avg.tech$m.pred.y[avg.tech$Type=="CALI"][p.cali],
                     avg.tech$m.obs.y[avg.tech$Type=="VALI"][p.vali] )
        
        # PERFORMANCE AND name technique update
        modl.perf.av <- rbind(modl.perf.av,
                              c(cali_slop= fit.cali[1],cali_std.slop=fit.cali[2] ,
                                cali_int= fit.cali[3], cali_std.int= fit.cali[4],
                                cali_R2= fit.cali[5], cali_RMSE=fit.cali[8],
                                cali_RMSPE= fit.cali[9],
                                cali_2sigma=si.cali, cali_m.prec = m.prec.cali,
                                cali_sd.prec = sd.prec.cali,
                                cali_m.resi= m.resi.cali,cali_sd.resi=sd.resi.cali,
                                ## vali set...
                                vali_slop=fit.vali[1] ,vali_std.slop=fit.vali[2] ,
                                vali_int=fit.vali[3] , vali_std.int= fit.vali[4],
                                vali_R2=fit.vali[5], vali_RMSE=fit.vali[8],
                                vali_RMSPE= fit.vali[9],
                                vali_2sigma=si.vali, vali_m.prec = m.prec.vali,
                                vali_sd.prec = sd.prec.vali,
                                vali_m.resi= m.resi.vali,vali_sd.resi=sd.resi.vali,
                                #model overfitting
                                bias2=bi, variance=var))
        
        mod.name.perf  <- rbind(mod.name.perf, 
                                c(species= list.species[i],var.sel= list.var[j],
                                  technic=list.tech[k]) )
        
        # write the average table for the tech. 
        write.table(avg.tech, file=paste(outpath ,"/",list.tech[k],"_cali.vali_mean.txt",sep=""),sep="\t", 
                    append=FALSE, row.names=FALSE,col.names=TRUE, quote=FALSE)
        
      } # end ist ML technic loop
      
      # make mean and sd of the data per id and CALI/VALI
      avg.model <-   
        mean.model %>%   # First step in the next string of statements
        group_by(id,Type) %>%   # Groups the summary file by Plot number
        summarize(           # Coding for how we want our CWMs summarized
          m.obs.y = mean(obs.y),
          sd.obs.y = sd(obs.y),
          m.pred.y = mean(pred.y),
          sd.pred.y = sd(pred.y)   )
      
      # evalated the average model per techn
      fit.cali <- asses.lm.(avg.model$m.pred.y[avg.model$Type=="CALI"],
                            avg.model$m.obs.y[avg.model$Type=="CALI"] )
      
      fit.vali <- asses.lm.(avg.model$m.pred.y[avg.model$Type=="VALI"],
                            avg.model$m.obs.y[avg.model$Type=="VALI"] )
      
      ##### compare the two model with ancova 
      test <- aov(m.obs.y~ m.pred.y+Type, data =avg.model)
        summary(test)      
      out.file <- paste(inpath,"/",list.species[i],"/",list.var[j],"/ancova_res.txt",sep="")
        
        sink(out.file) ## work only with R-studio!
        print(summary(test))
        sink()  # returns o
      
      # evaluated the prediction interval of the model for 95% confidence interval
      si.cali = qt(0.975,df=length(avg.model$m.obs.y[avg.model$Type=="CALI"])-1)*
        sigma(lm(avg.model$m.obs.y[avg.model$Type=="CALI"] ~
                   avg.model$m.pred.y[avg.model$Type=="CALI"]))
      
      si.vali = qt(0.975,df=length(avg.model$m.obs.y[avg.model$Type=="VALI"])-1)*
        sigma(lm(avg.model$m.obs.y[avg.model$Type=="VALI"] ~
                   avg.model$m.pred.y[avg.model$Type=="VALI"])) 
      
      # precision and accuracy of the model
      prec.cali <- avg.model$sd.pred.y[avg.model$Type=="CALI"]/
        avg.model$m.pred.y[avg.model$Type=="CALI"]*100
      m.prec.cali <-mean(prec.cali, na.rm=TRUE)
      sd.prec.cali <-sd(prec.cali, na.rm=TRUE)
      
      prec.vali <- avg.model$sd.pred.y[avg.model$Type=="VALI"]/
        avg.model$m.pred.y[avg.model$Type=="VALI"]*100
      m.prec.vali <-mean(prec.vali, na.rm=TRUE)
      sd.prec.vali <-sd(prec.vali, na.rm=TRUE)
      
      resi <- (avg.model$m.pred.y[avg.model$Type=="CALI"] - avg.model$m.obs.y[avg.model$Type=="CALI"])/
               avg.model$m.obs.y[avg.model$Type=="CALI"]*100
      resi <-resi[is.finite(resi)]
      m.resi.cali <- mean(resi,na.rm=TRUE)
      sd.resi.cali <- sd(resi,na.rm=TRUE)
      
      resi <- (avg.model$m.pred.y[avg.model$Type=="VALI"] - avg.model$m.obs.y[avg.model$Type=="VALI"])/
              avg.model$m.obs.y[avg.model$Type=="VALI"]*100
      resi <-resi[is.finite(resi)]
      m.resi.vali <- mean(resi,na.rm=TRUE)
      sd.resi.vali <- sd(resi,na.rm=TRUE)
      
      # bias and variance of the model
      bi <- bias2(avg.model$m.pred.y[avg.model$Type=="VALI"],
                  avg.model$m.obs.y[avg.model$Type=="VALI"])
      
      # set up the number of cali-vali dataset
      p.cali <- match(avg.model$id[avg.model$Type=="VALI"],
                      avg.model$id[avg.model$Type=="CALI"] )
      p.cali <-na.omit(p.cali)
      
      p.vali <- match(avg.model$id[avg.model$Type=="CALI"],
                      avg.model$id[avg.model$Type=="VALI"] )
      p.vali <- na.omit(p.vali)
      
      var <- vari2(avg.model$m.pred.y[avg.model$Type=="VALI"][p.vali],
                   avg.model$m.pred.y[avg.model$Type=="CALI"][p.cali],
                   avg.model$m.obs.y[avg.model$Type=="VALI"][p.vali] )
      
      # PERFORMANCE AND name technique update
      modl.perf.av <- rbind(modl.perf.av,
                            c(cali_slop= fit.cali[1],cali_std.slop=fit.cali[2] ,
                              cali_int= fit.cali[3], cali_std.int= fit.cali[4],
                              cali_R2= fit.cali[5], cali_RMSE=fit.cali[8],
                              cali_RMSPE= fit.cali[9],
                              cali_2sigma=si.cali, cali_m.prec = m.prec.cali,
                              cali_sd.prec = sd.prec.cali,
                              cali_m.resi= m.resi.cali,cali_sd.resi=sd.resi.cali,
                              ## vali set...
                              vali_slop=fit.vali[1] ,vali_std.slop=fit.vali[2] ,
                              vali_int=fit.vali[3] , vali_std.int= fit.vali[4],
                              vali_R2=fit.vali[5], vali_RMSE=fit.vali[8],
                              vali_RMSPE= fit.vali[9],
                              vali_2sigma=si.vali, vali_m.prec = m.prec.vali,
                              vali_sd.prec = sd.prec.vali,
                              vali_m.resi= m.resi.vali,vali_sd.resi=sd.resi.vali,
                              #model overfitting
                              bias2=bi, variance=var))
      
      mod.name.perf  <- rbind(mod.name.perf, 
                              c(species= list.species[i],var.sel= list.var[j],
                                technic="average_model") )
      
      outpath <- paste(inpath,"/",list.species[i],"/",list.var[j],sep="")
      
      pdf(paste(outpath,"/",list.species[i],"_OBSvsPRED_average.pdf",sep=""),width = 5.5, height = 5.5)
      
      west.col <- wes_palette("FantasticFox1", n=5, type = c("discrete"))
        ## Plot observed VS predicted 
        plot (NULL, NULL,
              xlab="Predict", ylab="Observed",
              las=1,
              xlim=c(min(avg.model$m.pred.y[avg.model$Type=="CALI"]),
                     max(avg.model$m.pred.y[avg.model$Type=="CALI"])),
              ylim=c(min(avg.model$m.obs.y[avg.model$Type=="CALI"]),
                     max(avg.model$m.obs.y[avg.model$Type=="CALI"])) )
        
        # regression
        abline(lm(avg.model$m.obs.y[avg.model$Type=="CALI"] ~avg.model$m.pred.y[avg.model$Type=="CALI"] ),
               col=west.col[1], lwd=2)

        abline(lm(avg.model$m.obs.y[avg.model$Type=="VALI"] ~avg.model$m.pred.y[avg.model$Type=="VALI"] ),
               col=west.col[3], lwd=2)
        
        points(avg.model$m.pred.y[avg.model$Type=="CALI"], 
                        avg.model$m.obs.y[avg.model$Type=="CALI"],
            cex.lab=1.5, pch = 21, cex = 1, bg=west.col[1])
      
        points(avg.model$m.pred.y[avg.model$Type=="VALI"], 
             avg.model$m.obs.y[avg.model$Type=="VALI"],
             cex.lab=1.5, pch = 23, cex = 1, bg=west.col[3])
      
        legend("topleft", list.species[i], bty="n") 
      
      dev.off()
            
    } # end loop if qualy only or quanty
    
 # } # end loop variable selection
  
} # end loop species
    
# calculate avearge values for model perf 
modl.perf.av <- data.frame(cbind(mod.name.perf, modl.perf.av)) 

#write table for average perf and var imp
write.table(modl.perf.av, file=paste(inpath,"/av.",list.var.sele,"_mod.perf.txt",sep=""),sep="\t", 
            append=FALSE, row.names=FALSE,col.names=TRUE, quote=FALSE)

#####################################################################################
#####################################################################################

