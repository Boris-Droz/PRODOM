####################################################################################################################################
###
### PARAFAC ANALYSIS -- ROBUST MODEL
#####============================================================
# Historic - version- authors
##############################
# v1.0 - Elena Fernandez December 2021
# v2.0 - Modif by Boris Droz - June 2022
# v3.0 - 

# PRODOM project - School of Biological, Earth and Environmental Sciences
# Environmental Research Institute (ERI)
# University College Cork
#################################################################################
## Descriptif
##############
# Performed PARAFAC model from corrected EEM following the tutorial of Pucher 2008.

## REFERENCE
############
## Pucher, M.; W?nsch, U.; Weigelhofer, G.; Murphy, K.; Hein, T.; Graeber, D., 
##        staRdom: Versatile software for analyzing spectroscopic data of dissolved organic matter in R. Water 2019, 11, (11), 2366.
## Stedmon, C. A.; et al. Limnol. Oceanogr. Methods 2008, 6, 572–579.

## Load library
##################
# work under R-4.2.0 [64bits]
# Script PARAFAC -- install in the first time use
#install.packages("devtools")
#devtools::install_github("MatthiasPucher/staRdom")
# install.packages("matlab")

library(eemR) # v1.0.1
library(ggplot2) #3.3.6
library(parallel) # 4.2.0
library(dplyr) # 1.0.9
library(tidyr) # 1.2.0
library(staRdom) # 1.1.25
library(plotly) # 4.10.0
library("matlab") # 1.0.4
library(devtools) #2.4.3

## set parameter
################
# Setting working directory
workdir <- "C:/Users/bodro/Documents/Actual/EPA 2019-W-MS43 - PRODOM/Data_analysis/PARAFAC analysis"

ctol <- 10^-8 # decrease tolerance in PARAFAC analysis
nstart = 1000 # increase number of random starts
maxit = 10000 # increase number of maximum interations

n.comp <- 6 ## number of component 

exc.list <- "YES" # Yes or no if a list of outlier as been created

################################################################################################################################################
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
##########################################################################
##########################################################################
########################### PARAFAC ANALYSIS ###########################
#################################################################################
#################################################################################
setwd(workdir)
# self-creat a PARAFAC folder
creat.subDir (paste(workdir,"/output",sep=""),"PARAFAC_MODEL")
creat.subDir (paste(workdir,"/output/PARAFAC_MODEL",sep=""),"residual_plot")

# cores <- parallel::detectCores(logical = FALSE) # detect available parrallel processing

################################################################################
#### 1) Loading of our DATA ####
################################
sample_dir <- paste(workdir,"/output/EEM_corr/",sep="")
eem_list <- eem_read(sample_dir, recursive = TRUE, import_function = eem_csv)

if( exc.list == "YES"){ 
    # exclusion list
    fns <- paste(workdir,"/output/EXPLOR_PARAFAC/exclude.list.txt",sep="")
    d.exclu <- read.table(fns, header=FALSE)
    
    exclude <- list("ex" = c(),
                    "em" = c(),
                    "sample" = unlist(as.vector(d.exclu))
                    )
    
    eem_list_ex <- eem_exclude(eem_list, exclude)
    }else{}

################################################################################
#### 2) Run model ####   
######################

pf4 <- eem_parafac(eem_list_ex, comps = n.comp, normalise = TRUE,
                    const = c("nonneg", "nonneg", "nonneg"),
                    maxit = maxit, nstart = nstart, ctol = ctol, 
                    output = "all", cores = cores)

pf4 <- lapply(pf4, eempf_rescaleBC, newscale = "Fmax")

# save the PARAFAC model - R file
file.name <- paste(workdir,"/output/PARAFAC_MODEL/PARAFAC_Rmodel",sep='')
save(list=c('pf4'),file=file.name)

# Check the convergence behaviour of the created models
sink(paste(workdir,"/output/PARAFAC_MODEL/model.conv.txt", sep=""))
  print(eempf_convergence(pf4[[1]]))
  print(paste("number of emm used: ",length(eem_list_ex),sep="") ) #number of used eem
sink()

# Check out the loading of each component 
## see recomendation in Stedmon, C. A.; et al. Limnol. Oceanogr. Methods 2008, 6, 572–579.
png(filename = paste(workdir,"/output/PARAFAC_MODEL/ex-em.spectra_",
                     n.comp,"-comp.png",
                     sep=""), width = 520, height = 520)
        eempf_compare(pf4, contour = TRUE)  
dev.off()

# check leverage to see if all sample are within the range wanted 
## --> check for oulier
png(filename = paste(workdir,"/output/PARAFAC_MODEL/leverage_",
                     n.comp,"-comp.png", sep=""),
                     width = 520, height = 520)
        eempf_leverage_plot(eempf_leverage(pf4[[1]])) # [[1]] means the 4th model in the list, 6 component model in that case
dev.off()

lev.samp <- eempf_leverage(pf4[[1]])
  
png(filename = paste(workdir,"/output/PARAFAC_MODEL/cor.plot_",
                     n.comp,"-comp.png", sep=""), 
                      width = 520, height = 520)
      eempf_corplot(pf4[[1]], progress = FALSE)
dev.off()
################################################################################
## 3) Plot the resulting components and loadings
# loading -- show the fluorescence differences of different 
#  components in different samples.
####################################

png(filename = paste(workdir,"/output/PARAFAC_MODEL/",
                     n.comp,"-comp_emm.png", sep=""), 
                width = 520, height = 520)
  print( eempf_comp_load_plot(pf4[[1]], contour = TRUE)[[1]] )
dev.off()

png(filename = paste(workdir,"/output/PARAFAC_MODEL/samp.load.png", sep=""), 
    width = 520, height = 520)
print( eempf_comp_load_plot(pf4[[1]], contour = TRUE)[[2]] )
dev.off()

# export a Fmax table
pf.load <- norm2A(pf4[[1]])$A

write.table(pf.load, file=paste(workdir,"/output/PARAFAC_MODEL/Fmax",n.comp,"-comp_model.txt",sep=""),
            sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)

# loading table in relative ratio of total comp.
pf.load.rel <- pf.load/apply(pf.load,1,sum)

write.table(pf.load.rel, file=paste(workdir,"/output/PARAFAC_MODEL/Prop.Comp",n.comp,"-comp_model.txt",sep=""),
            sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)

################################################################################
#### 4) Residuals and error of the model ####
## plot residual of all sample for inspection
## see recomendation in Stedmon, C. A.; et al. Limnol. Oceanogr. 
## Methods 2008, 6, 572–579.
############################################################

for (i in 1:length(eem_list_ex))
      {
        png(filename = paste(workdir,"/output/PARAFAC_MODEL/residual_plot/",
                             eem_list_ex[[i]]$sample,".png", sep=""), 
                      width = 520, height = 520)
      
          print(eempf_residuals_plot(pf4[[1]],eem_list_ex,spp=6,
                           select= eem_list_ex[[i]]$sample,
                           residuals_only = TRUE, cores = cores, 
                              contour = FALSE) )# 6-comp model
      dev.off() 
      }

# claculate residual
emm.res <- eempf_residuals(pf4[[1]], eem_list_ex, cores = cores)

write.table(emm.res, file=paste(workdir,"/output/PARAFAC_MODEL/all.residual.txt",sep=""),
            sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)


## emm.res <- read.table(choose.files())
## percentage residual from origin
res <- emm.res[emm.res$type=="residual",]
ori <- emm.res[emm.res$type=="sample",]

per.resi <- res$value/ori$value*100
length(per.resi)
#per.resi <- abs(per.resi)
per.resi <-per.resi[per.resi>-Inf]

min(per.resi,na.rm=TRUE)
max(per.resi,na.rm=TRUE)
mean(per.resi,na.rm=TRUE)
sd(per.resi,na.rm=TRUE)

## average per sample
emm.res <- emm.res[emm.res$type=="residual",]

mean.res <- aggregate( value ~ sample, emm.res, mean)
sd.res <- aggregate( value ~ sample, emm.res, sd)

sampl.resi <- cbind(mean.res,sd.res[,2])

colnames(sampl.resi)<-c("sample","resi.mean", "resi.sd")

write.table(sampl.resi, file=paste(workdir,"/output/PARAFAC_MODEL/residual.sampl.txt",sep=""),
            sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)

summary.res <- c(mean(emm.res$value),sd(emm.res$value)) # summary of the residual  
#print(summary.res)

residual.sampl <- eempf_residuals_metrics(emm.res,lev.samp)$sample

write.table(residual.sampl, file=paste(workdir,"/output/PARAFAC_MODEL/resi.metric.sampl.txt",sep=""),
            sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)

summary.res.metric <- rbind(apply(residual.sampl[,2:6],2, mean), apply(residual.sampl[,2:6],2, sd) )
################################################################################  
#### 5) Split-half analysis ####
# Intended to show the stability of the model (takes time)
##########################################################

#calculate split_half analysis
sh <- splithalf(eem_list, n.comp, normalise = TRUE, rand = TRUE, 
                cores = cores, nstart = nstart, maxit = maxit, ctol = ctol)

# We plot the results. Our model is stable, if the graphs look similar.
splithalf_plot(sh)

# We check Tucker's congruency coefficients
tcc_sh_table <- splithalf_tcc(sh)

write.table(tcc_sh_table, file=paste(workdir,"/output/PARAFAC_MODEL/tcc_table.txt",sep=""),
            sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)

summary.tcc <- c(mean(c(tcc_sh_table$tcc_em,tcc_sh_table$tcc_ex)),
                 sd(c(tcc_sh_table$tcc_em,tcc_sh_table$tcc_ex)) )
# print(summary.tcc)

# write summary table
sum.table <- rbind(summary.res, summary.tcc,t(summary.res.metric) )
colnames(sum.table) <- c("mean","sd")
row.names(sum.table)[1:2] <- c("residual","TCC")

write.table(sum.table, file=paste(workdir,"/output/PARAFAC_MODEL/summary_table.txt",sep=""),
            sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)

################################################################################
#### 6) Calculat the importance of each component ####
#######################################################
# Importance of components into the model
varimp <- eempf_varimp(pf4[[1]], eem_list_ex, cores = cores) # relative values

varimp <- varimp/sum(varimp) # importance relative

write.table(varimp, file=paste(workdir,"/output/PARAFAC_MODEL/imp.comp.txt",sep=""),
            sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)

################################################################################   
#### 7) Model export and interpretation ####
##############################################
# export your model for further Comparing your data using openfluor.org
eempf_openfluor(pf4[[1]], file = paste(workdir,"/output/PARAFAC_MODEL/my_modelforopenfluor.txt", sep="") )

# Exporting the model
# eempf_export(pf4[[1]], export = "PARAFAC_model_export.txt", Fmax = TRUE)

#################################
## parafac max peak table of the component!!!
peak.em <- as.numeric(row.names(pf4[[1]]$B) [apply(pf4[[1]]$B, 2, which.max)] )
peak.ex <- as.numeric(row.names(pf4[[1]]$C) [apply(pf4[[1]]$C, 2, which.max)] )

parafac.comp.peak <- rbind(peak.em,peak.ex)
colnames(parafac.comp.peak) <- paste(rep("Comp.",length(peak.em) ),
                                         seq(from=1, to= length(peak.em),by=1),
                                         sep="")

write.table(parafac.comp.peak, file=paste(workdir,"/output/PARAFAC_MODEL/peak.parfac.comp.txt",sep=""),
            sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)

###########################################################################################
###########################################################################################
###########################################################################################
###### GO FOR COFFE , maybe a dinner
#######################
#######################
