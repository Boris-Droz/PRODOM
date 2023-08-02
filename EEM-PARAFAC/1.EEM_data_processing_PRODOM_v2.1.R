####################################################################################################################################
###
### Excitation-emission matrice fluorescence spectra
##  preprocessing and indices calculation with staRdom Package 
####
#####============================================================
## Historic - version- authors
##############################
# v1.0 - Elena Fernandez December 2021
# v2.0 - Modif by Boris Droz - June 2022

# PRODOM project - School of Biological, Earth and Environmental Sciences
# Environmental Research Institute (ERI)
# University College Cork
#################################################################################
## REFERENCE
############
## Pucher, M.; W?nsch, U.; Weigelhofer, G.; Murphy, K.; Hein, T.; Graeber, D., 
## staRdom: Versatile software for analyzing spectroscopic data of dissolved 
## organic matter in R. Water 2019, 11, (11), 2366. 
## https://doi.org/10.3390/w11112366.

## Descriptif
##############
# script perform the following data processing prior to parafac analysis:
# Blank subtraction
# Scattering removal and interpolation of missing data
# Inner filter effect correction
# Raman normalization
# Wavelengths cut
# Interpolation
# Export processed EEMs to .cvs format
# Extraction of fluorescence peaks and indices
# Slope parameter calculation
# VLux fluorescence peaks calculation 
## ->> self made peak to compare with the sensorVLux Tpro from Chelsea tech.
###############################################################################

## NEEDS ###
############
# input data are in two folder at the workdir roots
# EEM data and ABS data are oragnized per folder (sampling time) and should have the same name
# Name of data are the follow "location_date"
# One blank (deio water) is needed per each EEM folder with a name "Blank_"
# the two folder needed are:
# - /input/EEM
# - /input/ABS

# self-creat a output folder with "EEM-corr" spectra.
# output folder will contain peak picking and indice and plot smooth data

## Load library version and requirement.
#######################################
# work under R-4.2.0 [64bits]

# Script PARAFAC -- install in the first time use
#install.packages("devtools")
#devtools::install_github("MatthiasPucher/staRdom")
# install.packages("matlab")

library(devtools) #2.4.3
library(dplyr) # 1.0.9
library(eemR) # v1.0.1
library(ggplot2) #3.3.6
library(magrittr) # 2.0.3
library(parallel) # 4.2.0
library(plotly) # 4.10.0
library(rlang) #1.0.2
# library(Rtools) --> should be install at the roots C:/rtools42
library(staRdom) # 1.1.25
library(tidyr) # 1.2.0
library(usethis) # 2.1.6

## set parameter
################
# Setting working directory
workdir <- "C:/Users/bodro/Documents/Actual/EPA 2019-W-MS43 - PRODOM/Data_analysis/PARAFAC analysis"

# range of eem interest
em_range = c(300, 600) # e.g. c(300,500), c(0,Inf) to use everything
ex_range = c(250, 460) # e.g. c(300,500), c(0,Inf) to use everything

cuv.siz = 1 # size of the cell for the spectro

remove_scatter_width <- c(15,10,16,12) # defaut parameter of the function

RtoQ<- 0.0085 #factor of Raman to quinine ratio specific to instrument

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
#################################################################################
########################### Preprocessing EEM spectra ###########################
#################################################################################
#################################################################################
setwd(workdir)
# self-creat a EEM-corrected folder
creat.subDir (paste(workdir,"/output",sep=""),"EEM_corr")
creat.subDir (paste(workdir,"/output",sep=""),"plot_EEM_corr")

cores <- detectCores(logical = FALSE) # define parrallel processing

# make folder list (should be similar in EEM and ABS folder)
list.folder <- list.dirs(path = paste(workdir,"/input/EEM",sep=""), full.names = FALSE, recursive = FALSE)

d.out <- NULL
for ( k in 1:length(list.folder))
    {
  
  cat("PROCESS",list.folder[k], "\n",append = F) 
    ################################################################################
    #### Read EEM data ####
    sample_dir <- paste(workdir,"/input/EEM/",list.folder[k],sep="")
    eem_list <- eem_read(sample_dir, recursive = TRUE, import_function = eem_csv)
    
    ## Read absorbance data ####
    absorbance_dir = paste(workdir,"/input/ABS/",list.folder[k],sep="")
    abs_data <- absorbance_read(absorbance_dir, verbose = FALSE)
    
    # blank subtraction
    eem_list <- eem_remove_blank(eem_list)
    
    ## inner filtration effect correction for 1 cm cuvette 
    eem_list <- eem_ife_correction(eem_list,abs_data, cuvl = cuv.siz)
    
    # raman normalization
    eem_list <- eem_raman_normalisation2(eem_list, blank = "blank")
    
    # remove blank from dataset
    eem_list <- eem_extract(eem_list, c("blank"),ignore_case = TRUE)
    
    # Remove and interpolate scattering
    remove_scatter <- c(TRUE, TRUE, TRUE, TRUE)
    
    eem_list <- eem_rem_scat(eem_list, remove_scatter = remove_scatter, 
                             remove_scatter_width = remove_scatter_width)
    
    # interpolate data - recommended to increase parafac co
    # (1) spline interpolation with mba.points (Lee, Wolberg, and Shin 1997)
    eem_list <- eem_interp(eem_list, cores = cores, type = 1, extend = FALSE)
    
   # eem_overview_plot(eem_list, spp=1, contour = TRUE) # plot one by one
    
    # smooth datat not --> not appropri for PARAFAC
   #eem4peaks <- eem_smooth(eem_list, n = 16)
    #eem4peaks <- eem4peaks %>% eem_range(ex=c(240,400),em=c(300,550))
    #mycolor <- colorRampPalette(brewer.pal(8, "Spectral"))(200)
    #eem_overview_plot(eem4peaks, spp=1, contour = TRUE, colpal= rev(mycolor)) # plot one by one
    
    ## Wavelengths cut ####
    eem_list <-eem_list %>% eem_range(ex=ex_range,em=em_range)
    
    for (n in 1:length(eem_list))
        {
          # save eem process one by one
          row_names <- c(NA,eem_list[[n]]$ex)
          col_names <- eem_list[[n]]$em
          
          Sample_1_col <- eem_list[[n]][["x"]] 
          Sample_1_row <- cbind (col_names, Sample_1_col)
          Sample_1 <- rbind(row_names, Sample_1_row)
          write.table(Sample_1, file=paste(workdir,"/output/EEM_corr/",eem_list[[n]]$sample,".csv" ,sep=""),
                                   sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
        }
    # save the plot one by one in a date folder
    creat.subDir (paste(workdir,"/output/plot_EEM_corr",sep=""),list.folder[k])
    
    png(filename = paste(workdir,"/output/plot_EEM_corr/",list.folder[k],"/Rplot%03d.png",
                         sep=""), width = 480, height = 480)
    
        eem_overview_plot(eem_list, spp=1, contour = TRUE)   
  
    dev.off()
    
    #### Peak picking and indices
    eem4peaks <- eem_smooth(eem_list, n = 4)
    
    bix <- eem_biological_index(eem4peaks)
    coble_peaks <- eem_coble_peaks(eem4peaks)
    fi <- eem_fluorescence_index(eem4peaks)
    hix <- eem_humification_index(eem4peaks, scale = TRUE)
    
    indices_peaks <- coble_peaks %>%
      full_join(bix, by = "sample") %>%
      full_join(fi, by = "sample") %>%
      full_join(hix, by = "sample")

    slope_parms <- abs_parms(abs_data, cuvl = 1, cores = cores, 
                             unit= c("absorption"))
    
    #### 19. VLux fluorescence peaks calculation ###
    interpolation <- TRUE
    smooth = 4 
    if(smooth != FALSE & interpolation) 
    {eem4peaks <- eem_list %>% eem_smooth(n=smooth)} else {eem4peaks <- eem_list}
    
    tlf.eem_qsu <-NULL
    hlf.eem_qsu <-NULL
    for (q in 1:length(eem4peaks))
        {
          eem4peaks[[q]][["x"]] <- eem4peaks[[q]][["x"]]/RtoQ ## convert data in QSU
          ## index m is instrument specific see detail in
          ##  Murphy, K. R.; et al. Environ. Sci. Technol. 2010, 44 (24), 9405–9412. https://doi.org/10.1021/es102362t.
          tlf.eem_qsu <-c(tlf.eem_qsu, max(pracma::interp2(eem4peaks[[q]][["ex"]], eem4peaks[[q]][["em"]], 
                                             eem4peaks[[q]][["x"]], rep(280, length(340:390)), 340:390)) )
          
          hlf.eem_qsu <-c(hlf.eem_qsu, max(pracma::interp2(eem4peaks[[q]][["ex"]], eem4peaks[[q]][["em"]], eem4peaks[[q]][["x"]], 
                                             rep(280, length(425:475)), 425:475)) )
        }
    
    d.out <- rbind(d.out, cbind(indices_peaks,slope_parms[,2:ncol(slope_parms)], 
                                tlf.eem_qsu=tlf.eem_qsu,hlf.eem_qsu=hlf.eem_qsu ))
} # end list folder  

# data table of all values
write.table(d.out, file=paste(workdir,"/output/fluo_indice.txt",sep=""),
            sep="\t", append=FALSE, row.names=TRUE,col.names=TRUE, quote=FALSE)

###########################################################################################
###########################################################################################
###########################################################################################
###### GO FOR COFFE
#######################
#######################

