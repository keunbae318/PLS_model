##PLS model###

##clean##
rm(list=ls()) 
graphics.off()

Path <- "C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/Mid_Infrared/Data"
setwd("C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/Mid_Infrared/")

#library#
library(prospectr)
library(ggplot2)
library(data.table)
library(reshape2)
library(gridExtra)
library(pls)
library(ChemoSpec)
library(caret)
library(ChemoSpec)
library(hyperSpec)
library(tidyverse)

#Source
source("Rcode/fun/misc.R")
source("Rcode/fun/C_Spectraaverage.R")
###################################################################################
##Average spectra based on the replicates##
###################################################################################

load("Data/Rawdata/MIRs.RData")


spectraMatrix=as.matrix(dataset$spc)
columns <- colnames(spectraMatrix) %in% c("4000":"650")
subsetMatrix <- spectraMatrix[, columns]
idMatrix=as.matrix(dataset$ID)

Averages <- SpecAver(subsetMatrix, 4, idMatrix, thr = .07)


#A list containing the average spectra's and ID's is constructed
spectraAverageslist           <- list() 
spectraAverageslist$ID        <- Averages$ID                      # The ID of the spectra
spectraAverageslist$Site      <- Averages$Site                    # The sitename
spectraAverageslist$Treatment <- Averages$Treatment               # The replica number of the soil core
spectraAverageslist$Depth     <- Averages$Depth                   # The upper depth of the sample
spectraAverages               <- data.frame(spectraAverageslist)  # The list is converted to a data frame
tmp                           <- data.frame(Averages$SpecMeans)   # The average spectra are copied
spectraAverages$spc           <- tmp                              # The average spectra are added
colnames(spectraAverages$spc) <- colnames(subsetMatrix) # The wavenumbers are used as column names for the spectra
#View(spectraAverages)

remove(list = c("spectraMatrix","tmp","idMatrix","subsetMatrix"))

spectraAverages$Depth=factor(spectraAverages$Depth,levels = c("T","M","B"),
                             labels = c("0-5 cm","15-20 cm","45-50 cm"))
spectraAverages$Treatment=factor(spectraAverages$Treatment,levels = c("U","R","D"),
                                 labels = c("Undrained","Rewetted","Drained"))

#####################################################################################
#GGPlot###
#####################################################################################

WSSL_hy <- new("hyperSpec",spc=as.matrix(spectraAverages$spc),
               wavelength = as.numeric(colnames(spectraAverages$spc)),
               data=spectraAverages[,1:4],
               label=list(.wavelength="nm",spc="Absorbance"))

WSSL_hy

###Plotting for whole data##
ggplot(as.long.df(WSSL_hy),aes(x=.wavelength,y=spc,color=Site,group=ID))+
  geom_line()+
  facet_grid(Depth~Treatment)+
  theme(legend.position = "top")+
  theme_bw()+
  xlab(bquote('Wavenumber '(cm^-1)))+
  ylab('Absorbance')+
  scale_x_reverse()

######################################################################################
############################
#Spectral data pre-processing
######################################################################################
for (SpectraPreProcessing in 1) {
  #------------------------------------------------------
  #Moving average
  #------------------------------------------------------
  mov_spc                        <- movav(spectraAverages$spc,11)
  tmp                            <- data.frame(mov_spc)
  colnames(tmp)                  <- colnames(mov_spc)
  spectraAvg_movav               <- subset(spectraAverages, select = -spc)
  spectraAvg_movav$spc           <- tmp
  ###################Sample_PLOT##########################
  sample_movav                   <- subset(spectraAvg_movav[1:5,])
  WSSL_hy_movav <- new("hyperSpec",spc=as.matrix(sample_movav$spc),
                       wavelength = as.numeric(colnames(sample_movav$spc)),
                       data=sample_movav[,1:4],
                       label=list(.wavelength="nm",spc="Absorbance"))
  
  grDevices::windows()
  ggplot(as.long.df(WSSL_hy_movav),aes(x=.wavelength,y=spc,color=ID,group=ID))+
    geom_line()+
    theme(legend.position = "top")+
    theme_bw()+
    xlab(bquote('Wavenumber '(cm^-1)))+
    ylab('Absorbance (MovAvg)')+
    scale_x_reverse()
  
  #------------------------------------------------------
  #Savitzky-Golay filtering
  #------------------------------------------------------
  sg_spc                         <- savitzkyGolay(spectraAverages$spc,p=3,w=21,m=0)
  tmp                            <- data.frame(sg_spc)
  colnames(tmp)                  <- colnames(sg_spc)
  spectraAvg_sg                  <- subset(spectraAverages, select = -spc)
  spectraAvg_sg$spc              <- tmp
  ###################Sample_PLOT##########################
  sample_sg                      <- subset(spectraAvg_sg[1:5,])
  WSSL_hy_sg    <- new("hyperSpec",spc=as.matrix(sample_sg$spc),
                       wavelength = as.numeric(colnames(sample_sg$spc)),
                       data=sample_sg[,1:4],
                       label=list(.wavelength="nm",spc="Absorbance"))
  
  grDevices::windows()
  ggplot(as.long.df(WSSL_hy_sg),aes(x=.wavelength,y=spc,color=ID,group=ID))+
    geom_line()+
    theme(legend.position = "top")+
    theme_bw()+
    xlab(bquote('Wavenumber '(cm^-1)))+
    ylab('Absorbance (Sav.-Golay)')+
    scale_x_reverse()
  
  #------------------------------------------------------
  #First  derivatives
  #------------------------------------------------------
  firstDer_spc                   <- t(diff(t(spectraAverages$spc), differences = 1)) 
  tmp                            <- data.frame(firstDer_spc)
  colnames(tmp)                  <- colnames(firstDer_spc)
  spectraAvg_firstDer            <- subset(spectraAverages, select = -spc)
  spectraAvg_firstDer$spc        <- tmp
  
  ###################Sample_PLOT##########################
  sample_firstDer                <- subset(spectraAvg_firstDer[1:5,])
  WSSL_hy_firstDer  <- new("hyperSpec",spc=as.matrix(sample_firstDer$spc),
                           wavelength = as.numeric(colnames(sample_firstDer$spc)),
                           data=sample_firstDer[,1:4],
                           label=list(.wavelength="nm",spc="Absorbance"))
  
  grDevices::windows()
  ggplot(as.long.df(WSSL_hy_firstDer),aes(x=.wavelength,y=spc,color=ID,group=ID))+
    geom_line()+
    theme(legend.position = "top")+
    theme_bw()+
    xlab(bquote('Wavenumber '(cm^-1)))+
    ylab('Absorbance (firstDer_spc)')+
    scale_x_reverse()
  
  #------------------------------------------------------
  #Second  derivatives
  #------------------------------------------------------
  secondDer_spc                  <- t(diff(t(spectraAverages$spc), differences = 2)) 
  tmp                            <- data.frame(secondDer_spc)
  colnames(tmp)                  <- colnames(secondDer_spc)
  spectraAvg_secondDer           <- subset(spectraAverages, select = -spc)
  spectraAvg_secondDer$spc       <- tmp
  
  ###################Sample_PLOT##########################
  sample_secondDer               <- subset(spectraAvg_secondDer[1:5,])
  WSSL_hy_secondDer  <- new("hyperSpec",spc=as.matrix(sample_secondDer$spc),
                            wavelength = as.numeric(colnames(sample_secondDer$spc)),
                            data=sample_secondDer[,1:4],
                            label=list(.wavelength="nm",spc="Absorbance"))
  grDevices::windows()
  ggplot(as.long.df(WSSL_hy_secondDer),aes(x=.wavelength,y=spc,color=ID,group=ID))+
    geom_line()+
    theme(legend.position = "top")+
    theme_bw()+
    xlab(bquote('Wavenumber '(cm^-1)))+
    ylab('Absorbance (SecondDer_spc)')+
    scale_x_reverse()
  
  #------------------------------------------------------
  #Savitzky-Golay derivatives
  #------------------------------------------------------
  sgSeconder_spc                 <- savitzkyGolay(spectraAverages$spc,p=3,w=21,m=2)
  tmp                            <- data.frame(sgSeconder_spc)
  colnames(tmp)                  <- colnames(sgSeconder_spc)
  spectraAvg_sgSeconder          <- subset(spectraAverages, select = -spc)
  spectraAvg_sgSeconder$spc      <- tmp
  ###################Sample_PLOT##########################
  sample_sgSeconder              <- subset(spectraAvg_sgSeconder[1:5,])
  WSSL_hy_sgSeconder<- new("hyperSpec",spc=as.matrix(sample_sgSeconder$spc),
                           wavelength = as.numeric(colnames(sample_sgSeconder$spc)),
                           data=sample_sgSeconder[,1:4],
                           label=list(.wavelength="nm",spc="Absorbance"))
  
  grDevices::windows()
  ggplot(as.long.df(WSSL_hy_sgSeconder),aes(x=.wavelength,y=spc,color=ID,group=ID))+
    geom_line()+
    theme(legend.position = "top")+
    theme_bw()+
    xlab(bquote('Wavenumber '(cm^-1)))+
    ylab('Absorbance (Sav.-Golay)')+
    scale_x_reverse()
  #------------------------------------------------------
  # MSC (Multiplicative Scatter Correction)
  # MSC is anther method for scatter correction
  #------------------------------------------------------
  msc_spc                        <- msc(X = as.matrix(spectraAverages$spc))
  tmp                            <- data.frame(msc_spc)
  colnames(tmp)                  <- colnames(msc_spc)
  spectraAvg_msc                 <- subset(spectraAverages, select = -spc)
  spectraAvg_msc$spc             <- tmp
  ###################Sample_PLOT##########################
  sample_msc                     <- subset(spectraAvg_msc[1:5,])
  WSSL_hy_msc   <- new("hyperSpec",spc=as.matrix(sample_msc$spc),
                       wavelength = as.numeric(colnames(sample_msc$spc)),
                       data=sample_msc[,1:4],
                       label=list(.wavelength="nm",spc="Absorbance"))
  
  grDevices::windows()
  ggplot(as.long.df(WSSL_hy_msc),aes(x=.wavelength,y=spc,color=ID,group=ID))+
    geom_line()+
    theme(legend.position = "top")+
    theme_bw()+
    xlab(bquote('Wavenumber '(cm^-1)))+
    ylab('Absorbance (MSC)')+
    scale_x_reverse()
  
  #------------------------------------------------------
  # SNV (Standard Normal Variate)
  # SNV is anther simple way for normalizing spectra that intends to correct for lgiht scatter
  #------------------------------------------------------
  snv_spc                        <- standardNormalVariate(X =spectraAverages$spc)
  tmp                            <- data.frame(snv_spc)
  colnames(tmp)                  <- colnames(snv_spc)
  spectraAvg_snv                 <- subset(spectraAverages, select = -spc)
  spectraAvg_snv$spc             <- tmp
  ###################Sample_PLOT##########################
  sample_snv                    <- subset(spectraAvg_snv[1:5,])
  WSSL_hy_snv   <- new("hyperSpec",spc=as.matrix(sample_snv$spc),
                       wavelength = as.numeric(colnames(sample_snv$spc)),
                       data=sample_snv[,1:4],
                       label=list(.wavelength="nm",spc="Absorbance"))
  
  grDevices::windows()
  ggplot(as.long.df(WSSL_hy_snv),aes(x=.wavelength,y=spc,color=ID,group=ID))+
    geom_line()+
    theme(legend.position = "top")+
    theme_bw()+
    xlab(bquote('Wavenumber '(cm^-1)))+
    ylab('Absorbance (SNV)')+
    scale_x_reverse()
  
  #------------------------------------------------------
  # SNV-Detrend (Standard Normal Variate)
  # SNV-Detrend further account for wavelength-dependent scattering effects
  #------------------------------------------------------
  dt_spc                         <- detrend(X = spectraAverages$spc, wav = as.numeric(colnames(spectraAverages$spc)))
  tmp                            <- data.frame(dt_spc)
  colnames(tmp)                  <- colnames(dt_spc)
  spectraAvg_dt                  <- subset(spectraAverages, select = -spc)
  spectraAvg_dt$spc              <- tmp
  ###################Sample_PLOT##########################
  sample_dt                    <- subset(spectraAvg_dt[1:5,])
  WSSL_hy_dt    <- new("hyperSpec",spc=as.matrix(sample_dt$spc),
                       wavelength = as.numeric(colnames(sample_dt$spc)),
                       data=sample_dt[,1:4],
                       label=list(.wavelength="nm",spc="Absorbance"))
  
  grDevices::windows()
  ggplot(as.long.df(WSSL_hy_dt),aes(x=.wavelength,y=spc,color=ID,group=ID))+
    geom_line()+
    theme(legend.position = "top")+
    theme_bw()+
    xlab(bquote('Wavenumber '(cm^-1)))+
    ylab('Absorbance (SNV Detrend)')+
    scale_x_reverse()
  
  #------------------------------------------------------
  # Centering and scaling
  # Should be applied after data transformation!
  #------------------------------------------------------
  Xscaled                        <- scale(spectraAverages$spc, center = T, scale =T)
  tmp                            <- data.frame(Xscaled)
  colnames(tmp)                  <- colnames(Xscaled)
  spectraAvg_scaled              <- subset(spectraAverages, select = -spc)
  spectraAvg_scaled$spc          <- tmp
  ###################Sample_PLOT##########################
  sample_scaled                  <- subset(spectraAvg_scaled[1:5,])
  WSSL_hy_scaled <- new("hyperSpec",spc=as.matrix(sample_scaled$spc),
                        wavelength = as.numeric(colnames(sample_scaled$spc)),
                        data=sample_scaled[,1:4],
                        label=list(.wavelength="nm",spc="Absorbance"))
  
  grDevices::windows()
  ggplot(as.long.df(WSSL_hy_scaled),aes(x=.wavelength,y=spc,color=ID,group=ID))+
    geom_line()+
    theme(legend.position = "top")+
    theme_bw()+
    xlab(bquote('Wavenumber '(cm^-1)))+
    ylab('Absorbance (Centering-Scaling)')+
    scale_x_reverse()
  
}






#############################################################################################
#############################################################################################
GCMS.dat<-read.csv("Data/Proportion_pyrolysis_V2.csv")

#----------------------------------------------------------------------
# Partial Least Squares regression (PLS)
#----------------------------------------------------------------------
str(GCMS.dat)

#----------------------------------------------------------------------
# Leave one out cross validation
#----------------------------------------------------------------------

spectrometer    = "MIR"
Compounds_group = c("Aliphatics","Alkanes","Alkenes","Benzenes","Carbohydrates","Lignins",       
                    "Methyl.ketones","N.compounds","Phenols","Unidentified")
#"Guaiacols","Syringols","Intact_lignin","Sedges","Cellulose")
PLSfactors      = 1:20
nObs            = 117   #Idunno

Nit =1;
#Create matrix to write the result

n = 1
regResults           <- data.frame(model=NA,
                                   compound=NA,
                                   RsqCV=NA, 
                                   p=NA, 
                                   RMSECV=NA, 
                                   RPD=NA, 
                                   SSE=NA,
                                   AIC=NA)

for (w in 1: length(Compounds_group)){
  #----------------------------
  #The measured index is loaded
  #----------------------------
  Index         <- Compounds_group[w]
  measured      <- GCMS.dat[,c("ID",Index)]
  
  #----------------------------------------
  #Partial Least Squares regression (PLS)
  #----------------------------------------
  number_of_preproc           = 10
  number_of_factor            = length(PLSfactors)
  number_of_total_iterations  = length(Compounds_group) * number_of_preproc * number_of_factor
  
  

  
  
  for (p in 1: number_of_preproc){
      preproc = p
      #the preprocessed spectrum is selected 
      if(preproc ==1){
        spectrum <- spectraAvg_movav
        preprocName <- "movav"
      }
      if(preproc ==2){
        spectrum <- spectraAvg_msc
        preprocName <- "msc"
      }
      if(preproc ==3){
        spectrum <- spectraAvg_dt
        preprocName <- "dt"
      }
      if(preproc ==4){
        spectrum <- spectraAvg_sg
        preprocName <- "sg"
      }
      if(preproc ==5){
        spectrum <- spectraAvg_snv
        preprocName <- "snv"
      }
      if(preproc ==6){
        spectrum <- spectraAvg_firstDer
        preprocName <- "firstDer"
      }
      if(preproc ==7){
        spectrum <- spectraAvg_secondDer
        preprocName <- "secondDer"
      }
      if(preproc ==8){
        spectrum <- spectraAvg_sgSeconder
        preprocName <- "sgSeconDer"
      }
      if(preproc ==9){
        spectrum <- spectraAvg_scaled
        preprocName <- "sclaed"
      }
      if(preproc ==10){
        spectrum <- spectraAverages
        preprocName <- "no preproc"
      }
      
      #Create the "LOO" cross-validation dataset
      calibData <- spectrum
      calibData$measured <- measured[,2]

      X                 <-  as.matrix(calibData$spc)
      Y                 <-  as.matrix(calibData$measured)
      
      pls.model         <-  mvr(Y ~ X, ncomp = 20, method = "kernelpls", validation = "LOO")
      
      

        
      for (ci in 1: number_of_factor){
        
        c = PLSfactors[ci]
        tableStorageName = paste(spectrometer,"_",Index,"_",preprocName, c, sep="")
        
        #-validation-
        results <- as.data.frame(plot(pls.model, ncomp = c, asp = 1, line = TRUE))
        
      
      
    
        #linear regression between measured and calculated index
          reg      <- lm(results$predicted ~results$measured)
          Rvalue   <- summary(reg)$r.squared
          pvalue   <- summary(reg)$coefficient[2,4]
      
        #calculate RMSE, RPD, and AIC
          temp  <-  (results$measured -results$predicted)^2
          RMSE  <-  (mean(temp, na.rm = TRUE))^0.5
          RPD   <-  sd(results$measured)/RMSE
          SSE   <-  sum((results$measured - results$predicted)^2, na.rm =TRUE)
          AIC   <-  nObs*log(SSE/nObs)+2*c
      
        #the R and p-value for each regression combination are stored
          
          
          regResults[n,1] <- tableStorageName
          regResults[n,2] <- Index
          regResults[n,3] <- Rvalue
          regResults[n,4] <- pvalue
          regResults[n,5] <- RMSE
          regResults[n,6] <- RPD
          regResults[n,7] <- SSE
          regResults[n,8] <- AIC
          
          
          
          n   = n+1
          Nit = Nit+1
      
      # #write table with results for each calibration
      # Name      = paste(Index, format(nObs, digits =2), sep ="")
      # setwd("C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/Mid_Infrared/Data/PLS_model/")
      # dir.create(Name, showWarnings = FALSE)
      # tableName = paste(Path,"/PLS_model/",Name,"/",tableStorageName,".txt",sep ="")
      # write.table(results, file = tableName,row.names = FALSE)
      # 
      Percent   = paste(format((Nit/number_of_total_iterations)*100,digits=4), "%",sep ="")
      print(paste(tableStorageName," (",Percent,"complelted )", sep =""))
    }
  }
  
}

#write.csv(regResults, "C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/Mid_Infrared/PLS_model/LOO_result.csv")
PLS.result.dat <- read.csv("PLS_model/LOO_result.csv")

lowest_aic_models <- PLS.result.dat %>%
  group_by(compound) %>%
  slice(which.min(AIC))

View(lowest_aic_models)
  
  


