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
library(ggpubr)
#Source
source("Rcode/fun/misc.R")
source("Rcode/fun/C_Spectraaverage.R")
source("Rcode/Reference/C_Areabased_Normalization.R")
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

#####################################################################################
#####################################################################################

spectraAverages$spc <- savitzkyGolay(spectraAverages$spc,p=3,w=21,m=0)
spectraAverages$spc <- Normalization100(spectraAverages$spc)

#dataset4$spc2 <- binning(dataset4$spc2, bin.size = 10)

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

#############################################################################################
#############################################################################################
#Labile_C_prediction##
# 
Cummulative.dat=read.csv("C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/CO2_Respiration/Rawdata/Cummulativeflux_105.csv")
EA.dat = read.csv("C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/Elemental_Analysis/Rdata/EA_data.csv")

Cummulative.dat$Flux <- Cummulative.dat$Flux/(EA.dat[EA.dat$Depth %in% c("Top","Bottom"), ]$C*0.01)

PLSsubset.dat <- subset(spectraAverages, Depth %in% c("0-5 cm","45-50 cm")) %>% 
  mutate(flux = Cummulative.dat$Flux) %>% 
  select(1:4,6,5) %>% 
  subset(Depth == "0-5 cm")


########################################################################################
##########################################################################################
#PLS_model#
PLS_flux.model <- plsr(flux ~ spc, ncomp =10, data=PLSsubset.dat,validation ="LOO", scale= FALSE)
ncomp.permut <- selectNcomp(PLS_flux.model, method = "randomization", plot = TRUE)
ncomp.onesigma <- selectNcomp(PLS_flux.model, method = "onesigma", plot = TRUE)

summary(PLS_flux.model)
plot(PLS_flux.model, which = "validation")


plot(RMSEP(PLS_flux.model), legendpos = "topright")
plot(PLS_flux.model, ncomp = 1, asp = 1, line = TRUE)
explvar(PLS_flux.model)

plot(PLS_flux.model, plottype = "scores", comps = 1:4)
plot(PLS_flux.model, "loadings", comps = 1:2, legendpos = "topleft", labels = "numbers", xlab = "nm")


a <- as.data.frame(plot(PLS_flux.model, ncomp = 3, asp = 1, line = TRUE))
b <- as.data.frame(predict(PLS_flux.model,ncomp = 3, newdata =PLSsubset.dat$spc))
a$predicted2 <- b$`flux.3 comps`

stats_result <- data.frame(n=39,
                           nf=3,
                           Rsq=1-(sum((a$predicted2-a$measured)^2)/sum((a$measured-mean(a$measured))^2)),
                           RMSE=(mean((a$measured -a$predicted2)^2))^0.5,
                           Rsq_CV= 1-(sum((a$predicted-a$measured)^2)/sum((a$measured-mean(a$measured))^2)),
                           RMSE_CV=(mean((a$measured -a$predicted)^2))^0.5,
                           RPD=sd(a$measured)/((mean((a$measured -a$predicted)^2))^0.5)
) 


a <- a %>% 
  gather(predicted:predicted2, key= test, value= Predicted) %>% 
  mutate(measured = measured) %>% 
  mutate(test = factor(test, levels = c("predicted","predicted2"),
                       labels =c("Calibration","Validation"))) %>%
  data.frame() %>% 
  mutate(n    = stats_result$n,
         nf   = stats_result$nf,
         Rsq  = stats_result$Rsq,
         RMSE = stats_result$RMSE,
         Rsq_CV = stats_result$Rsq_CV,
         RMSE_CV =stats_result$RMSE_CV,
         RPD= stats_result$RPD) %>% 
  mutate_if(is.numeric,round,digits=2)


Pls_plot<-ggplot(aes(measured, Predicted, group =test, color = test, shape = test), data = a) +
  stat_smooth(method='lm', 
              formula = y ~ x,
              se = F)+
  geom_point(size=2)+
  geom_abline(linetype="dashed", color="darkgray",size =1.5)+
  scale_color_grey()+
  theme_bw()+
  theme(legend.position = "top",
        legend.title = element_blank())+
  labs(x = expression(Measured~cummulative~flux~on~day~105*phantom(x)*(mu*g~CO[2]-C~g~soil^{-1})),
       y = expression(Predicted~cummulative~flux~on~day~105*phantom(x)*(mu*g~CO[2]-C~g~soil^{-1})))+
  geom_text(aes(x = Inf, y = -Inf, label = paste("n =", n)), vjust = -13, hjust = 1, size = 5,show_guide = F)+
  geom_text(aes(x = Inf, y = -Inf, label = paste("nf =", nf)), vjust = -11, hjust = 1, size = 5,show_guide = F)+
  geom_text(aes(x = Inf, y = -Inf, label = paste("Rsq =", Rsq)), vjust = -9, hjust = 1, size = 5,show_guide = F)+
  geom_text(aes(x = Inf, y = -Inf, label = paste("RMSE =", RMSE)), vjust = -7, hjust = 1, size = 5,show_guide = F)+
  geom_text(aes(x = Inf, y = -Inf, label = paste("Rsq_CV =", Rsq_CV)), vjust = -5, hjust = 1, size = 5,show_guide = F)+
  geom_text(aes(x = Inf, y = -Inf, label = paste("RMSE_CV =", RMSE_CV)), vjust = -3, hjust = 1, size = 5,show_guide = F)+
  geom_text(aes(x = Inf, y = -Inf, label = paste("RPD =", RPD)), vjust = -1, hjust = 1, size = 5,show_guide = F)


Pls_plot

plot(PLS_flux.model, plottype = "coef", ncomp=1:3, legendpos = "bottomleft",
     labels = "numbers", xlab = "nm")
coefficent.dat <- as.data.frame(PLS_flux.model$coefficients)
coefficent2.dat <- coefficent.dat %>% 
  mutate(Wavenumer = rownames(.),
         Coefficient = coefficent.dat$`flux.4 comps`) %>% 
  select(Wavenumer,Coefficient)




loading.spc <- as_tibble(coefficent2.dat[,1:2]) 

loading.spc$Wavenumer <- as.numeric(loading.spc$Wavenumer)
loading.spc$axis <- "Regression coefficient"
# loading.spc$PC_axis <- factor(loading.spc$PC_axis, levels = c("comp1","comp2"),
#                               labels = c("Latent variable 1","Latent variable 2"))

PLS_coefficients_Plot<-ggplot(loading.spc,aes(x=Wavenumer,y=Coefficient, color = axis))+ 
  geom_line(size=0.8)+
  theme_bw()+
  theme(legend.position = "top",
        legend.title =element_blank(),
        panel.grid.major.x = element_line(color = "black",
                                          size = 0.1,
                                          linetype = 2),
        panel.grid.major.y = element_line(color = "grey",
                                          size = 0.1,
                                          linetype = 2))+
  xlab(bquote('Wavenumber '(cm^-1)))+
  ylab('PLS regression coefficients')+
  scale_x_reverse(n.breaks = 10)+
  grids(linetype = "dashed")+
  scale_color_manual(values=c("black"))


ggarrange(Pls_plot, PLS_coefficients_Plot, ncol = 2, labels = c("A", "B"), widths = c(1,1.5))


Inspection.dat<-ggplot(loading.spc,aes(x=Wavenumer,y=Coefficient, color = axis))+ 
  geom_line(size=0.8)+
  theme_bw()+
  scale_x_reverse(n.breaks = 10)+
  grids(linetype = "dashed")+
  scale_color_manual(values=c("black"))

ggplotly(Inspection.dat)




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
# Choose spectrometer and type of iron-data to rate
#----------------------------------------------------------------------

spectrometer    = "MIR"
Compounds_group = c("Aliphatics","Alkanes","Alkenes","Benzenes","Carbohydrates","Lignins",       
                    "Methyl.ketones","N.compounds","Phenols","Unidentified")
#"Guaiacols","Syringols","Intact_lignin","Sedges","Cellulose")
PLSfactors      = 1:20
nObs            = 39   #Idunno

Nit =1;
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
  
  #Create matrix to write the result
  number_of_iterations <- number_of_preproc
  regResults           <- data.frame(data=NA,Rsquared=NA, p=NA, RMSD=NA, RPD=NA, SSE=NA,AIC=NA)
  n                    <- 1
  
  for (p in 1: number_of_preproc){
    for (ci in 1: number_of_factor){
      c = PLSfactors[ci]
      
      #the preprocessed spectrum is selected 
      preproc = p 
      
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
      tableStorageName = paste(spectrometer,"_",Index,"_",preprocName, c, sep="")
      
      #A subset of the spectra with measure 
      Xdata     <- spectrum
      XdataSort = Xdata
      # Sort measured index data
      names = colnames(measured)
      names = names[2]
      measured = measured[order(measured[[names]]),]
      
      # This subset is sorted in the same way the measured index is ordered 
      Measured = measured
      names = as.character(measured$ID)
      for(s in 1:nrow(measured)){
        nameTmp       <- measured$ID[s]
        row           <- match(nameTmp, Xdata$ID)
        XdataSort[s,] <- Xdata[row,]
      }
      
      #Part of this dataset is stored for calibration
      calib   = seq(1,117, by=3)
      
      #Create the calibration dataset
      calibData <- XdataSort[-c(calib),]
      
      #The rest is kept for validation
      validData <- XdataSort[c(calib),]
      
      #In this subset, an additional variable is created in which the measured index are stored
      calibData$meas    =  NA
      validData$meas    =  NA
      totalSpectraCalib =  nrow(calibData)
      totalSpectraValid =  nrow(validData)
      
      #a for loop to include the measured index into the dataframe with th spectra for the calibration set
      for (i in 1: totalSpectraCalib){
        nameTmp            <- calibData$ID[i]
        row                <- match(nameTmp, Measured$ID)
        calibData$meas[i]  =  measured[[Index]][row]
      }
      for (i in 1: totalSpectraValid){
        nameTmp            <- validData$ID[i]
        row                <- match(nameTmp, Measured$ID)
        validData$meas[i]  =  Measured[[Index]][row]
      }
      
      # PLS model
      #-calibration-
      X                 <-  as.matrix(calibData$spc)
      Y                 <-  as.matrix(calibData$meas)
      pls.model         <-  mvr(Y ~ X, ncomp = c, method = "kernelpls")
      loadingweight     =   loading.weights(pls.model)
      XScores           =   scores(pls.model)
      XLoading          =   loadings(pls.model)
      YScores           =   Yscores(pls.model)
      YLoading          =   Yloadings(pls.model)
      ExpVar            =   explvar(pls.model)
      
      #-validation-
      restSpectra       <-  as.matrix(validData$spc)
      pls.pred          <-  predict(pls.model, newdata = restSpectra, ncomp =c)
      new_values        <-  as.matrix(pls.pred)
      
      results           <-  subset(validData, select = -spc) 
      results$predicted <-  NA
      results$predicted <- new_values
      
      #linear regression between measured and calculated index
      reg      <- lm(results$predict ~results$meas)
      Rvalue   <- summary(reg)$r.squared
      pvalue   <- summary(reg)$coefficient[2,4]
      ytextR   =  min(results$predicted, na.rm = TRUE)+((max(results$predicted, na.rm = TRUE)
                                                         -min(results$predicted, na.rm = TRUE))*(9/10))
      xtextR   =  min(results$meas, na.rm = TRUE)+((max(results$meas, na.rm = TRUE)
                                                    -min(results$meas, na.rm = TRUE))*(1/5))
      ytextP   =  ytextR
      xtextP   =  xtextR
      
      #calculate RMSE, RPD, and AIC
      temp  <-  (results$meas -results$predicted)^2
      RMSD  <-  (mean(temp, na.rm = TRUE))^0.5
      RPD   <-  sd(results$meas)/RMSD
      SSE   <-  sum((results$meas - results$predicted)^2, na.rm =TRUE)
      AIC   <-  nObs*log(SSE/nObs)+2*c
      
      #the R and p-value for each regression combination are stored
      
      regResults[n,1] <- tableStorageName
      regResults[n,2] <- Rvalue
      regResults[n,3] <- pvalue
      regResults[n,4] <- RMSD
      regResults[n,5] <- RPD
      regResults[n,6] <- SSE
      regResults[n,7] <- AIC
      
      n   = n+1
      Nit = Nit+1
      
      #write table with results for each calibration
      Name      = paste(Index, format(nObs, digits =2), sep ="")
      setwd("C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/Mid_Infrared/Data/PLS_model/")
      dir.create(Name, showWarnings = FALSE)
      tableName = paste(Path,"/PLS_model/",Name,"/",tableStorageName,".txt",sep ="")
      write.table(results, file = tableName,row.names = FALSE)
      
      Percent   = paste(format((Nit/number_of_total_iterations)*100,digits=4), "%",sep ="")
      print(paste(tableStorageName," (",Percent,"complelted )", sep =""))
    }
  }
  Path
  #write the results to textfile
  Name      = paste(Index, format(nObs, digits =2), sep ="")
  tableName = paste(Path,"/PLS_model/",Name,"/",spectrometer,Index,"diffFactors.txt", sep = "")
  write.table(regResults, file= tableName,row.names = FALSE)
}








all_new_data <- list()
for(mirs in 1:length(Compounds_group)){
  
  
  #data import
  Index <- Compounds_group[mirs]
  data_import <- read.table(paste0(Index,"39/MIR",Index,"diffFactors.txt"),header = T)
  
  #best_model_selection
  model_name <- data_import[data_import$AIC == sort(data_import$AIC)[1],]$data
  
  #import plot
  validation.dat <-read.table(paste0(Index,"39/",model_name,".txt"),header=T)
  validation.dat$Index <- Compounds_group[mirs]
  split_strings <- strsplit(model_name, "_")
  validation.dat$model_name <- sapply(split_strings,'[',3)
  validation.dat$Rsquared <- data_import[data_import$AIC == sort(data_import$AIC)[1],]$Rsquared
  validation.dat$RMSD <- data_import[data_import$AIC == sort(data_import$AIC)[1],]$RMSD
  validation.dat$AIC <- data_import[data_import$AIC == sort(data_import$AIC)[1],]$AIC
  
  
  all_new_data[[mirs]] <- validation.dat
}

final_new_data <- do.call(rbind, all_new_data) %>% 
  mutate_if(is.numeric,round,digits=2)

final_new_data$Index <- factor(final_new_data$Index,
                               levels = c("Aliphatics", "Alkanes","Alkenes","Benzenes","Carbohydrates",
                                          "Lignins","Methyl.ketones","N.compounds","Phenols","Unidentified"),
                               labels = c("Aliphatics", "Alkanes","Alkenes","Benzenes","Carbohydrates",
                                          "Lignins","Methyl_ketones","N_compounds","Phenols","Unidentified"))



ggplot(final_new_data,aes(x = meas, y=predicted))+
  geom_point()+
  geom_abline(color="grey",size = 1)+
  theme_bw()+
  facet_wrap(Index ~ ., scales = "free",  ncol= 5)+
  geom_text(aes(x = Inf, y = -Inf, label = paste("Best_fit =", model_name)), vjust = -5, hjust = 1, size = 4)+
  geom_text(aes(x = Inf, y = -Inf, label = paste("Rsq =", Rsquared)), vjust = -3.5, hjust = 1, size = 4)+
  geom_text(aes(x = Inf, y = -Inf, label = paste("RMSD =", RMSD)), vjust = -2, hjust = 1, size = 4)+
  geom_text(aes(x = Inf, y = -Inf, label = paste("AIC =", AIC)), vjust = -0.5, hjust = 1, size = 4)+
  xlab("Reference (%)")+
  ylab("Predicted (%)")





#VIP graph##
spectrum      <- spectraAvg_msc
measured      <- GCMS.dat[,c("ID","Carbohydrates")]


#A subset of the spectra with measure 
Xdata     <- spectrum
XdataSort = Xdata
# Sort measured index data
names = colnames(measured)
names = names[2]
measured = measured[order(measured[[names]]),]

# This subset is sorted in the same way the measured index is ordered 
Measured = measured
names = as.character(measured$ID)
for(s in 1:nrow(measured)){
  nameTmp       <- measured$ID[s]
  row           <- match(nameTmp, Xdata$ID)
  XdataSort[s,] <- Xdata[row,]
}

#Part of this dataset is stored for calibration
calib   = seq(1,117, by=3)

#Create the calibration dataset
calibData <- XdataSort[-c(calib),]

#The rest is kept for validation
validData <- XdataSort[c(calib),]

#In this subset, an additional variable is created in which the measured index are stored
calibData$meas    =  NA
validData$meas    =  NA
totalSpectraCalib =  nrow(calibData)
totalSpectraValid =  nrow(validData)

#a for loop to include the measured index into the dataframe with th spectra for the calibration set
for (i in 1: totalSpectraCalib){
  nameTmp            <- calibData$ID[i]
  row                <- match(nameTmp, Measured$ID)
  calibData$meas[i]  =  measured[["Carbohydrates"]][row]
}
for (i in 1: totalSpectraValid){
  nameTmp            <- validData$ID[i]
  row                <- match(nameTmp, Measured$ID)
  validData$meas[i]  =  Measured[["Carbohydrates"]][row]
}
# PLS model
X          <- as.matrix(calibData$spc) #the spectra
Y          <- as.matrix(calibData$meas) #the measured iron
pls.model  <- mvr(Y ~ X, ncomp = 7, method = "kernelpls") #mvr funtion performs mean-centering automatically (manual)
loadWeight = loading.weights(pls.model)
XScores    = scores(pls.model)
XLoadings  = loadings(pls.model)
YScores    = Yscores(pls.model)
YLoadings  = Yloadings(pls.model)
ExpVar     = explvar(pls.model)
coefficients = coef(pls.model)


plot(pls.model, plottype = "coef", ncomp=1:3, 
     legendpos = "bottomleft",  labels = "numbers", xlab = "nm")      
coefficients = coef(pls.model)
sum.coef = sum(sapply(coefficients, abs))
coefficients = coefficients * 100 / sum.coef
coefficients = sort(coefficients[, 1 , 1])
plot(coefficients)


as.data.frame(coefficients)

barplot(coefficients)
str(as.data.frame(coefficients)
)
View(coefficients)

measured_VIP      <- GCMS.dat[,c("ID","Carbohydrates")]

calibData<- spectrum[XdataSort,]
calibData_Y<- measured[XdataSort,]  

X             <- as.matrix(calibData$spc)
Y             <- as.matrix(calibData_Y$Carbohydrates)
c             <- 10

pls.model     <- mvr(Y~ X, ncomp = 7, method = "kernelpls")



VIPs <- varImp(pls.model)

str(VIPs)
ggplot(VIPs, aes(x = rownames(VIPs), y=Overall))+
  geom_point()