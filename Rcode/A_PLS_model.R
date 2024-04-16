#PLS model#

#Clean##
rm(list= is())
graphics.off()

##Working directory##
Path <- "C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/Mid_Infrared/Data"
setwd(Path)
##Package management##

installpkg <- function(x){
  if(x %in% rownames(installed.packages())==FALSE) {
    if (x %in% rownames(available.packages())== FALSE) {
      paste(x, "is not a valid package - please check again..")
    } else {
      install.packages(x)
    }
  } else {
    paste(x, "package already installed...")
  }
}

#install necessary packages 
required_packages <-c("prospectr","ggplot2","data.table","reshape2","gridExtra","pls","ChemoSpec","caret")
lapply(required_packages, installpkg)

source("Rcode/fun/misc.R")
source("Rcode/fun.R")

library(prospectr)
library(ggplot2)
library(data.table)
library(reshape2)
library(gridExtra)
library(pls)
library(ChemoSpec)
library(caret)
###################################################################################
##Average spectra based on the replicates##
###################################################################################
spectraMatrix <- as.matrix(dataset$spc)
idMatrix <- as.matrix(dataset$ID)



Averages <- SpecAver(spectraMatrix, 4, idMatrix)

spectraAveragelist           <- list()
spectraAveragelist$ID        <- Averages$ID
spectraAveragelist$Site      <- Averages$Site
spectraAveragelist$Treatment <- Averages$Treatment
spectraAveragelist$Depth     <- Averages$Depth
spectraAverage               <- data.frame(spectraAveragelist)
spectraAverage$spc           <- data.frame(Averages$SpecMeans)
colnames(spectraAverage$spc) <- colnames(dataset$spc)

#
rm(spectraMatrix)

#####################################################################################
#GGPlot###
#####################################################################################

WSSL_hy <- new("hyperSpec",spc=as.matrix(spectraAverage$spc),
               wavelength = as.numeric(colnames(spectraAverage$spc)),
               data=spectraAverage[,1:4],
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
  mov_spc                        <- movav(spectraAverage$spc,11)
  tmp                            <- data.frame(mov_spc)
  colnames(tmp)                  <- colnames(mov_spc)
  spectraAvg_movav               <- subset(spectraAverage, select = -spc)
  spectraAvg_movav$spc           <- tmp
  ###################Sample_PLOT##########################
  #windows()
  sample_movav                   <- subset(spectraAvg_movav[1:5,])
  WSSL_hy_movav <- new("hyperSpec",spc=as.matrix(sample_movav$spc),
                 wavelength = as.numeric(colnames(sample_movav$spc)),
                 data=sample_movav[,1:4],
                 label=list(.wavelength="nm",spc="Absorbance"))
  
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
  sg_spc                         <- savitzkyGolay(spectraAverage$spc,p=3,w=21,m=0)
  tmp                            <- data.frame(sg_spc)
  colnames(tmp)                  <- colnames(sg_spc)
  spectraAvg_sg                  <- subset(spectraAverage, select = -spc)
  spectraAvg_sg$spc              <- tmp
  ###################Sample_PLOT##########################
  #windows()
  sample_sg                      <- subset(spectraAvg_sg[1:5,])
  WSSL_hy_sg    <- new("hyperSpec",spc=as.matrix(sample_sg$spc),
                       wavelength = as.numeric(colnames(sample_sg$spc)),
                       data=sample_sg[,1:4],
                       label=list(.wavelength="nm",spc="Absorbance"))
  
  ggplot(as.long.df(WSSL_hy_sg),aes(x=.wavelength,y=spc,color=ID,group=ID))+
    geom_line()+
    theme(legend.position = "top")+
    theme_bw()+
    xlab(bquote('Wavenumber '(cm^-1)))+
    ylab('Absorbance (Sav.-Golay)')+
    scale_x_reverse()
  
  #------------------------------------------------------
  #First and second derivatives
  #------------------------------------------------------
  firstDer_spc                   <- t(diff(t(spectraAverage$spc), differences = 1)) 
  tmp                            <- data.frame(firstDer_spc)
  colnames(tmp)                  <- colnames(firstDer_spc)
  spectraAvg_firstDer            <- subset(spectraAverage, select = -spc)
  spectraAvg_firstDer$spc        <- tmp
  
  secondDer_spc                  <- t(diff(t(spectraAverage$spc), differences = 2)) 
  tmp                            <- data.frame(secondDer_spc)
  colnames(tmp)                  <- colnames(secondDer_spc)
  spectraAvg_secondDer           <- subset(spectraAverage, select = -spc)
  spectraAvg_secondDer$spc       <- tmp
  ###################Sample_PLOT##########################
  #windows()
  sample_firstDer                <- subset(spectraAvg_firstDer[1:5,])
  WSSL_hy_firstDer  <- new("hyperSpec",spc=as.matrix(sample_firstDer$spc),
                       wavelength = as.numeric(colnames(sample_firstDer$spc)),
                       data=sample_firstDer[,1:4],
                       label=list(.wavelength="nm",spc="Absorbance"))
  
  ggplot(as.long.df(WSSL_hy_firstDer),aes(x=.wavelength,y=spc,color=ID,group=ID))+
    geom_line()+
    theme(legend.position = "top")+
    theme_bw()+
    xlab(bquote('Wavenumber '(cm^-1)))+
    ylab('Absorbance (firstDer_spc)')+
    scale_x_reverse()
  ###################Sample_PLOT##########################
  #windows()
  sample_secondDer               <- subset(spectraAvg_secondDer[1:5,])
  WSSL_hy_secondDer  <- new("hyperSpec",spc=as.matrix(sample_secondDer$spc),
                           wavelength = as.numeric(colnames(sample_secondDer$spc)),
                           data=sample_secondDer[,1:4],
                           label=list(.wavelength="nm",spc="Absorbance"))
  
  ggplot(as.long.df(WSSL_hy_secondDer),aes(x=.wavelength,y=spc,color=ID,group=ID))+
    geom_line()+
    theme(legend.position = "top")+
    theme_bw()+
    xlab(bquote('Wavenumber '(cm^-1)))+
    ylab('Absorbance (firstDer_spc)')+
    scale_x_reverse()
  
  #------------------------------------------------------
  #Savitzky-Golay derivatives
  #------------------------------------------------------
  sgSeconder_spc                 <- savitzkyGolay(spectraAverage$spc,p=3,w=21,m=2)
  tmp                            <- data.frame(sgSeconder_spc)
  colnames(tmp)                  <- colnames(sgSeconder_spc)
  spectraAvg_sgSeconder          <- subset(spectraAverage, select = -spc)
  spectraAvg_sgSeconder$spc      <- tmp
  ###################Sample_PLOT##########################
  #windows()
  sample_sgSeconder              <- subset(spectraAvg_sgSeconder[1:5,])
  WSSL_hy_sgSeconder<- new("hyperSpec",spc=as.matrix(sample_sgSeconder$spc),
                           wavelength = as.numeric(colnames(sample_sgSeconder$spc)),
                           data=sample_sgSeconder[,1:4],
                           label=list(.wavelength="nm",spc="Absorbance"))
  
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
  msc_spc                        <- msc(X = as.matrix(spectraAverage$spc))
  tmp                            <- data.frame(msc_spc)
  colnames(tmp)                  <- colnames(msc_spc)
  spectraAvg_msc                 <- subset(spectraAverage, select = -spc)
  spectraAvg_msc$spc             <- tmp
  ###################Sample_PLOT##########################
  #windows()
  sample_msc                     <- subset(spectraAvg_msc[1:5,])
  WSSL_hy_msc   <- new("hyperSpec",spc=as.matrix(sample_msc$spc),
                       wavelength = as.numeric(colnames(sample_msc$spc)),
                       data=sample_msc[,1:4],
                       label=list(.wavelength="nm",spc="Absorbance"))
  
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
  snv_spc                        <- standardNormalVariate(X =spectraAverage$spc)
  tmp                            <- data.frame(snv_spc)
  colnames(tmp)                  <- colnames(snv_spc)
  spectraAvg_snv                 <- subset(spectraAverage, select = -spc)
  spectraAvg_snv$spc             <- tmp
  ###################Sample_PLOT##########################
  #windows()
  sample_snv                    <- subset(spectraAvg_snv[1:5,])
  WSSL_hy_snv   <- new("hyperSpec",spc=as.matrix(sample_snv$spc),
                       wavelength = as.numeric(colnames(sample_snv$spc)),
                       data=sample_snv[,1:4],
                       label=list(.wavelength="nm",spc="Absorbance"))
  
  ggplot(as.long.df(WSSL_hy_snv),aes(x=.wavelength,y=spc,color=ID,group=ID))+
    geom_line()+
    theme(legend.position = "top")+
    theme_bw()+
    xlab(bquote('Wavenumber '(cm^-1)))+
    ylab('Absorbance (SNV)')+
    scale_x_reverse()
  
  #------------------------------------------------------
  # SNV-Detrend (Standard Normal Variate)
  # SNV-Detrend further accoutn for wavelength-dependent scattering effects
  #------------------------------------------------------
  dt_spc                         <- detrend(X = spectraAverage$spc, wav = as.numeric(colnames(spectraAverage$spc)))
  tmp                            <- data.frame(dt_spc)
  colnames(tmp)                  <- colnames(dt_spc)
  spectraAvg_dt                  <- subset(spectraAverage, select = -spc)
  spectraAvg_dt$spc              <- tmp
  ###################Sample_PLOT##########################
  #windows()
  sample_dt                    <- subset(spectraAvg_dt[1:5,])
  WSSL_hy_dt    <- new("hyperSpec",spc=as.matrix(sample_dt$spc),
                       wavelength = as.numeric(colnames(sample_dt$spc)),
                       data=sample_dt[,1:4],
                       label=list(.wavelength="nm",spc="Absorbance"))
  
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
  Xscaled                        <- scale(spectraAvg_sgSeconder$spc, center = T, scale =T)
  tmp                            <- data.frame(Xscaled)
  colnames(tmp)                  <- colnames(Xscaled)
  spectraAvg_scaled              <- subset(spectraAverage, select = -spc)
  spectraAvg_scaled$spc          <- tmp
  ###################Sample_PLOT##########################
  #windows()
  sample_scaled                  <- subset(spectraAvg_scaled[1:5,])
  WSSL_hy_scaled <- new("hyperSpec",spc=as.matrix(sample_scaled$spc),
                       wavelength = as.numeric(colnames(sample_scaled$spc)),
                       data=sample_scaled[,1:4],
                       label=list(.wavelength="nm",spc="Absorbance"))
  
  ggplot(as.long.df(WSSL_hy_scaled),aes(x=.wavelength,y=spc,color=ID,group=ID))+
    geom_line()+
    theme(legend.position = "top")+
    theme_bw()+
    xlab(bquote('Wavenumber '(cm^-1)))+
    ylab('Absorbance (Centering-Scaling)')+
    scale_x_reverse()
  
}



















##########################################################################################
##########################################################################################
#Load data for Partial least squares regression (PLS)
getwd()
GCMS.dat<-read.csv("C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/Mid_Infrared/Data/Proportion_pyrolysis_V2.csv")

#----------------------------------------------------------------------
# Partial Least Squares regression (PLS)
#----------------------------------------------------------------------
str(GCMS.dat)

Compounds_group = c("Aliphatics","Alkanes","Alkenes","Benzenes","Carbohydrates","Lignins",       
                    "Methyl.ketones","N.compounds","Phenols","Unidentified","Guaiacols","Syringols",
                    "Intact_lignin","Sedges","Cellulose")

PLSfactors      = 1:23
nObs            = 27
spectrometer = "MIR"

Nit =1;

for (w in 1: length(Compounds_group)){
  
  #----------------------------------------
  #THe measured Index is loaded
  #---------------------------------------
  Index         <- Compounds_group[w]
  measured      <- GCMS.dat[,c("ID",Index)]
  
  
  
  #----------------------------------------
  #Partial Least Squares regression (PLS)
  #---------------------------------------
  number_of_preproc    <- 9
  number_of_factors    <- length(PLSfactors)
  number_of_total_iterations <- length(Compounds_group) * number_of_preproc * number_of_factors
  
  #creat matrix to write the result
  number_of_iterations <- number_of_preproc*number_of_factors
  reg_result           <- data.frame(data=NA, Rsquared=NA, p=NA, RMSD=NA, RPD=NA, SSE=NA, AIC=NA)
  n                    <- 1
  
  
  for (p in 1:number_of_preproc){   #loop for the difference preprocessing methods
    for(ci in 1:number_of_factors){ #loop for different amout of pls factors
      c = PLSfactors[ci]
      
      #the preprocessed spectrum is selected
      preproc = p
      
      if(preproc==1){
        spectrum <- spectraAvg_movav
        preprocName <-"movav"}
      if(preproc==2){
        spectrum <- spectraAvg_msc
        preprocName <-"msc"}
      if(preproc==3){
        spectrum <- spectraAvg_dt
        preprocName <-"dt"}
      if(preproc==4){
        spectrum <- spectraAvg_sg
        preprocName <-"sg"}
      if(preproc==5){
        spectrum <- spectraAvg_snv
        preprocName <-"snv"}
      if(preproc==6){
        spectrum <- spectraAvg_firstDer
        preprocName <-"firstDer"}
      if(preproc==7){
        spectrum <- spectraAvg_secondDer
        preprocName <-"secondDer"}
      if(preproc==8){
        spectrum <- spectraAvg_sgSeconder
        preprocName <-"sgSeconDer"}
      if(preproc==9){
        spectrum <- spectraAverage
        preprocName <-"no preproc"
      }
      tableStorageName = paste(spectrometer,"_",Index,"_",preprocName,c,sep="")
      
      #Data partitioning
      XdataSort <- createDataPartition(spectrum$Treatment, p = 0.75, list =FALSE)
      calibData <- spectrum[XdataSort,]
      validData <- spectrum[-XdataSort,]
      
      
      #
      calibData$meas = NA
      validData$meas = NA
      totalSpectraCalib = nrow(calibData)
      totalSpectraValid = nrow(validData)
      
      for (i in 1: totalSpectraCalib) {
        nametmp           <- calibData$ID[i]
        row               <- match(nametmp, measured$ID)
        calibData$meas[i] <- measured[[Index]][row]
      }
      
      for (i in 1: totalSpectraValid) {
        nametmp           <- validData$ID[i]
        row               <- match(nametmp, measured$ID)
        validData$meas[i] <- measured[[Index]][row]
      }

      #PLS Model
      X             <- as.matrix(calibData$spc)
      Y             <- as.matrix(calibData$meas)
      pls.model     <- mvr(Y~ X, ncomp = c, method = "kernelpls", validation = "LOO")
      loadWeight    <- loading.weights(pls.model)
      XScores       <- scores(pls.model)
      XLoadings     <- loadings(pls.model)
      YScores       <- Yscores(pls.model)
      YLoadings     <- Yloadings(pls.model)
      Expvar        <- explvar(pls.model)
      
      
      
      #validation value from model created by calibration data
      restSpectra       <- as.matrix(validData$spc)
      pls.pred          <- predict(pls.model, newdata= restSpectra, ncomp=c)
      new_val           <- as.matrix(pls.pred)
      
      
      results           <-subset(validData, select = -spc)
      results$predicted <- NA
      results$predicted <- new_val
      
      #linear regression between measured and estimated
      reg                <- lm(results$predicted ~results$meas)
      Rvalue             <- summary(reg)$r.squared
      pvalue             <- summary(reg)$coefficients[2,4]
      
      #calculate RMSE, RPD and AIC
      temp <- (results$meas -results$predicted)^2
      RMSD <- (mean(temp, na.rm =TRUE))^0.5
      RPD  <- sd(results$meas)/RMSD
      SSE  <- sum((results$meas-results$predicted)^2, na.rm = TRUE)
      AIC  <- nObs*log(SSE/nObs)+2*c
      
      #the R and p-value for each regression combination are stored
      reg_result[n,1] <- tableStorageName
      reg_result[n,2] <- Rvalue
      reg_result[n,3] <- pvalue
      reg_result[n,4] <- RMSD
      reg_result[n,5] <- RPD
      reg_result[n,6] <- SSE
      reg_result[n,7] <- AIC
      
      n = n+1
      Nit = Nit +1
      
      Name <- paste(Index, format(nObs,digits =2), sep = "")
      dir.create(paste("PLS_model","/",Name, sep = ""), showWarnings = FALSE)
      tableName = paste(Path,"/PLS_model/", Name,"/", tableStorageName,".txt",sep = "")
      write.table(results, file = tableName, row.names = FALSE)
      
      Percent = paste(format((Nit/number_of_total_iterations)*100,digits =4), "%", sep ="")
      print(paste(tableStorageName,"(",Percent,"completed)",sep = ""))
      }
  }
  Name <- paste(Index, format(nObs,digits =2), sep = "")
  tableName = paste(Path,"/PLS_model/", Name,"/", spectrometer,Index,"diffFactors.txt",sep = "")
  write.table(reg_result, file = tableName, row.names = FALSE)
}

reg_result[min(reg_result$RMSD),]


fit_model.selection <- read.table("PLS_model/Aliphatics27/MIRAliphaticsdiffFactors.txt", header = T)
fit_model.selection <- read.table("PLS_model/Carbohydrates27/MIRCarbohydratesdiffFactors.txt", header = T)
fit_model.selection <- read.table("PLS_model/Alkenes27/MIRAlkenesdiffFactors.txt", header = T)
fit_model.selection <- read.table("PLS_model/Benzenes27/MIRBenzenesdiffFactors.txt", header = T)
fit_model.selection <- read.table("PLS_model/Lignins27/MIRLigninsdiffFactors.txt", header = T)
fit_model.selection <- read.table("PLS_model/Methyl.ketones27/MIRMethyl.ketonesdiffFactors.txt", header = T)
fit_model.selection <- read.table("PLS_model/N.compounds27/MIRN.compoundsdiffFactors.txt", header = T)
fit_model.selection <- read.table("PLS_model/Phenols27/MIRPhenolsdiffFactors.txt", header = T)
fit_model.selection <- read.table("PLS_model/Unidentified27/MIRUnidentifieddiffFactors.txt", header = T)
fit_model.selection <- read.table("PLS_model/Intact_lignin27/MIRIntact_lignindiffFactors.txt", header = T)
fit_model.selection <- read.table("PLS_model/Sedges27/MIRSedgesdiffFactors.txt", header = T)



fit_model.selection[fit_model.selection$AIC == sort(fit_model.selection$AIC)[1],]

  
plot.dat <- read.table("PLS_model/Phenols27/", header = T)


library(ggplot2)
ggplot(plot.dat, aes(x = meas, y=predicted))+
  geom_point()+
  geom_abline()



#VIP graph##
spectrum <- spectraAvg_msc
XdataSort <- createDataPartition(spectrum$Treatment, p = 0.75, list =FALSE)
measured      <- GCMS.dat[,c("ID","Carbohydrates")]

calibData<- spectrum[XdataSort,]
calibData_Y<- measured[XdataSort,]  

X             <- as.matrix(calibData$spc)
Y             <- as.matrix(calibData_Y$Carbohydrates)
c             <- 10

pls.model     <- mvr(Y~ X, ncomp = c, method = "kernelpls", validation = "LOO")

VIPs <- varImp(pls.model)

str(VIPs)
ggplot(VIPs, aes(x = rownames(VIPs), y=Overall))+
  geom_point()
