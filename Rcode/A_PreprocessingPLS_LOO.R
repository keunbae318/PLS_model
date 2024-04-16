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


