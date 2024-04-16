
# The workspace is cleared
rm(list=ls()) # The workspace is cleared
graphics.off()
#initial variable
getwd()
setwd("C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/Mid_Infrared/")
Sys.setenv(LANGUAGE="EN")

#source("C_Areabased_Normalization.R")
load("Data/Rawdata/MIRs.RData")
#source("fun/misc.R")
source("Rcode/fun/C_Spectraaverage.R")
#install.packages("ChemoSpec")
# source the file containing other functions
library(stats)
library(dplyr)
library(grDevices)
library(ggpubr)
library(prospectr)
library(reshape2)
library(gridExtra)
library(hyperSpec)
library(ChemoSpec)
library(ChemoSpecUtils)
library(ggrepel)
library(factoextra)
library(mvoutlier)
# install.packages("plotly")
library(plotly)

# The workspace is cleared
rm(list=ls()) # The workspace is cleared
graphics.off()
#initial variable
getwd()
setwd("C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/Mid_Infrared/")
Sys.setenv(LANGUAGE="EN")

#source("C_Areabased_Normalization.R")
load("Data/Rawdata/MIRs.RData")
#source("fun/misc.R")
source("Rcode/fun/C_Spectraaverage.R")
#install.packages("ChemoSpec")
# source the file containing other functions
library(stats)
library(dplyr)
library(grDevices)
library(ggpubr)
library(prospectr)
library(reshape2)
library(gridExtra)
library(hyperSpec)
library(ChemoSpec)
library(ChemoSpecUtils)
library(ggrepel)
library(factoextra)
library(mvoutlier)
library(caret)
#install.packages("doParallel")
library(doParallel)
library(tidyr)

#------------------------------------------------------------------
#Calibratio average spectra are caluacated based on the replicates
#------------------------------------------------------------------
spectraMatrix=as.matrix(dataset$spc)
columns <- colnames(spectraMatrix) %in% c("4000":"650")
subsetMatrix <- spectraMatrix[, columns]
idMatrix=as.matrix(dataset$ID)

Averages <- SpecAver(subsetMatrix,4,idMatrix, thr = .07)

#View(Averages$Check)
## I need to build a loop for check data
# Averages$Check
# 
# Check = as.data.frame(Averages$Check)
# 
# View(Check)
# for (i in seq(from =1, to = 117)){
#   
# }

# ##############################################################
# 
# WSSL_hy <- new("hyperSpec",spc=as.matrix(outlier $spc),
#                wavelength = as.numeric(colnames(outlier $spc)),
#                data=outlier [,1:6],
#                label=list(.wavelength="nm",spc="Absorbance"))
# 
# WSSL_hy
# ###Plotting for whole data##
# ggplot(as.long.df(WSSL_hy),aes(x=.wavelength,y=spc,color=ID,group=ID))+
#   geom_line()+
#   facet_grid(Depth~Treatment)+
#   theme(legend.position = "top")+
#   theme_bw()+
#   xlab(bquote('Wavelength '(cm^-1)))+
#   ylab('Absorbance')+
#   scale_x_reverse()
# 
# 

#A list containing the average spectra's and ID's is constructed
spectraAverageslist           <- list() 
spectraAverageslist$ID        <- Averages$ID                      # The ID of the spectra
spectraAverageslist$Site      <- Averages$Site                    # The sitename
spectraAverageslist$Treatment    <- Averages$Treatment                  # The replica number of the soil core
spectraAverageslist$Depth     <- Averages$Depth                   # The upper depth of the sample
spectraAverages               <- data.frame(spectraAverageslist)  # The list is converted to a data frame
tmp                           <- data.frame(Averages$SpecMeans)   # The average spectra are copied
spectraAverages$spc           <- tmp                              # The average spectra are added
colnames(spectraAverages$spc) <- colnames(subsetMatrix) # The wavenumbers are used as column names for the spectra

rm(spectraMatrix,tmp)



dataset2<-spectraAverages
# dataset2$depth=factor(dataset2$depth,levels = c("T","M","B"),
#                       labels = c("0-5 cm","15-20 cm","45-50 cm"))
# 
# dataset2$treatment=factor(dataset2$treatment,levels = c("U","R","D"),
#                           labels = c("Undrained","Rewetted","Drained"))


EA.dat=read.csv("Data/EA_data.csv")
#Hodgkins.dat=read.csv("Data/Outcomes_Hodgkins.csv")
Microbial.dat=read.csv("Data/Microbial_data.csv")


#Hodgkins.dat <- subset(Hodgkins.dat, select = -c(X))
#Variables.dat<-left_join(EA.dat, Hodgkins.dat)
Variables.dat<-left_join(EA.dat, Microbial.dat)


dataset4 <-left_join(Variables.dat,dataset2)


dataset4$Depth=factor(dataset4$Depth,levels = c("T","M","B"),
                      labels = c("0-5 cm","15-20 cm","45-50 cm"))

dataset4$Treatment=factor(dataset4$Treatment,levels = c("U","R","D"),
                          labels = c("Undrained","Rewetted","Drained"))

dataset4$TsR=factor(dataset4$TsR,levels = c("L","M","S"),
                          labels = c(">25","10-25","<10"))


WSSL_hy <- new("hyperSpec",spc=as.matrix(dataset4$spc),
               wavelength = as.numeric(colnames(dataset4$spc)),
               data=dataset4[,1:11],
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



########################################################################
#Noise removal######
#########################################################################

dataset4$spc2 <- savitzkyGolay(dataset4$spc,p=3,w=21,m=0)
source("Rcode/Reference/C_Areabased_Normalization.R")
dataset4$spc2<-Normalization100(dataset4$spc2)

#dataset4$spc2 <- binning(dataset4$spc2, bin.size = 10)

WSSL_hy <- new("hyperSpec",spc=as.matrix(dataset4$spc2),
               wavelength = as.numeric(colnames(dataset4$spc2)),
               data=dataset4[,1:11],
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



#------------------------------------------------------------------
#Calibratio average spectra are caluacated based on the replicates
#------------------------------------------------------------------
spectraMatrix=as.matrix(dataset$spc)
columns <- colnames(spectraMatrix) %in% c("4000":"650")
subsetMatrix <- spectraMatrix[, columns]
idMatrix=as.matrix(dataset$ID)

Averages <- SpecAver(subsetMatrix,4,idMatrix, thr = .07)

#View(Averages$Check)
## I need to build a loop for check data
# Averages$Check
# 
# Check = as.data.frame(Averages$Check)
# 
# View(Check)
# for (i in seq(from =1, to = 117)){
#   
# }

# ##############################################################
# 
# WSSL_hy <- new("hyperSpec",spc=as.matrix(outlier $spc),
#                wavelength = as.numeric(colnames(outlier $spc)),
#                data=outlier [,1:6],
#                label=list(.wavelength="nm",spc="Absorbance"))
# 
# WSSL_hy
# ###Plotting for whole data##
# ggplot(as.long.df(WSSL_hy),aes(x=.wavelength,y=spc,color=ID,group=ID))+
#   geom_line()+
#   facet_grid(Depth~Treatment)+
#   theme(legend.position = "top")+
#   theme_bw()+
#   xlab(bquote('Wavelength '(cm^-1)))+
#   ylab('Absorbance')+
#   scale_x_reverse()
# 
# 

#A list containing the average spectra's and ID's is constructed
spectraAverageslist           <- list() 
spectraAverageslist$ID        <- Averages$ID                      # The ID of the spectra
spectraAverageslist$Site      <- Averages$Site                    # The sitename
spectraAverageslist$Treatment    <- Averages$Treatment                  # The replica number of the soil core
spectraAverageslist$Depth     <- Averages$Depth                   # The upper depth of the sample
spectraAverages               <- data.frame(spectraAverageslist)  # The list is converted to a data frame
tmp                           <- data.frame(Averages$SpecMeans)   # The average spectra are copied
spectraAverages$spc           <- tmp                              # The average spectra are added
colnames(spectraAverages$spc) <- colnames(subsetMatrix) # The wavenumbers are used as column names for the spectra

rm(spectraMatrix,tmp)



dataset2<-spectraAverages
# dataset2$depth=factor(dataset2$depth,levels = c("T","M","B"),
#                       labels = c("0-5 cm","15-20 cm","45-50 cm"))
# 
# dataset2$treatment=factor(dataset2$treatment,levels = c("U","R","D"),
#                           labels = c("Undrained","Rewetted","Drained"))


EA.dat=read.csv("Data/EA_data.csv")
#Hodgkins.dat=read.csv("Data/Outcomes_Hodgkins.csv")
Microbial.dat=read.csv("Data/Microbial_data.csv")


#Hodgkins.dat <- subset(Hodgkins.dat, select = -c(X))
#Variables.dat<-left_join(EA.dat, Hodgkins.dat)
Variables.dat<-left_join(EA.dat, Microbial.dat)


dataset4 <-left_join(Variables.dat,dataset2)


dataset4$Depth=factor(dataset4$Depth,levels = c("T","M","B"),
                      labels = c("0-5 cm","15-20 cm","45-50 cm"))

dataset4$Treatment=factor(dataset4$Treatment,levels = c("U","R","D"),
                          labels = c("Undrained","Rewetted","Drained"))

dataset4$TsR=factor(dataset4$TsR,levels = c("L","M","S"),
                    labels = c(">25","10-25","<10"))


WSSL_hy <- new("hyperSpec",spc=as.matrix(dataset4$spc),
               wavelength = as.numeric(colnames(dataset4$spc)),
               data=dataset4[,1:11],
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



########################################################################
#Noise removal######
#########################################################################

dataset4$spc2 <- savitzkyGolay(dataset4$spc,p=3,w=21,m=0)
source("Rcode/Reference/C_Areabased_Normalization.R")
dataset4$spc2<-Normalization100(dataset4$spc2)

#dataset4$spc2 <- binning(dataset4$spc2, bin.size = 10)

WSSL_hy <- new("hyperSpec",spc=as.matrix(dataset4$spc2),
               wavelength = as.numeric(colnames(dataset4$spc2)),
               data=dataset4[,1:11],
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
##PLSDA##
###################################################################################
X <- dataset4[,c(1,3,4,5)]
Y <- as.data.frame(dataset4$spc2)
Y$ID <- dataset4$ID

merged_df <- merge(X, Y, suffixes = c('', ''), by= "ID")
str(merged_df)
merged_df <- merged_df[,-c(1,4)] 
#merged_df <- merged_df[merged_df$Treatment %in% c("Drained","Undrained")]
#merged_df$Treatment <- droplevels(merged_df$Treatment)

# #Delete all factors
# indexTreatment <- createDataPartition(merged_df$Treatment, p= 0.7, list = FALSE)
# trainData      <- merged_df[indexTreatment,]
# testData       <- merged_df[-indexTreatment,]


set.seed(123)

control1 <- caret::trainControl(method = "LOOCV", number = 5, classProbs =TRUE,verboseIter = TRUE)
#cl <- makePSOCKcluster(5)
#registerDoParallel(cl)                                

model_asRecieved <- caret::train(Treatment~.,
                                 data=merged_df[,-1],
                                 method ="pls",
                                 trControl =control1,
                                 preProcess= "center",
                                 tuneLength =20,
                                 summaryFunction =multiClassSummary)

#stopCluster(cl)

confusionMatrix(predict(model_asRecieved,testData),as.factor(testData$Treatment))

plot(model_asRecieved)
model_tune <-data.frame(cbind(1:length(model_asRecieved$results$Accuracy),model_asRecieved$results$Accuracy)) %>% 
  rename(Components=1, Accuracy =2)


  
model_tune %>% ggplot(aes(Components,Accuracy))+
  geom_line(color="steelblue")+
  geom_point(aes(Components,Accuracy), color="steelblue", size=2, alpha=0.5)+
  geom_point(aes(y=max(Accuracy),x=Components[which.max(Accuracy)]), color="red", size=3)+
  ggtitle("Cross-validation tuning", subtitle = "Best accuracy marked with red")+
  scale_x_continuous(breaks=seq(1,25,2))

design <- data.frame(sample=merged_df$Site)

plsda_model <- mixOmics::plsda(merged_df[,-c(1,2)],
                              merged_df$Treatment, 
                              multilevel =  design,
                              ncomp = 17)


#####################################################################################
##PLSDA## Treatment 
###################################################################################
X <- dataset4[,c(1,3,4,5)]
Y <- as.data.frame(dataset4$spc2)
Y$ID <- dataset4$ID

merged_df <- merge(X, Y, suffixes = c('', ''), by= "ID")
str(merged_df)
merged_df <- merged_df[,-c(1,4)] 


plsda_model <- mixOmics::plsda(merged_df[,-c(1,2)],
                               merged_df$Treatment, 
                               #multilevel =  design,
                               ncomp = 20)

#undergo performance evaluation in oder to tune the number of componets to use

perf.plsda_model <- mixOmics::perf(plsda_model, validation ="Mfold", 
                         nrepeat =100, folds=10,
                         progressBar = TRUE, auc =TRUE)

plot(perf.plsda_model,col= mixOmics::color.mixo(5:7),sd =TRUE,
     legend.position = "horizontal")

perf.plsda_model$choice.ncomp
design <- data.frame(sample=merged_df$Site)


plsda_model <- mixOmics::plsda(merged_df[,-c(1,2)],
                               merged_df$Treatment, 
                               multilevel =  design,
                               ncomp = 17)


mixOmics::plotIndiv(plsda_model, group = merged_df$Treatment, 
          ind.names = merged_df$Treatment,
          legend = TRUE, legend.title = 'Treatment',
          title = 'Sample Plot of sPLS-DA on Vac18 data')

score_plot<- as.data.frame(plsda_model$variates$X[,1:2])
score_plot<- as.data.frame(plsda_model$variates$X[,c("comp1","comp2")])
colnames(score_plot) <- c("PC1","PC2")


score_plot$Treatment <-  merged_df$Treatment
score_plot$Site <-  merged_df$Site

PLSDA_Treatment<-ggplot(score_plot, aes(PC1,PC2 , color =Treatment, shape = Treatment))+
                      geom_point(size =2)+
                      scale_color_manual(values=c("#00AFBB","#E7B800","#FC4E07"))+
                      theme_bw()+
                      stat_ellipse(type = "t")+
                      xlab(paste0("Latent variable 1 (45.4 %)"))+
                      ylab(paste0("Latent variable 2 (6.0 %)"))+
                      theme(legend.position = "top",
                            legend.title=element_blank())#+

plsda_model$prop_expl_var


loadings_data <- as_tibble(plsda_model$loadings$X) %>% 
  mutate(INDEX = row.names(plsda_model$loadings$X)) %>% 
  relocate(INDEX)
loadings_data$INDEX <- as.numeric(loadings_data$INDEX)


loading.spc <- as_tibble(loadings_data[,1:3]) %>% 
  gather(comp1:comp2,
         key = PC_axis,
         value = Loadings)

loading.spc$INDEX <- as.numeric(loading.spc$INDEX)
loading.spc$PC_axis <- factor(loading.spc$PC_axis, levels = c("comp1","comp2"),
                              labels = c("Latent variable 1","Latent variable 2"))

Loading1=ggplot(loading.spc,aes(x=INDEX,y=Loadings,group=PC_axis,colour=PC_axis))+ 
  geom_line(size=0.8)+
  theme_bw()+
  #stat_peaks(geom="text",span=11,color="red",parse=TRUE, position = position_nudge(y = 0.01))+
  #stat_valleys(geom="text",span=11,color="red",parse=TRUE, position = position_nudge(y = -0.01))+
  theme(legend.position = "top",
        legend.title =element_blank(),
        panel.grid.major.x = element_line(color = "black",
                                          size = 0.1,
                                          linetype = 2),
        panel.grid.major.y = element_line(color = "grey",
                                          size = 0.1,
                                          linetype = 2))+
  #xlab(bquote('Wavenumber '(cm^-1)))+
  ylab('Loadings')+
  scale_x_reverse(n.breaks = 10)+
  grids(linetype = "dashed")+
  scale_color_manual(values=c("#0073C2FF", "#EFC000FF"))
Loading1



q <- ggplotly(Loading1)

q <- q  %>% layout(dragmode = "pan")
q



###########################################################################################3
#PLSDA_Depth
############################################################################################

X <- dataset4[,c(1,3,4,5)]
Y <- as.data.frame(dataset4$spc2)
Y$ID <- dataset4$ID

merged_df <- merge(X, Y, suffixes = c('', ''), by= "ID")
str(merged_df)
merged_df <- merged_df[,-c(1,3)] 


plsda_model <- mixOmics::plsda(merged_df[,-c(1,2)],
                               merged_df$Depth, 
                               #multilevel =  design,
                               ncomp = 20)

#undergo performance evaluation in oder to tune the number of componets to use

perf.plsda_model <- mixOmics::perf(plsda_model, validation ="Mfold", 
                                   nrepeat =100, folds=10,
                                   progressBar = TRUE, auc =TRUE)

plot(perf.plsda_model,col= mixOmics::color.mixo(5:7),sd =TRUE,
     legend.position = "horizontal")

perf.plsda_model$choice.ncomp
design <- data.frame(sample=merged_df$Site)


plsda_model <- mixOmics::plsda(merged_df[,-c(1,2)],
                               merged_df$Depth, 
                               multilevel =  design,
                               ncomp = 12)


mixOmics::plotIndiv(plsda_model, group = merged_df$Depth, 
                    ind.names = merged_df$Depth,
                    legend = TRUE, legend.title = 'Treatment',
                    title = 'Sample Plot of sPLS-DA on Depth data')


score_plot<- as.data.frame(plsda_model$variates$X[,1:2])
score_plot<- as.data.frame(plsda_model$variates$X[,c("comp1","comp2")])

colnames(score_plot) <- c("PC1","PC2")


score_plot$Depth <-  merged_df$Depth
score_plot$Site <-  merged_df$Depth


plsda_model$prop_expl_var

PLSDA_Depth<-ggplot(score_plot, aes(PC1,PC2 , color =Depth, shape = Depth))+
                geom_point(size =2)+
                ggsci::scale_color_jco()+
                #scale_color_manual(values=c("#00AFBB","#E7B800","#FC4E07"))+
                theme_bw()+
                stat_ellipse(type = "t")+
                xlab(paste0("Latent variable 1 (45.6 %)"))+
                ylab(paste0("Latent variable 2 (20.3 %)"))+
                theme(legend.position = "top",
                      legend.title=element_blank())#+

plsda_model$loadings


loadings_data <- as_tibble(plsda_model$loadings$X) %>% 
  mutate(INDEX = row.names(plsda_model$loadings$X)) %>% 
  relocate(INDEX)
loadings_data$INDEX <- as.numeric(loadings_data$INDEX)


loading.spc <- as_tibble(loadings_data[,1:3]) %>% 
  gather(comp1:comp2,
         key = PC_axis,
         value = Loadings)

loading.spc$INDEX <- as.numeric(loading.spc$INDEX)
loading.spc$PC_axis <- factor(loading.spc$PC_axis, levels = c("comp1","comp2"),
                              labels = c("Latent variable 1","Latent variable 2"))


Loading2=ggplot(loading.spc,aes(x=INDEX,y=Loadings,group=PC_axis,colour=PC_axis))+ 
  geom_line(size=0.8)+
  theme_bw()+
  #stat_peaks(geom="text",span=11,color="red",parse=TRUE, position = position_nudge(y = 0.01))+
  #stat_valleys(geom="text",span=11,color="red",parse=TRUE, position = position_nudge(y = -0.01))+
  theme(legend.position = "top",
        legend.title =element_blank(),
        panel.grid.major.x = element_line(color = "black",
                                          size = 0.1,
                                          linetype = 2),
        panel.grid.major.y = element_line(color = "grey",
                                          size = 0.1,
                                          linetype = 2))+
  #xlab(bquote('Wavenumber '(cm^-1)))+
  ylab('Loadings')+
  scale_x_reverse(n.breaks = 10)+
  grids(linetype = "dashed")+
  scale_color_manual(values=c("#0073C2FF", "#EFC000FF"))
Loading2



ggarrange(ggarrange(PLSDA_Treatment, PLSDA_Depth, ncol = 2, labels = c("A", "B")),
          Loading1,Loading2,nrow = 3,labels = c("","C","D")) 

ggarrange(ggarrange(PLSDA_Treatment, Loading1, ncol = 2, labels = c("A", "B"), widths = c(1,1.5)),
          ggarrange(PLSDA_Depth, Loading2, ncol = 2, labels = c("C", "D"), widths = c(1,1.5)),nrow = 2) 

?ggarrange()

q <- ggplotly(Loading1)
q <- ggplotly(Loading2)

q <- q  %>% layout(dragmode = "pan")
q

