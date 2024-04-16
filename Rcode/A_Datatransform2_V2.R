
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
library(plotly)

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



ggplotly(ggplot(as.long.df(WSSL_hy),aes(x=.wavelength,y=spc,color=Site,group=ID))+
           geom_line()+
           facet_grid(Depth~Treatment)+
           theme(legend.position = "top")+
           theme_bw()+
           #xlab(bquote('Wavenumber '(cm^-1)))+
           ylab('Absorbance')+
           scale_x_reverse())

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

# # Robust PCA based on Projection Pursuit (outlier removal)
# res2 <- pcout(x=dataset4$spc2,makeplot = T)   #inspection of distance                  
# outliers<-as.character(dataset4$ID[res2$x.dist1>50])
# outliers
# dataset4 <- subset(dataset4,!(ID %in% outliers))















#Moving average
#dataset4$spc2 <- movav(dataset4$spc,11)
#dataset4$spc2 <-dataset4$spc2*100
#dataset4$spc2<- 2-log10(dataset4$spc)

#Baseline Correction
# wav <- as.numeric(colnames(dataset4$spc))
# dataset4$spc2 <- baseline(dataset4$spc, wav)

#dataset4$spc2<- -log10(dataset4$spc2)



dataset4$spc2 <- savitzkyGolay(dataset4$spc,p=3,w=7,m=0)
#area-based normalization#
dataset4$spc2 <- Normalization100(dataset4$spc2)

#dataset4$spc2  <- standardNormalVariate(X = dataset4$spc2)

#dataset4$spc2<- binning(dataset4$spc2, bin.size = 10) 


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
  scale_x_reverse()#+
  xlim(1600,700)#+
  # geom_vline(xintercept = 1050, linetype="dotted", 
  #            color = "black", size=0.5)+
  # geom_vline(xintercept = 815, linetype="dotted", 
  #            color = "black", size=0.5)+
  # geom_vline(xintercept = 780, linetype="dotted", 
  #            color = "black", size=0.5)+
  # geom_vline(xintercept = 1030, linetype="dotted", 
  #            color = "black", size=0.5)+
  # geom_vline(xintercept = 1510, linetype="dotted", 
  #            color = "black", size=0.5)+
  # geom_vline(xintercept = 1630, linetype="dotted", 
  #            color = "black", size=0.5)

ggplot(as.long.df(WSSL_hy),aes(x=.wavelength,y=spc,color=Site,group=ID))+
  geom_line()+
  facet_grid(Depth~Treatment)+
  theme(legend.position = "top")+
  theme_bw()+
  xlab(bquote('Wavenumber '(cm^-1)))+
  ylab('Absorbance')+
  scale_x_reverse()#+
  #xlim(700,2000)+
  geom_vline(xintercept = 790, linetype="dotted", 
             color = "red", size=1)+
  geom_vline(xintercept = 1050, linetype="dotted", 
             color = "blue", size=1)+
  geom_vline(xintercept = 1510, linetype="dotted", 
             color = "green", size=1)+
  geom_vline(xintercept = 1630, linetype="dotted", 
             color = "green", size=1)+
  geom_vline(xintercept = 1720, linetype="dotted", 
             color = "black", size=1)#+
  geom_vline(xintercept = 2850, linetype="dotted", 
             color = "purple", size=1)+
  geom_vline(xintercept = 2920, linetype="dotted", 
             color = "purple", size=1)
  




########################################################################
#Binning######and continuumRemoval
#########################################################################
# dataset4$spc3<- binning(dataset4$spc2, bin.size = 10)
# 
# dataset4$spc3 <- continuumRemoval(X = dataset4$spc2, 
#                                   wav = as.numeric(colnames(dataset4$spc2)),
#                                   type= "A")
# 
# 
# WSSL_hy <- new("hyperSpec",spc=as.matrix(dataset4$spc3),
#                wavelength = as.numeric(colnames(dataset4$spc3)),
#                data=dataset4[,1:11],
#                label=list(.wavelength="nm",spc="Absorbance"))
# 
# WSSL_hy
# 
# ###Plotting for whole data##
# ggplot(as.long.df(WSSL_hy),aes(x=.wavelength,y=spc,color=Site,group=ID))+
#   geom_line()+
#   facet_grid(Depth~Treatment)+
#   theme(legend.position = "top")+
#   theme_bw()+
#   xlab(bquote('Wavenumber '(cm^-1)))+
#   ylab('Absorbance')+
#   scale_x_reverse()
# 
# 

########################################################################
#Savitzky-Golay######
#########################################################################

#Option 1
#dataset4$spc4 <- savitzkyGolay(dataset4$spc3,p=3,w=21,m=0)

#Option 2

# 
# 
# WSSL_hy <- new("hyperSpec",spc=as.matrix(dataset4$spc4),
#                wavelength = as.numeric(colnames(dataset4$spc4)),
#                data=dataset4[,1:11],
#                label=list(.wavelength="nm",spc="Absorbance"))
# 
# WSSL_hy
# 
# ###Plotting for whole data##
# ggplot(as.long.df(WSSL_hy),aes(x=.wavelength,y=spc,color=Site,group=ID))+
#   geom_line()+
#   facet_grid(Depth~Treatment)+
#   theme(legend.position = "top")+
#   theme_bw()+
#   xlab(bquote('Wavenumber '(cm^-1)))+
#   ylab('Absorbance')+
#   scale_x_reverse()



########################################################################
#Center and auto-scaled matric######
#########################################################################

# dataset4$spc5 <- scale(dataset4$spc4,center = T, scale = F)
# dataset4$spc5 <- scale(dataset4$spc2,center = T, scale = F)
# 
# dataset4$spc5<- binning(dataset4$spc2, bin.size = 10) 
# # dataset4$spc5 <- scale(dataset4$spc4, center = apply(dataset4$spc4,2, median),
# #                        scale = F )
# 
# WSSL_hy <- new("hyperSpec",spc=as.matrix(dataset4$spc5),
#                wavelength = as.numeric(colnames(dataset4$spc5)),
#                data=dataset4[,1:11],
#                label=list(.wavelength="nm",spc="Absorbance"))
# 
# WSSL_hy
# 
# ###Plotting for whole data##
# ggplot(as.long.df(WSSL_hy),aes(x=.wavelength,y=spc,color=Site,group=ID))+
#   geom_line()+
#   facet_grid(Depth~Treatment)+
#   theme(legend.position = "top")+
#   theme_bw()+
#   xlab(bquote('Wavenumber '(cm^-1)))+
#   ylab('Absorbance')+
#   scale_x_reverse()
# 

# row.names(loading_spc )
# loadings_spc
# View(loadings)
# loadings$ID <- row.names(loadings)
# 
# loadings2 <- loadings[,-which(names(loadings) %in% c("ID"))]

#   
#   
# summary(res.pca)
# ?apply
# 
# t(loadings)
# ?ggplot
# 
# View(loadings)
# 
# svDecomp <- svd(x=dataset5$spc5)
# length(svDecomp$d)
# dim(svDecomp$u)
# dim(svDecomp$v)
# svDecomp$v
# pcLoadings <-svDecomp$v
# pcLoadings[1,]
# diagM <- diag(svDecomp$d)
# pcScores <-svDecomp$u %*% diagM
# pairs(pcScores[,1:3])
# colnames(pcScores) <- paste("PC",1:ncol(pcScores),sep="")
# rownames(pcLoadings) <-paste("PC",1:nrow(pcLoadings),sep="")
# 
# 
# ##plot####
# matplot(x=as.numeric(colnames(dataset4$spc5)), y=pcLoadings[,1],type = "l",
#         xlab=bquote('Wavelength '(cm^-1)),ylab="Loading value",
#         lty=1)
# 
# fviz_eig(res.pca)
