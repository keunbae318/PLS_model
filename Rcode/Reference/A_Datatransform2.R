# The workspace is cleared
#initial variable
getwd()
setwd("C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/MIR/Rcode/")
Sys.setenv(LANGUAGE="EN")

#source("fun/misc.R")
#source("fun/SpectraAverage")
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
#---------------------------------------------------------
#fucntion to check whether or not packages are installed
#---------------------------------------------------------
installpkg <- function (x) {
  if(x %in% rownames(installed.packages())==FALSE) {
    if(x %in% rownames(available.packages())==FALSE) {
      paste(x,"is not a valid package -please check again...")
    } else {
      install.packages(x)
    } 
    
  } else {
    paste(x,"package already installed...")
  }
}

# install necessary packages

required_packages <- c("prospectr", "ggplot2", "data.table", "reshpae2", "gridExtra", "pls", "ChemoSpec")
lapply(required_packages, installpkg)


#------------------------------------------------------------------
#Calibratio average spectra are caluacated based on the replicates
#------------------------------------------------------------------
spectraMatrix=as.matrix(dataset$spc)
idMatrix=as.matrix(dataset$ID)

Averages <- SpecAver(spectraMatrix,4,idMatrix, thr = .07)

## I need to build a loop for check data
# Averages$Check
# 
# Check = as.data.frame(Averages$Check)
# 
# View(Check)
# for (i in seq(from =1, to = 117)){
#   
# }
# 91*4+(4*9)
# outlier <- dataset[364:396,]
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



#A list containing the average spectra's and ID's is constructed
spectraAverageslist           <- list() 
spectraAverageslist$ID        <- Averages$ID                      # The ID of the spectra
spectraAverageslist$Site      <- Averages$Site                    # The sitename
spectraAverageslist$Treatment    <- Averages$Treatment                  # The replica number of the soil core
spectraAverageslist$Depth     <- Averages$Depth                   # The upper depth of the sample
spectraAverages               <- data.frame(spectraAverageslist)  # The list is converted to a data frame
tmp                           <- data.frame(Averages$SpecMeans)   # The average spectra are copied
spectraAverages$spc           <- tmp                              # The average spectra are added
colnames(spectraAverages$spc) <- colnames(dataset$spc) # The wavenumbers are used as column names for the spectra

rm(spectraMatrix,tmp)


dataset2<-spectraAverages
# dataset2$depth=factor(dataset2$depth,levels = c("T","M","B"),
#                       labels = c("0-5 cm","15-20 cm","45-50 cm"))
# 
# dataset2$treatment=factor(dataset2$treatment,levels = c("U","R","D"),
#                           labels = c("Undrained","Rewetted","Drained"))

getwd()
EA.dat=read.csv("C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/MIR/Data/EA_data.csv")
str(EA.dat)
str(dataset2)
dataset3<-left_join(dataset2,EA.dat)

dataset4<- dataset3[,c(1,6,2,3,4,7,8,9,10,11,12,5)]

dataset4$Depth=factor(dataset4$Depth,levels = c("T","M","B"),
                      labels = c("0-5 cm","15-20 cm","45-50 cm"))

dataset4$Treatment=factor(dataset4$Treatment,levels = c("U","R","D"),
                          labels = c("Undrained","Rewetted","Drained"))


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

dataset4$spc2 <- movav(dataset4$spc,11)

dataset4$spc2<- msc(X = dataset4$spc2, ref_spectrum = colMeans(dataset4$spc2))


wav <- as.numeric(colnames(dataset4$spc2))
dataset4$spc2 <- baseline(dataset4$spc2, wav)

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

#######################################################################
Binning######and continuumRemoval
########################################################################
#dataset4$spc3<- binning(dataset4$spc2, bin.size = 10)

dataset4$spc3 <- continuumRemoval(X = dataset4$spc2,
                                  wav = as.numeric(colnames(dataset4$spc2)),
                                  type= "A")

WSSL_hy <- new("hyperSpec",spc=as.matrix(dataset4$spc3),
               wavelength = as.numeric(colnames(dataset4$spc3)),
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
#Savitzky-Golay######
#########################################################################
#Option 0
dataset4$spc4 <- savitzkyGolay(dataset4$spc2,p=2,w=11,m=0)

# #Option 1
# dataset4$spc4 <- savitzkyGolay(dataset4$spc2,p=3,w=21,m=0)

#Option 2
#dataset4$spc4 <- savitzkyGolay(dataset4$spc2,p=3,w=41,m=1)

dataset4$spc5  <- standardNormalVariate(X = dataset4$spc4)


wav <- as.numeric(colnames(dataset4$spc5))
dataset4$spc5 <- baseline(dataset4$spc5, wav)

dataset4$spc5<- binning(dataset4$spc5, bin.size = 5) 

WSSL_hy <- new("hyperSpec",spc=as.matrix(dataset4$spc5),
               wavelength = as.numeric(colnames(dataset4$spc5)),
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
#Center and auto-scaled matric######
#########################################################################

#dataset4$spc5 <- scale(dataset4$spc5,center = T, scale = T)
#sum(dataset4$spc5^2)
# dataset4$spc5 <- scale(dataset4$spc4,center = T, scale = F)
# dataset4$spc5 <- scale(dataset4$spc5,
#                        center = apply(dataset4$spc5,2, median), 
#                        scale = apply(dataset4$spc5, 2, mad))
# dataset4$spc5<- binning(dataset4$spc4, bin.size = 5) 
# dataset4$spc5 <- scale(dataset4$spc5, center = apply(dataset4$spc5,2, median),
#                         scale = apply(dataset4$spc5, 2, mad)) 

WSSL_hy <- new("hyperSpec",spc=as.matrix(dataset4$spc5),
               wavelength = as.numeric(colnames(dataset4$spc5)),
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

#install.packages("mvoutlier")
library(mvoutlier)
x.out <- pcout(dataset4$spc5,makeplot=TRUE)
outliers <- as.character(WSSL_hy$ID[x.out$x.dist1>10])


dataset4<- subset(dataset4,!ID %in% (outliers))

