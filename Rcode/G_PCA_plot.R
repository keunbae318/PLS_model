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
library(FactoMineR)
library(ggpmisc)
library(factoextra)
library(ggsci)
library(plotly)
library(dplyr)
library(tidyr)
#------------------------------------------------------------------
#Calibratio average spectra are caluacated based on the replicates
#------------------------------------------------------------------
spectraMatrix=as.matrix(dataset$spc)
columns <- colnames(spectraMatrix) %in% c("4000":"650")
subsetMatrix <- spectraMatrix[, columns]
idMatrix=as.matrix(dataset$ID)

Averages <- SpecAver(subsetMatrix,4,idMatrix, thr = .07)


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
  #xlab(bquote('Wavenumber '(cm^-1)))+
  ylab('Absorbance')+
  scale_x_reverse()


#################################################################################3
#PCA#############################
##################################################################################

#dataset5 <- subset(dataset4, Depth=="0-5 cm")
dataset5 <- dataset4
Meta.dat <- dataset4[,1:26]

res.pca <- prcomp(dataset5$spc2,center = TRUE,scale = FALSE)

# #axis#
# options(digits=3)
# summary(res.pca)$importance[2,1]*100
# paste(options(summary(res.pca)$importance[2,1], digits=2))*100,"%",sep="")
# summary(res.pca)$importance[2,2]*100


basic_plot <- as.data.frame(res.pca$x)
basic_plot$ID <-Meta.dat$ID
basic_plot$ID <-as.factor(basic_plot$ID)
basic_plot_ordispider <- left_join(as.data.frame(basic_plot), Meta.dat, by = "ID")

cent <- aggregate(cbind(PC1, PC2) ~ Depth, data = basic_plot_ordispider, FUN = mean) #Parameter (Depth vs Treatment)
segs <- merge(basic_plot_ordispider, setNames(cent, c('Depth','PC1.avg','PC2.avg')),     #Parameter (Depth vs Treatment)
              by = 'Depth', sort = FALSE)

fviz_pca_ind(res.pca,
             axes = c(1,2),
             label = "none", # hide individual labels
             habillage = Meta.dat$Treatment, # color by groups
             #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE # Concentration ellipses
)

?fviz_pca_ind
PCA_Treatment <- ggplot(basic_plot_ordispider, aes(x = PC1, y = PC2, col = Treatment, shape= Treatment)) +
  geom_point(size =2) +
  scale_color_manual(values=c("#00AFBB","#E7B800","#FC4E07"))+
  stat_ellipse(type = "t")+
  # geom_segment(data = segs,
  #              mapping = aes(xend = PC1.avg, yend = PC2.avg),
  #              col = "gray",
  #              size =0.2) + # spiders
  # #geom_point(data = cent, size = 2) +                 # centroids
  # sample scores
  #coord_fixed()+
  theme_bw()+
  #scale_color_jco()+
  # geom_text_repel(label=basic_plot_ordispider$Site,
  #                 #nudge_x = 0,
  #                 #nudge_y = -1,
  #                 size= 2.0,
  #                 box.padding = unit(0.1, "lines"), force = 2,
  #                 segment.colour = NA,
  #                 max.overlaps = 20,
  #                 show.legend = FALSE)+
  #scale_color_brewer(palette="Dark2")+
  #scale_color_manual(values = wes_palette("GrandBudapest1", n = 3))
  
  xlab(paste0("PC1 (57.6 %)"))+
  ylab(paste0("PC2 (16.8 %)"))+
  theme(legend.position = "top",
        legend.title=element_blank())#+

PCA_Treatment



PCA_Depth <- ggplot(basic_plot_ordispider, aes(x = PC1, y = PC2, col = Depth, shape= Depth)) +
  geom_point(size =2) +
  stat_ellipse(type = "t")+
  # geom_segment(data = segs,
  #              mapping = aes(xend = PC1.avg, yend = PC2.avg),
  #              col = "gray",
  #              size =0.2) + # spiders
  # # #geom_point(data = cent, size = 2) +                 # centroids
  # sample scores
  #coord_fixed()+
  theme_bw()+
  ggsci::scale_color_jco()+
  # geom_text_repel(label=basic_plot_ordispider$Site,
  #                 #nudge_x = 0,
  #                 #nudge_y = -1,
  #                 size= 2.0,
  #                 box.padding = unit(0.1, "lines"), force = 2,
  #                 segment.colour = NA,
  #                 max.overlaps = 20,
  #                 show.legend = FALSE)+
  #scale_color_brewer(palette="Dark2")+
  #scale_color_manual(values = wes_palette("GrandBudapest1", n = 3))
  
  xlab(paste0("PC1 (57.6 %)"))+
  ylab(paste0("PC2 (16.8 %)"))+
  theme(legend.position = "top",
        legend.title=element_blank())#+

PCA_Depth

ggarrange(PCA_Treatment,PCA_Depth, nrow = 1,align = "h", 
          common.legend = FALSE,labels = c("(A)","(B)"))


##loading_plot


loadings_data <- as_tibble(res.pca$rotation) %>% 
  mutate(INDEX = row.names(res.pca$rotation)) %>% 
  relocate(INDEX)
loadings_data$INDEX <- as.numeric(loadings_data$INDEX)


loading.spc <- as_tibble(loadings_data[,1:3]) %>% 
  gather(PC1:PC2,
         key = PC_axis,
         value = Loadings)

loading.spc$INDEX <- as.numeric(loading.spc$INDEX)


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
  xlab(bquote('Wavenumber '(cm^-1)))+
  ylab('Loadings')+
  scale_x_reverse(n.breaks = 10)+
  grids(linetype = "dashed")+
  scale_color_manual(values=c("#0073C2FF", "#EFC000FF"))


q <- ggplotly(Loading1)



ggarrange(ggarrange(PCA_Treatment, PCA_Depth, ncol = 2, labels = c("A", "B")),
          Loading1,nrow = 2,labels = c("","C")) 


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
  # xlab(bquote('Wavenumber '(cm^-1)))+
   ylab('Loadings')+
   scale_x_reverse(n.breaks = 10)
   grids(linetype = "dashed")+
   scale_color_manual(values=c("#0073C2FF", "#EFC000FF"))


q <- ggplotly(Loading1)
q
          