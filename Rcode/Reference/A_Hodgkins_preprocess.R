dataset4<- dataset3[,c(1,6,2,3,4,7,8,9,10,11,12,5)]

dataset4$Depth=factor(dataset4$Depth,levels = c("T","M","B"),
                      labels = c("0-5 cm","15-20 cm","45-50 cm"))

dataset4$Treatment=factor(dataset4$Treatment,levels = c("U","R","D"),
                          labels = c("Undrained","Rewetted","Drained"))

########################################################################
#Rawdata#
#########################################################################

WSSL_hy <- new("hyperSpec",spc=as.matrix(dataset4$spc),
               wavelength = as.numeric(colnames(dataset4$spc)),
               data=dataset4[,1:11],
               label=list(.wavelength="nm",spc="Absorbance"))


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
#Noise removal n baseline correction#
#########################################################################
#moving average
dataset4$spc2 <- movav(dataset4$spc,11)
#baseline correction
wav <- as.numeric(colnames(dataset4$spc2))
dataset4$spc2 <- baseline(dataset4$spc2, wav)

#etc#
#dataset4$spc2 <- blockScale(X = dataset4$spc2, type = "soft")$Xscaled
#dataset4$spc2  <- scale(dataset4$spc2)
#dataset4$spc2 <-normalize01(dataset4$spc2)
#dataset4$spc2<- binning(dataset4$spc2, bin.size = 10) 
#dataset4$spc2 <- blockNorm(X = dataset4$spc2, targetnorm = 1)$Xscaled


WSSL_hy <- new("hyperSpec",spc=as.matrix(dataset4$spc2),
               wavelength = as.numeric(colnames(dataset4$spc2)),
               data=dataset4[,1:11],
               label=list(.wavelength="nm",spc="Absorbance"))


###Plotting for whole data##
ggplot(as.long.df(WSSL_hy),aes(x=.wavelength,y=spc,color=Site,group=ID))+
  geom_line()+
  facet_grid(Depth~Treatment)+
  theme(legend.position = "top")+
  theme_bw()+
  xlab(bquote('Wavenumber '(cm^-1)))+
  ylab('Absorbance')+
  scale_x_reverse()



#############################################################
#data export for Hodgkins#
############################################################
New.dat<-dataset4$spc2
rownames(New.dat)<-dataset4$ID
New.dat<-t(New.dat)


write.csv(New.dat,"C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/MIR/Data/Prep_Hodgkins.csv")
############################################################



