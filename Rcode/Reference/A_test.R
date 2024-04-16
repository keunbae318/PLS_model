dataset5 <- subset(dataset4, Depth=="15-20 cm")

WSSL_hy <- new("hyperSpec",spc=as.matrix(dataset5$spc),
               wavelength = as.numeric(colnames(dataset5$spc)),
               data=dataset5[,1:11],
               label=list(.wavelength="nm",spc="Absorbance"))


WSSL_hy <- new("hyperSpec",spc=as.matrix(dataset4$spc5),
               wavelength = as.numeric(colnames(dataset4$spc5)),
               data=dataset4[,1:11],
               label=list(.wavelength="nm",spc="Absorbance"))

WSSL_hy
 


###Plotting for whole data##
ggplot(as.long.df(WSSL_hy),aes(x=.wavelength,y=spc,color=Site,group=ID))+
  geom_line()+
  facet_grid(Treatment~Depth)+
  theme(legend.position = "top")+
  theme_bw()+
  xlab(bquote('Wavenumber '(cm^-1)))+
  ylab('Absorbance')+
  scale_x_reverse()+
  xlim(3800, 3650)


ggplot(as.long.df(WSSL_hy),aes(x=.wavelength,y=spc,color=Site,group=ID))+
  geom_line()+
  facet_grid(Treatment~Depth)+
  theme(legend.position = "top")+
  theme_bw()+
  xlab(bquote('Wavenumber '(cm^-1)))+
  ylab('Absorbance')+
  scale_x_reverse()+
  xlim(2000,1700)


ggplot(as.long.df(WSSL_hy),aes(x=.wavelength,y=spc,color=Site,group=ID))+
  geom_line()+
  facet_grid(Treatment~Depth)+
  theme(legend.position = "top")+
  theme_bw()+
  xlab(bquote('Wavenumber '(cm^-1)))+
  ylab('Absorbance')+
  scale_x_reverse()+
  xlim(900,1550)

