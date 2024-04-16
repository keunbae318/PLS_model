######################################################################
#Model developing to confirm if it susscessfully infers the parameters
######################################################################
library(chemometrics)
library(mvoutlier)
library(resemble)
WSSL_trans_spc <- movav(dataset4$spc,11)
WSSL_trans_spc <- savitzkyGolay(WSSL_trans_spc,m = 1,p = 5,w = 21)
WSSL_trans_spc <- scale(WSSL_trans_spc,
                        center= apply(WSSL_trans_spc,2,median),
                        scale=apply(WSSL_trans_spc,2,mad))


Xpca <- princomp(WSSL_trans_spc)
res <- pcaDiagplot(WSSL_trans_spc,Xpca,a=5,quantile = .999,plot = F)

#######################################################################
WSSL_spc<- movav(dataset4$spc,11)
scv<- cov(WSSL_spc)
scvInv=solve(scv)

differ <- as.matrix(sweep(x=WSSL_spc,MARGIN = 2,FUN = "-",
                          STATS = colMeans(WSSL_spc)))

mDist <- (rowSums(((differ) %*% scvInv )* (differ)))^0.5
mDist <- (rowSums(((differ) %*% scvInv) * (differ)))^0.5

######################################################################
PC <- ortho_projection(Xr=WSSL_spc,method="pca",center=TRUE,scales=TRUE)
dim(PC$scores)
tree <- hclust(dist(PC$scores))
plot(tree)


cl <- as.factor(cutree(tree,k=7))
g1 <- plotSpectra(spc=WSSL_spc,by=cl,by.wrap=F)


#######################################################################
dataset4$test <- movav(dataset4$spc,11)
dataset4$test <- savitzkyGolay(dataset4$test,p=3,w=21,m=0)
WSSL_spc <- msc(X=as.matrix(dataset4$test))

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
