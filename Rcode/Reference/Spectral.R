citation(package="prospectr")

#Moving average or running mean

View(dataset$spc)
#plot the first spectrum
plot(x=as.numeric(colnames(dataset$spc)),
     y=dataset$spc[5,],
     type = "l",
     lwd=1.5,
     xlab = "Wavelegnth",
     ylab ="Absorbance")

X<-movav(X=dataset$spc,w=4) #window size of 11 bands

View(X)
plot(x=as.numeric(colnames(X$spc)),
     y=X$spc[1,],
     type = "l",
     lwd=1.5,
     xlab = "Wavelegnth",
     ylab ="Absorbance")

lines(x=as.numeric(colnames(X)),y=X[1,], lwd=1.5, col="red")
grid()
legend("topleft",
       legend=c("raw","moving average"),
       lty=c(1,1),col=c("black","red"))

#p=plynomial order w=window size (must be odd) m =m-th derivative (0=smoothing)
#the function accepts vectors, data.frames or matrices. For a matric input, 
#observations should be arranged row-wise


sgvec <- savitzkyGolay(X=dataset$spc[1,],p=3,w=11,m=0)
sg <- savitzkyGolay(X=dataset$spc,p=3,w=11,m=0)

dim(dataset$spc)
dim(sg)

#X=wavelength
#Y=spectral matrix
#n=order

d1=t(diff(t(dataset$spc),differences = 1)) #first derivative
d2=t(diff(t(dataset$spc),differences = 2)) #second derivative

plot(as.numeric(colnames(d1)),
     d1[1,],
     type = "l",
     lwd=1.5,
     xlab= "Wavelength",
     ylab="")
lines(as.numeric(colnames(d2)),d2[1,],lwd=1.5,col="red")
grid()
legend("topleft",
       legend = c("1st der","2nd der"),
       lty = c(1,1),
       col = c("black","red"))

#first derivative with a gap of 10 bands 
gd1<- t(diff(t(dataset$spc),differences = 1,lag=10))

#m = order of the derivative
#w = gap size
#s = segment size
#first derivative with a gap of 5 bands
gsd1 <-gapDer(X=dataset$spc,m=1,w=11,s=5)
plot(as.numeric(colnames(d1)),
     d1[1,],
     type="l",
     lwd =1.5,
     xlab = "Wavelength",
     ylab = "")
lines(as.numeric(colnames(gsd1)),gsd1[1,], lwd = 1.5, col= "red")     
grid()
legend("topleft",
       legend = c("1st der","gap-segment 1st der"),
       lty = c(1,1),
       col = c("black","red"))

snv <- standardNormalVariate(X=dataset$spc)

# X= input spectral matrix
msc_spc <- msc(X= dataset$spc,ref_spectrum = colMeans(dataset$spc))

plot(as.numeric(colnames(dataset$spc)),
     dataset$spc[1,],
     type = "l",
     xlab = "Wavelength, nm", ylab = "Absorbance",
     lwd= 1.5)
lines(as.numeric(colnames(dataset$spc)),
      msc_spc[1,],
      lwd = 1.5, col ="red")
grid()
legend("topleft",
       legend = c("raw","MSC signal"),
       lty= c (1,1),
       col= c("black","red"))

#X=input spectral matrix
#Wav= band centers
dt<-detrend(X=dataset$spc,wav=as.numeric(colnames(dataset$spc)))
plot(as.numeric(colnames(dataset$spc)),
                dataset$spc[1,],
                type="l",
                xlab="Wavelength",
                ylab="Absorbance",
                lwd=1.5)
par(new=TRUE)
plot(dt[1,],
     xaxt="n",
     yaxt="n",
     xlab = "",
     ylab = "",
     lwd =1.5,
     col= "red",
     type = "l")
axis(4,col = "red")
grid()
legend("topleft",
       legend=c("raw", "detrend signal"),
       lty=c(1,1),
       col=c("black","red"))
par(new=FALSE)

wav <- as.numeric(colnames(dataset$spc))
#plot of the 3 first absorbance spectra

matplot(wav,
        t(dataset$spc[1:3,]),
        type ="l",
        #ylim =c(0,0.3),
        xlab ="Wavelength /nm",
        ylab = "Absorbance")
grid()
bs <- baseline(dataset$spc,wav)
matlines(wav,t(bs[1:3,]))

fitted_baselines <- attr(bs,"baselines")
matlines(wav,t(fitted_baselines[1:3,]))

#X-spectral matric type= 'soft' or 'hard' The output is a list with the 
#scaled matrix (Xscaled) and the divisor (f)
bs <- blockScale(X=dataset$spc, type = "hard")$Xscaled
sum(apply(bs,2,var))# this works!

#X=spectral matrix targetnorm=desired norm for X
bn <- blockNorm(X=dataset$spc, targetnorm =1)$Xscaled
sum(bn^2)

sum(is.na(dataset$spc[437,]))

View(dataset$spc[437,])
plot(as.numeric(colnames(dataset$spc)),dataset$spc[436,])

is.na(dataset$spc[437,])

View(dataset$spc[437,])
matplot(dataset$spc[437,],as.numeric(colnames(dataset$spc)))
#type of data: 'R' for reflectance (default), 'A' for absorbance
cr <- continuumRemoval(X = dataset$spc[1:430,], type = "A")
cr <- continuumRemoval(X = NIRsoil$spc, type = "A")
View(NIRsoil$spc)
?continuumRemoval
matplot(as.numeric(colnames(dataset$spc[1:3,])),
        t(dataset$spc[1:3,]),
        type="l",
        lty=1,
        ylim=c(0,6),
        xlab="Wavelength /nm",
        ylab="Absorbance")

matlines(as.numeric(colnames(dataset$spc)),lty = 1, t(cr[1:3,]))


kms<-naes(X=dataset$spc,k=5,pc=2,iter.max = 100)
plot(kms$pc,col=rgb(0,0,0,0.3),pch=19,main="k-mean")
grid()
points(kms$pc[kms$model,],col="red",pch=19)


#Create a dataset for illustrating how the calibration sampling
#algorithms work
X<-data.frame(x1=rnorm(1000), x2=rnorm(1000))
plot(X,col=rgb(0,0,0,0.3),pch=19,main="Kennard-Stone (synthetic)")
grid()
#KenStone produces a list with row index of the points selected for calibration 
ken <- kenStone(X,k=40)
#plot selected points 
points(X[ken$model,], col="red",pch=19,cex=1.4)

Ken_mahal <- kenStone(X=dataset$spc, k=20,metric = "mahal",pc=2)
plot(Ken_mahal$pc[,1],
     Ken_mahal$pc[,2],
     col=rgb(0,0,0,0.3),
     pch =19,
     xlab="PC1",
     ylab="PC2",
     main="Kennard-Stone")

grid()
#This is the selected points in the pc space
points(Ken_mahal$pc[Ken_mahal$model,1],
       Ken_mahal$pc[Ken_mahal$model,2],
       pch = 19, col ="red")

#Indices of the initialization samples
initialization_ind <- c(486, 702, 722, 728)
Ken_mahal_init <- kenStone(X=dataset$spc, k= 20, metric = "mahal", pc = 2, init =initialization_ind)

plot(Ken_mahal_init$pc[,1],
     Ken_mahal_init$pc[,2],
     col=rgb(0,0,0,0.3),
     pch =19,
     xlab ="PC1",
     ylab ="PC2",
     main = "Kennard-Stone with 4 initialization samples")
