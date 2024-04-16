Dataset.dat <- dataset4

###Reanalyiss#######
# rm(list = ls())
# getwd()
# setwd("R:/Guillermo's Lab/POM/Research/POM/POM Data/Data statical analysis/Sequence_POM")
# Dataset.dat=read.csv("Dataset.csv")

##Explore data#######################
str(Dataset.dat)


###Orgonize dataframe###################
Dataset.dat$Site=as.factor(Dataset.dat$Site)

levels(Dataset.dat$Treatment)
table(Dataset.dat$Treatment)
table(Dataset.dat$Treatment,Dataset.dat$Depth)
table(droplevels(Dataset.dat)$Treatment)
##library file#######################################
#install.packages("corrplot")
#install.packages("PerformanceAnalytics")
#install.packages("ggcorrplot")

library(corrplot)
library(PerformanceAnalytics)
library(ggcorrplot)
str(Dataset.dat)
###############################################################
names(Dataset.dat)[names(Dataset.dat) == "CN"] <- "C:N"
################################################################
# mat : is a matrix of data
# ... : further arguments to pass to the native R cor.test function
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}
# matrix of the p-value of the correlation
str(Sub.dat)

########################################################
Sub.dat=subset(Dataset.dat,Depth=="0-5 cm")
Sub.dat=subset(Sub.dat,select = -c(ID,Country,Site,Treatment,Depth,spc,spc2,spc3,spc4,spc5))
p.mat <- cor.mtest(Sub.dat)
M=cor(na.omit(Sub.dat),method="spearman")

a=corrplot(M,type = "upper",
           p.mat = p.mat, sig.level = 0.05, insig = "blank",diag = FALSE)


#######################################################
Sub.dat=subset(Dataset.dat,Depth=="15-20 cm")
Sub.dat=subset(Sub.dat,select = -c(ID,Country,Site,Treatment,Depth,spc,spc2,spc3,spc4,spc5))
p.mat <- cor.mtest(Sub.dat)
M=cor(na.omit(Sub.dat),method="spearman")

b=corrplot(M,type = "upper",
           p.mat = p.mat, sig.level = 0.05, insig = "blank",diag = FALSE)

#######################################################
Sub.dat=subset(Dataset.dat,Depth=="45-50 cm")
Sub.dat=subset(Sub.dat,select = -c(ID,Country,Site,Treatment,Depth,spc,spc2,spc3,spc4,spc5))
p.mat <- cor.mtest(Sub.dat)
M=cor(na.omit(Sub.dat),method="spearman")

c=corrplot(M,type = "upper",
           p.mat = p.mat, sig.level = 0.05, insig = "blank",diag = FALSE)


#######################################################
Sub.dat=subset(Dataset.dat,select = -c(ID,Country,Site,Treatment,Depth,spc,spc2,spc3,spc4,spc5))
p.mat <- cor.mtest(Sub.dat)
M=cor(na.omit(Sub.dat),method="spearman")

d=corrplot(M,type = "upper",
           p.mat = p.mat, sig.level = 0.05, insig = "blank",diag = FALSE)


a
b
c
d
