#install.packages("FactoMineR")
library(FactoMineR)
library(ggpmisc)
library(factoextra)
#install.packages("QuantPsyc")
#install.packages("readxl")
library("readxl")
library(QuantPsyc)
library(ggcorrplot)
library(ggrepel)
library("xlsx")
library(MVN)
library(tidyr)
library(ggpubr)
library(ggplot2)
library (vegan)
library(FactoMineR)
library(factoextra)
library(paran)
library(DESeq2)
library(S4Vectors)
library(stats4)
library(BiocGenerics)
library(sva)
library("ggsci")
library(wesanderson)
library(tidyr)
library(mixOmics)
library(sva) 
library(tidyr)

#install.packages("spectratrait")

###Call_pyrolysis##


setwd("C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/Pyrolysis_GCMS/Data/")
MZmine.dat<-read.csv("Merged_V2.csv",check.names = FALSE,row.names = 1)

#Meta_data##
z=c("Arlon","Zwarte_Beek","Binnenveld","Drentse_Aa","Gutzkow","Peene_mouth","Kiel","Recknitz","Biebrza",
    "Rospuda","Suwalszczyzna","Mazury","Anglesey")
Meta.dat <- data.frame(Sample_ID=paste0(rep(1:39,each=3),c("A","B","C")),
                       Treatment= rep(c("Undrained","Drained","Rewetted"),each=3),
                       Depth = rep (c("0-5 cm","15-20 cm","45-50 cm")),
                       Site = as.factor(rep (z, each=9)))
Meta.dat$Treatment <- factor(Meta.dat$Treatment, levels=c("Undrained","Rewetted","Drained"))
Meta.dat$Batch <-ifelse(Meta.dat$Sample_ID %in% c("2C","4C","5B","5C","6B","7A",
                                                  "7C","10A","12B","13A","14C","15A","15C","16C","17B","17C","19A","20A","23C",
                                                  "25C","26B","27C","32A","32C","34A","34B","36B","38C","39A"),"Batch2","Batch1") 
row.names(Meta.dat) <- Meta.dat$Sample_ID
####################################################################
###Data_transforms####
Compound.dat<-as.data.frame(t(MZmine.dat[,c(6:122)]))#whole data
##1## replacezero
Compound.dat[Compound.dat < 1000] <-0
replacezero <- function (x) {
  x[!x|is.na(x)] <- min(x[x>0],na.rm=TRUE) / 2
  return(x)
}
imputed.dat <-t(apply(Compound.dat,1,replacezero))
##2## Proportional data
Proportion.dat <- prop.table(as.matrix(imputed.dat), margin = 1) * 100
Proportion.dat <- Proportion.dat[order(match(rownames(Proportion.dat),Meta.dat$Sample_ID)),] #reorder
### Combat batch effect correction##
Combat.mod <- model.matrix(~as.factor(Depth)*as.factor(Treatment),data=as.data.frame(Meta.dat)) #model
Proportion.combat <- t(ComBat(t(Proportion.dat), batch = Meta.dat$Batch, mod = Combat.mod, par.prior = F, ref.batch="Batch1"))
row_sums <- rowSums(Proportion.combat)
rows_to_scale <- row_sums != 100
min_value <- abs(min(Proportion.combat, na.rm = TRUE))+0.001
Proportion.combat[rows_to_scale, ] <- Proportion.combat[rows_to_scale, ] + min_value
Proportion.combat[rows_to_scale, ] <- t(t(Proportion.combat[rows_to_scale, ]) / row_sums[rows_to_scale] * 100)
############################
Group.dat <- merge(Meta.dat,Proportion.combat, by ="row.names", sort=FALSE)
Group1.dat <- Group.dat %>% 
  gather(7:81, key = Group, value=Percentage)
for (i in 1:nrow(Group1.dat)){
  compound_value <- Group1.dat[i,"Group"]
  if (compound_value %in% rownames(MZmine.dat)){
    updated_value <-MZmine.dat[compound_value, "Group"]
    Group1.dat[i,"Group"] <-updated_value
  }
}
Group2.dat <- aggregate(Percentage ~ Group*Sample_ID*Treatment*Depth*Site, data= Group1.dat, FUN = sum) 
Group2.dat <- Group2.dat[,-c(3,4,5)]
Group3.dat <- Group2.dat %>% 
  tidyr::spread(key = Group, value =Percentage)
Group3.dat<-Group3.dat[match(Group.dat$Sample_ID,Group3.dat$Sample_ID),]
rownames(Group3.dat)<-1:nrow(Group3.dat)
Group3.dat <- Group3.dat[,-1]
Group3.dat

export.dat <- data.frame(Meta.dat,Group3.dat,Group.dat[,-(1:6)],check.names = FALSE)
########################################################################################

export.dat$Guaiacols <- rowSums(export.dat[, c("160", "200", "224", "235", "252","267","271")], na.rm = TRUE)
export.dat$Syringols <- rowSums(export.dat[, c("242", "288", "305", "309", "315")], na.rm = TRUE)

export.dat$Intact_lignin <- export.dat$'267'/export.dat$Guaiacols #G3/G
export.dat$Sedges <- export.dat$'235'/export.dat$Guaiacols #4-vG/G
export.dat$Cellulose <- export.dat$'274'/export.dat$Carbohydrates #levoglucosan/carbohydrate
export.dat <-export.dat[,c(1:15, 91:95)]
#write.csv(export.dat,"C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/Mid_Infrared/Data/Proportion_pyrolysis_V2.csv")

str(export.dat)

########################################################################################
#devtools::install_github(repo = "TESTgroup-BNL/spectratrait", dependencies=TRUE)
### Setup other functions and options
# not in
`%notin%` <- Negate(`%in%`)
opar <- par(no.readonly = T)
colnames(Pyrolysis)



################################################
Start.wave <- 3997
End.wave <- 403
wv <- seq(Start.wave,End.wave,-1)
###################################################
MIR <-as.matrix(dataset4$spc2)
colnames(MIR) <-c(paste0("Wave_",3997:403))
rownames(MIR)<-rownames(dataset4)
Pyrolysis<-as.matrix(Group3.dat)
######################################################


Spectra.dat<-data.frame(Meta.dat,Pyrolysis,MIR,check.names = FALSE)

str(Spectra.dat)

Sample_info <- Spectra.dat %>% 
  dplyr::select(Sample_ID,Treatment,Depth,Carbohydrates)  ####choose compounds

inVar <- "Carbohydrates"
plsr_data <-data.frame(Sample_info,MIR)

#####



str(plsr_data)
# store matrices into training and test set:
# rownames(Meta.dat) <- 1:nrow((Meta.dat))
# twospectra.dat <-data.frame(Pyrolysis = I(Pyrolysis), MIR = I(MIR), Meta = I(Meta.dat))



split_data <- spectratrait::create_data_split(dataset=plsr_data, approach="dplyr", 
                                              split_seed=7529075, prop=0.8, 
                                              group_variables="Treatment")


cal.plsr.data <- split_data$cal_data
val.plsr.data <- split_data$val_data
str(cal.plsr.data)

print(paste("Cal observations: ",dim(cal.plsr.data)[1],sep=""))
print(paste("Val observations: ",dim(val.plsr.data)[1],sep=""))
#inVar <- "Carbohydrates"

text_loc <- c(max(hist(cal.plsr.data[,paste0(inVar)], plot=FALSE)$counts),
              max(hist(cal.plsr.data[,paste0(inVar)], plot=FALSE)$mids))
cal_hist_plot <- qplot(cal.plsr.data[,paste0(inVar)],geom="histogram",
                       main = paste0("Calibration Histogram for ",inVar),
                       xlab = paste0(inVar),ylab = "Count",fill=I("grey50"),col=I("black"),
                       alpha=I(.7))# + 
  annotate("text", x=text_loc[2], y=text_loc[1], label= "1.",size=10)
val_hist_plot <- qplot(val.plsr.data[,paste0(inVar)],geom="histogram",
                       main = paste0("Validation Histogram for ",inVar),
                       xlab = paste0(inVar),ylab = "Count",fill=I("grey50"),col=I("black"),
                       alpha=I(.7))
histograms <- grid.arrange(cal_hist_plot, val_hist_plot, ncol=2)

dev.off()

################################################################################

cal_spec <- as.matrix(cal.plsr.data[, which(names(cal.plsr.data) %in% 
                                              paste0("Wave_",wv))])
cal.plsr.data <- data.frame(cal.plsr.data[, which(names(cal.plsr.data) %notin% 
                                                    paste0("Wave_",wv))], 
                            Spectra=I(cal_spec))
val_spec <- as.matrix(val.plsr.data[, which(names(val.plsr.data) %in% 
                                              paste0("Wave_",wv))])
val.plsr.data <- data.frame(val.plsr.data[, which(names(val.plsr.data) %notin% 
                                                    paste0("Wave_",wv))],
                            Spectra=I(val_spec))


f.plot.spec(Z=cal.plsr.data$Spectra,wv=wv,
                          plot_label="Calibration",
                          type = "Absorbance")
######################################################################################
library(pls)
library(spectratrait)

if(grepl("Windows", sessionInfo()$running)){
  pls.options(parallel = NULL)
} else {
  pls.options(parallel = parallel::detectCores()-1)
}

method <- "pls" #pls, firstPlateau, firstMin
random_seed <- 1245565
seg <- 50
maxComps <- 16
iterations <- 80
prop <- 0.70


if (method=="pls") {
  # pls package approach - faster but estimates more components....
  nComps <- spectratrait::find_optimal_components(dataset=cal.plsr.data, targetVariable=inVar,
                                                  method=method, 
                                                  maxComps=maxComps, seg=seg, 
                                                  random_seed=random_seed)
  print(paste0("*** Optimal number of components: ", nComps))
} else {
  nComps <- spectratrait::find_optimal_components(dataset=cal.plsr.data, targetVariable=inVar,
                                                  method=method, 
                                                  maxComps=maxComps, iterations=iterations, 
                                                  seg=seg, prop=prop, 
                                                  random_seed=random_seed)
}


plsr.out <- plsr(as.formula(paste(inVar,"~","Spectra")),scale=FALSE,ncomp=nComps,validation="LOO",
                 trace=FALSE,data=cal.plsr.data)
fit <- plsr.out$fitted.values[,1,nComps]
pls.options(parallel = NULL)

# External validation fit stats
par(mfrow=c(1,2)) # B, L, T, R
pls::RMSEP(plsr.out, newdata = val.plsr.data)

plot(pls::RMSEP(plsr.out,estimate=c("test"),newdata = val.plsr.data), main="MODEL RMSEP",
     xlab="Number of Components",ylab="Model Validation RMSEP",lty=1,col="black",cex=1.5,lwd=2)
box(lwd=2.2)

R2(plsr.out, newdata = val.plsr.data)



plot(pls::R2(plsr.out,estimate=c("test"),newdata = val.plsr.data), main="MODEL R2",
     xlab="Number of Components",ylab="Model Validation R2",lty=1,col="black",cex=1.5,lwd=2)
box(lwd=2.2)



cal.plsr.output <- data.frame(cal.plsr.data[, which(names(cal.plsr.data) %notin% "Spectra")], 
                             PLSR_Predicted=fit,
                             PLSR_CV_Predicted=as.vector(plsr.out$validation$pred[,,nComps]))
cal.plsr.output <- cal.plsr.output %>%
  mutate(PLSR_CV_Residuals = PLSR_CV_Predicted-get(inVar))
cal.plsr.output

cal.R2 <- round(pls::R2(plsr.out,intercept=F)[[1]][nComps],2)
cal.RMSEP <- round(sqrt(mean(cal.plsr.output$PLSR_CV_Residuals^2)),2)

val.plsr.output <- data.frame(val.plsr.data[, which(names(val.plsr.data) %notin% "Spectra")],
                              PLSR_Predicted=as.vector(predict(plsr.out, 
                                                               newdata = val.plsr.data, 
                                                               ncomp=nComps, type="response")[,,1]))
val.plsr.output <- val.plsr.output %>%
  mutate(PLSR_Residuals = PLSR_Predicted-get(inVar))
head(val.plsr.output)



val.R2 <- round(pls::R2(plsr.out,newdata=val.plsr.data,intercept=F)[[1]][nComps],2)
val.RMSEP <- round(sqrt(mean(val.plsr.output$PLSR_Residuals^2)),2)

rng_quant <- quantile(cal.plsr.output[,inVar], probs = c(0.001, 0.999))
cal_scatter_plot <- ggplot(cal.plsr.output, aes(x=PLSR_CV_Predicted, y=get(inVar))) + 
  theme_bw() + geom_point() + geom_abline(intercept = 0, slope = 1, color="dark grey", 
                                          linetype="dashed", size=1.5) + xlim(rng_quant[1], 
                                                                              rng_quant[2]) + 
  ylim(rng_quant[1], rng_quant[2]) +
  labs(x=paste0("Predicted ", paste(inVar), " (units)"),
       y=paste0("Observed ", paste(inVar), " (units)"),
       title=paste0("Calibration: ", paste0("Rsq = ", cal.R2), "; ", paste0("RMSEP = ", 
                                                                            cal.RMSEP))) +
  theme(axis.text=element_text(size=18), legend.position="none",
        axis.title=element_text(size=20, face="bold"), 
        axis.text.x = element_text(angle = 0,vjust = 0.5),
        panel.border = element_rect(linetype = "solid", fill = NA, size=1.5))

cal_resid_histogram <- ggplot(cal.plsr.output, aes(x=PLSR_CV_Residuals)) +
  geom_histogram(alpha=.5, position="identity") + 
  geom_vline(xintercept = 0, color="black", 
             linetype="dashed", size=1) + theme_bw() + 
  theme(axis.text=element_text(size=18), legend.position="none",
        axis.title=element_text(size=20, face="bold"), 
        axis.text.x = element_text(angle = 0,vjust = 0.5),
        panel.border = element_rect(linetype = "solid", fill = NA, size=1.5))

rng_quant <- quantile(val.plsr.output[,inVar], probs = c(0.001, 0.999))
val_scatter_plot <- ggplot(val.plsr.output, aes(x=PLSR_Predicted, y=get(inVar))) + 
  theme_bw() + geom_point() + geom_abline(intercept = 0, slope = 1, color="dark grey", 
                                          linetype="dashed", size=1.5) + xlim(rng_quant[1], 
                                                                              rng_quant[2]) + 
  ylim(rng_quant[1], rng_quant[2]) +
  labs(x=paste0("Predicted ", paste(inVar), " (units)"),
       y=paste0("Observed ", paste(inVar), " (units)"),
       title=paste0("Validation: ", paste0("Rsq = ", val.R2), "; ", paste0("RMSEP = ", 
                                                                           val.RMSEP))) +
  theme(axis.text=element_text(size=18), legend.position="none",
        axis.title=element_text(size=20, face="bold"), 
        axis.text.x = element_text(angle = 0,vjust = 0.5),
        panel.border = element_rect(linetype = "solid", fill = NA, size=1.5))

val_resid_histogram <- ggplot(val.plsr.output, aes(x=PLSR_Residuals)) +
  geom_histogram(alpha=.5, position="identity") + 
  geom_vline(xintercept = 0, color="black", 
             linetype="dashed", size=1) + theme_bw() + 
  theme(axis.text=element_text(size=18), legend.position="none",
        axis.title=element_text(size=20, face="bold"), 
        axis.text.x = element_text(angle = 0,vjust = 0.5),
        panel.border = element_rect(linetype = "solid", fill = NA, size=1.5))

# plot cal/val side-by-side
scatterplots <- grid.arrange(cal_scatter_plot, val_scatter_plot, cal_resid_histogram, 
                             val_resid_histogram, nrow=2,ncol=2)

vips <- spectratrait::VIP(plsr.out)[nComps,]

par(mfrow=c(2,1))
plot(plsr.out$coefficients[,,nComps], x=wv,xlab="Wavelength (nm)",
     ylab="Regression coefficients",lwd=2,type='l')
box(lwd=2.2)
plot(wv, vips, xlab="Wavelength (nm)",ylab="VIP",cex=0.01)
lines(wv, vips, lwd=3)
abline(h=0.8, lty=2, col="dark grey")
box(lwd=2.2)
par(opar)






rmsep_percrmsep <- spectratrait::percent_rmse(plsr_dataset = val.plsr.output, 
                                              inVar = inVar, 
                                              residuals = val.plsr.output$PLSR_Residuals, 
                                              range="full")
RMSEP <- rmsep_percrmsep$rmse
perc_RMSEP <- rmsep_percrmsep$perc_rmse
r2 <- round(pls::R2(plsr.out, newdata = val.plsr.data, intercept=F)$val[nComps],2)
expr <- vector("expression", 3)
expr[[1]] <- bquote(R^2==.(r2))
expr[[2]] <- bquote(RMSEP==.(round(RMSEP,2)))
expr[[3]] <- bquote("%RMSEP"==.(round(perc_RMSEP,2)))
rng_vals <- c(min(val.plsr.output$LPI), max(val.plsr.output$UPI))
par(mfrow=c(1,1), mar=c(4.2,5.3,1,0.4), oma=c(0, 0.1, 0, 0.2))
plotrix::plotCI(val.plsr.output$PLSR_Predicted,val.plsr.output[,inVar], 
                li=val.plsr.output$LPI, ui=val.plsr.output$UPI, gap=0.009,sfrac=0.000, 
                lwd=1.6, xlim=c(rng_vals[1], rng_vals[2]), ylim=c(rng_vals[1], rng_vals[2]), 
                err="x", pch=21, col="black", pt.bg=scales::alpha("grey70",0.7), scol="grey80",
                cex=2, xlab=paste0("Predicted ", paste(inVar), " (units)"),
                ylab=paste0("Observed ", paste(inVar), " (units)"),
                cex.axis=1.5,cex.lab=1.8)

abline(0,1,lty=2,lw=2)
plotrix::plotCI(val.plsr.output$PLSR_Predicted,val.plsr.output[,inVar], 
                li=val.plsr.output$LCI, ui=val.plsr.output$UCI, gap=0.009,sfrac=0.004, 
                lwd=1.6, xlim=c(rng_vals[1], rng_vals[2]), ylim=c(rng_vals[1], rng_vals[2]), 
                err="x", pch=21, col="black", pt.bg=scales::alpha("grey70",0.7), scol="black",
                cex=2, xlab=paste0("Predicted ", paste(inVar), " (units)"),
                ylab=paste0("Observed ", paste(inVar), " (units)"),
                cex.axis=1.5,cex.lab=1.8, add=T)
legend("topleft", legend=expr, bty="n", cex=1.5)
legend("bottomright", legend=c("Prediction Interval","Confidence Interval"), 
       lty=c(1,1), col = c("grey80","black"), lwd=3, bty="n", cex=1.5)
box(lwd=2.2)


f.coef.valid(val.plsr.output, val.plsr.data, ncomp, inVar)




spectratrait::f.plot.coef(Z = t(bootstrap_coef), wv = wv, 
                          plot_label="Bootstrap regression coefficients", 
                          position = 'bottomleft')
abline(h=0,lty=2,col="grey50")
box(lwd=2.2)
method <- "firstMin" #pls, firstPlateau, firstMin
random_seed <- 7529075
seg <- 80
maxComps <- 16
iterations <- 50
prop <- 0.70
paste(inVar,"~","Spectra")
as.formula(paste(inVar,"~","Spectra"))
pls::plsr(as.formula("Aliphatics~Spectra"), scale=FALSE, center=TRUE, ncomp=maxComps, 
          validation="CV", segments = seg, segment.type="interleaved", trace=FALSE, 
          jackknife=TRUE, data=cal.plsr.data)
str(cal.plsr.data)

str(cal.plsr.data$Spectra)[1:5]

plsr(as.formula("Aliphatics~Spectra"))

if (method=="pls") {
  nComps <- spectratrait::find_optimal_components(dataset=cal.plsr.data, method=method, 
                                                  maxComps=maxComps, seg=seg, 
                                                  random_seed=random_seed)
  print(paste0("*** Optimal number of components: ", nComps))
} else {
  nComps <- spectratrait::find_optimal_components(dataset=cal.plsr.data, method=method, 
                                                  maxComps=maxComps, iterations=iterations, 
                                                  seg=seg, prop=prop, 
                                                  random_seed=random_seed)
}

spectratrait::find_optimal_components(dataset=cal.plsr.data, method=method, 
                        maxComps=maxComps, iterations=iterations, 
                        seg=seg, prop=prop, 
                        random_seed=random_seed)

?find_optimal_components

spectratrait::find_optimal_components(dataset=cal.plsr.data, method="firstMin", 
                                      maxComps=16, seg=80, 
                                      random_seed=7529075)

str(cal.plsr.data)
cal.plsr.data[lapply(cal.plsr.data,length)>0]

train <- sample(1:nrow(twospectra.dat), 80) # randomly select 50 samples in training
test <- setdiff(1:nrow(twospectra.dat), train) # rest is part of the test set

train.dat <- twospectra.dat[train, ]
test.dat <- twospectra.dat[test,]


MIR_PLSR <-plsr(Pyrolysis~MIR, ncomp = 10, 
                data = train.dat,validation ="LOO")

summary(MIR_PLSR)
validationplot(MIR_PLSR)

validationplot(MIR_PLSR, val.type="MSEP")
validationplot(MIR_PLSR, val.type="R2")
nmr.pred.test <- as.data.frame(predict(MIR_PLSR, ncomp = 6,newdata = test.dat))
nmr.pred.test

validationplot(MIR_PLSR, val.type="RMSEP")
plot(MIR_PLSR, ncomp=11, asp=1, line=TRUE)


sqrt(mean((nmr.pred.test - test.dat)^2))

library(mdatools)
plot(MIR_PLSR)

splotRegcoeffs(MIR_PLSR)
