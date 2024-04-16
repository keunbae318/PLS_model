library(mixOmics)
library(pls)
library(prospectr)
set.seed(5249)





spectra.data<-dataset4
spectra.data$spc <- standardNormalVariate(dataset4$spc)

spectra.data$spc2 <- movav(spectra.data$spc,11)
spectra.data$spc2<- binning(spectra.data$spc2, bin.size = 10)
spectra.data$spc2 <- standardNormalVariate(spectra.data$spc2)


# spectra.data_trained <- subset(spectra.data,spectra.data$Treatment == "Drained"|
#                                  spectra.data$Treatment == "Undrained")
# 
# spectra.data_trained$Treatment <- factor(spectra.data_trained$Treatment, levels = c("Undrained","Drained"))

#rownames(spectra.data_trained) <-1:nrow(spectra.data_trained)

X <- spectra.data$spc2 
Y2 <- spectra.data$Treatment

design <- data.frame(sample = spectra.data$Site)


###############################
pca.spectra <- pca(X, ncomp = 10, center = TRUE, scale = TRUE)
pca.multilevel.spectra.depth <- pca(X, ncomp =10, center = TRUE, scale = TRUE,
                              multilevel = design)
plot(pca.spectra)

#PCA plot##
plotA_graph <-
plotIndiv(pca.spectra, group = spectra.data$Treatment, ind.names = FALSE,
          ellipse = TRUE,
          legend = TRUE, title = "PCA on MIR comp 1-2")



#PCA plot with multilevel##
plotIndiv(pca.multilevel.spectra.depth, group = spectra.data$Treatment,
          ind.names = FALSE,
          ellipse = TRUE,
          legend = TRUE, legend.title = 'Stimulation',
          title = '(b) multilevel PCA on MIR comp 1-2')

##############################
spectra.splsda.depth <- splsda(X,Y2, ncomp = 15)
spectra.multilevel.splsda.depth <- splsda(X,Y2, ncomp = 15, multilevel = design)

#PLSDA plot##
plotIndiv(spectra.splsda.depth, comp = 1:2,
          group = spectra.data$Treatment, ind.names = FALSE,
          ellipse = TRUE,
          legend = TRUE, title = "(c) PLSDA with confidence ellipses")

background.depth <- background.predict(spectra.splsda.depth, comp.predicted = 2, dist = "max.dist")

plotIndiv(spectra.splsda.depth, comp= 1:2,
          group = spectra.data$Treatment, ind.names = FALSE,
          background = background.depth,
          legend = TRUE, title = "(d) PLSDA with prediction background")

#PLSDA plot with multilevel##
plotIndiv(spectra.multilevel.splsda.depth, comp =1:2,
          group = spectra.data$Treatment, ind.names = FALSE,
          ellipse =TRUE,
          legend = TRUE, title = "(e) PLSDA multilevl with confidence ellipses")

background.multilevel.depth <- background.predict(spectra.multilevel.splsda.depth, comp.predicted = 2, dist = "max.dist")

plotIndiv(spectra.multilevel.splsda.depth, comp = 1:2, 
          group = spectra.data$Treatment, ind.names = FALSE,
          background = background.multilevel.depth,
          #ellipse = TRUE,
          legend = TRUE, title = "(f) PLSDA multilevel with prediction background")
#################################
#Tuning sPLS-DA
#################################
if(grepl("Windows", sessionInfo()$running)){
  pls.options(parallel = NULL)
} else {
  pls.options(parallel = parallel::detectCores()-1)
}

#Option1
perf.splsda.spectra.depth2 <- perf(spectra.splsda.depth, validation = "Mfold",
                            folds = 5, nrepeat = 50,
                            progressBar = TRUE, auc = TRUE,
                            cpus = 10)

perf.splsda.multilevel.spectra.depth2 <- perf(spectra.multilevel.splsda.depth, validation = "Mfold",
                                       folds = 5, nrepeat =50,
                                       progressBar = TRUE, auc =TRUE,
                                       cpus = 10)

#Option2
# perf.splsda.spectra <- perf(spectra.splsda, validation = "loo",
#                             progressBar = FALSE, auc = TRUE,
#                             cpus = 15)

perf.splsda.spectra.depth


plot(perf.splsda.spectra.depth, sd = TRUE, 
     col = color.mixo(5:7),
     legend.position ='horizontal')

plot(perf.splsda.multilevel.spectra.depth2, sd =TRUE,
     col = color.mixo(5:7),
     legend.position = 'horizontal')

perf.splsda.spectra.depth$choice.ncomp
perf.splsda.multilevel.spectra.depth2$choice.ncomp

perf.splsda.spectra.depth$error.rate$overall[,'max.dist']
perf.splsda.multilevel.spectra.depth2$error.rate$overall[,'max.dist']


#Final PLS-DA model
final.plsda.spectra.depth2 <-plsda(X,Y2, ncomp = 3)

plotIndiv(final.plsda.spectra.depth, ind.names = FALSE, legend = TRUE,
          comp = c(1,2), ellipse = TRUE,
          title = 'PLS-DA on Spectra comp 1-2',
          X.label = 'PLS-DA comp 1', Y.label = 'PLS-DA comp 2')

plotIndiv(final.plsda.spectra.depth2, ind.names = FALSE, legend = TRUE,
          comp = c(2,3), ellipse = TRUE,
          title = 'PLS-DA on Spectra comp 2-3',
          X.label = 'PLS-DA comp 2', Y.label = 'PLS-DA comp 3')

final.plsda.multilevel.spectra.depth2 <-plsda(X,Y2, ncomp = 3,
                                       multilevel = design)


plotIndiv(final.plsda.multilevel.spectra.depth2, ind.names = FALSE, legend = TRUE,
          comp = c(1,2), 
          ellipse = TRUE,
          title = '(c) multilevel PLS-DA on Spectra comp 1-2',
          #X.label = 'PLS-DA comp 1', Y.label = 'PLS-DA comp 2'
)

plotIndiv(final.plsda.multilevel.spectra.depth2, ind.names = FALSE, legend = TRUE,
          comp = c(2,3), ellipse = TRUE,
          title = 'PLS-DA on Spectra comp 2-3',
          X.label = 'PLS-DA comp 2', Y.label = 'PLS-DA comp 3')
#Selecting the number of variables

list.keepX <- c(1:10, seq(20,300,10))
tune.splsda.spectra <-tune.splsda( X, Y2, ncomp = 10,
                                   validation = 'Mfold',
                                   folds =5, nrepeat = 50,
                                   dist = 'max.dist',
                                   measure = "overall",
                                   test.keepX = list.keepX,
                                   cpus = 10)

tune.splsda.multilevel.spectra2 <-tune.splsda( X, Y2, ncomp = 3,
                                              multilevel = design,
                                              validation = 'Mfold',
                                              folds =5, nrepeat = 50,
                                              dist = 'max.dist',
                                              measure = "overall",
                                              test.keepX = list.keepX,
                                              cpus = 10,
                                              progressBar = TRUE)


# tune.splsda.spectra <-tune.splsda( X, Y, ncomp = 10,
#                                    validation = 'Mfold',
#                                    folds =10, nrepeat = 100,
#                                    dist = 'mahalanobis.dist',
#                                    measure = "overall",
#                                    test.keepX = list.keepX,
#                                    cpus = 15)

plot(tune.splsda.multilevel.spectra,col=color.jet(3))

tune.splsda.multilevel.spectra$choice.ncomp$ncomp
tune.splsda.multilevel.spectra$choice.keepX

optimal.ncomp <- tune.splsda.multilevel.spectra$choice.ncomp$ncomp
optimal.keepX <- tune.splsda.multilevel.spectra$choice.keepX[1:optimal.ncomp]


final.splsda <- splsda (X,Y2,
                        ncomp = optimal.ncomp,
                        keepX = optimal.keepX)



final.splsda.multilevel2 <- splsda(X,Y2,
                                  ncomp = optimal.ncomp,
                                  keepX = optimal.keepX,
                                  multilevel = design
                                  )

plotIndiv(final.splsda, comp = c(1,2),group = spectra.data$Treatment, ind.names = FALSE,
          ellipse = TRUE, legend = TRUE,
          title = '(a) sPLS-DA on MIR, comp 1 & 2')


plotb_graph <-
plotIndiv(final.splsda.multilevel2, comp = c(1,2), group = spectra.data$Treatment,
          ind.names = FALSE,
          ellipse = TRUE, 
          legend = TRUE, 
          legend.title = 'Treatment',
          title = '(D) Multilevel sPLS-DA on MIR, comp 1 & 2')

ggarrange(plotA_graph,plotb_graph, ncol = 2)

dev.off()

legend = list(legend = levels(Y2),
              col=unique(color.mixo(Y2)),
              title = "Degradation",
              cex = 0.7)

cim <- cim(final.splsda.multilevel, row.sideColors = color.mixo(Y2),
           legend = legend)

plotVar(final.splsda, comp = c(1,2))
plotVar(final.splsda.multilevel, comp =c(1,2))


plotLoadings(final.splsda.multilevel, comp = 1) 
plotLoadings(final.splsda.multilevel2, comp = 1 ,contrib = 'max', method = 'median',study="global")
?plotLoadings


