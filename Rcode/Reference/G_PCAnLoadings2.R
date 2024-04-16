######################
####Data subset
#######################
dataset5 <- dataset4
dataset5 <- subset(dataset4, Depth=="0-5 cm")
dataset5 <- subset(dataset4, Depth=="15-20 cm")
dataset5 <- subset(dataset4, Depth=="45-50 cm")

dataset5 <- subset(dataset4, Treatment=="Undrained")
dataset5 <- subset(dataset4, Treatment=="Rewetted")
dataset5 <- subset(dataset4, Treatment=="Drained")

#############################################
###################################################################################

dataset5 <- subset(dataset4, Depth=="0-5 cm")
env.dat <- dataset5[,7:41] 
#####################################################################################
colnames(env.dat)[which(names(env.dat) == "carb")] <- "Carbohyrates_rel_abund"
colnames(env.dat)[which(names(env.dat) == "arom15")] <- "arom15_rel_abund"
colnames(env.dat)[which(names(env.dat) == "arom16")] <- "arom16_or_COO_rel_abund"
colnames(env.dat)[which(names(env.dat) == "acids")] <- "Acids_rel_abund"
colnames(env.dat)[which(names(env.dat) == "aliph28")] <- "Aliph28_PI"
colnames(env.dat)[which(names(env.dat) == "aliph29")] <- "Aliph29_PI"
colnames(env.dat)[which(names(env.dat) == "sum_arom")] <- "Sum_Aromatics_rel_abund"
colnames(env.dat)[which(names(env.dat) == "sum_aliph")] <- "Aliphatic_rel_abund"
#######################################################################################

res.pca <- prcomp(dataset5$spc2,center = TRUE,scale = FALSE)
# quanti.sup <- subset(env.dat, select =c(Carbohyrates_rel_abund,
#                                         Aromatics_rel_abund,
#                                         Aliphatic_rel_abund,
#                                         C,CN,d15N,BD,Mineral))
str(env.dat)
quanti.sup <- subset(env.dat, select =c(Cellulose,
                                        Lignin,
                                        Carbohyrates_rel_abund,
                                        arom15_rel_abund,
                                        arom16_or_COO_rel_abund,
                                        Acids_rel_abund,
                                        Aliphatic_rel_abund,
                                        C,N,CN,d13C,d15N,BD,Mineral))

quanti.sup <- subset(env.dat, select =c(Cellulose,
                                        Lignin,
                                        Carbohyrates_rel_abund,
                                        arom15_rel_abund,
                                        arom16_or_COO_rel_abund,
                                        
                                        Aliphatic_rel_abund,
                                        C,CN,d15N,BD,Mineral))

q

ef <- envfit(res.pca ~ ., data=quanti.sup,perm=999)
#ef <- envfit(res.pca ~ ., data=env.dat,perm=999)

ef

quanti.coord <- cor(quanti.sup, res.pca$x)
quanti.cos2 <- quanti.coord^2




PCA_T<- fviz_pca_ind(res.pca,
                       axes=c(1,2),
                       label="none",pointsize=2,
                       habillage = dataset5$Treatment,
                       #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                       addEllipses = TRUE,
                       ellipse.type = "convex",
                       invisible = "quali")+
  # geom_text_repel(label=dataset5$Site,
  #                 nudge_x = 0,
  #                 nudge_y = -1,
  #                 size= 3.0,
  #                 box.padding = unit(0.1, "lines"), force = 2,
  #                 segment.colour = NA,
  #                 max.overlaps = 20)+
  ggtitle("")+
  # xlab("PC1 (79.8%)") +
  # ylab("PC2 (11.7%)")+
  theme_bw()+
  theme(legend.position = "right",
        legend.title=element_blank())+
  xlab("PC1 (57.0% of variance)")+
  ylab("PC2 (26.1% of variance)")



PCA_T<- fviz_add(PCA_T,quanti.coord*0.08,color = "blue",geom = "arrow",linetype = "solid",
         labelsize = 2.5,repel=TRUE)


PCA_T



###################################################################################

dataset5 <- subset(dataset4, Depth=="15-20 cm")
env.dat <- dataset5[,7:41] 
colnames(env.dat)[which(names(env.dat) == "carb")] <- "Carbohyrates_rel_abund"
colnames(env.dat)[which(names(env.dat) == "arom15")] <- "arom15_rel_abund"
colnames(env.dat)[which(names(env.dat) == "arom16")] <- "arom16_or_COO_rel_abund"
colnames(env.dat)[which(names(env.dat) == "acids")] <- "Acids_rel_abund"
colnames(env.dat)[which(names(env.dat) == "aliph28")] <- "Aliph28_PI"
colnames(env.dat)[which(names(env.dat) == "aliph29")] <- "Aliph29_PI"
colnames(env.dat)[which(names(env.dat) == "sum_arom")] <- "Sum_Aromatics_rel_abund"
colnames(env.dat)[which(names(env.dat) == "sum_aliph")] <- "Aliphatic_rel_abund"


res.pca <- prcomp(dataset5$spc2,center = TRUE,scale = FALSE)
quanti.sup <- subset(env.dat, select =c(Cellulose,
                                        Lignin,
                                        Carbohyrates_rel_abund,
                                        arom15_rel_abund,
                                        arom16_or_COO_rel_abund,
                                        Acids_rel_abund,
                                        Aliphatic_rel_abund,
                                        C,N,CN,d13C,d15N,BD,Mineral))

quanti.sup <- subset(env.dat, select =c(Cellulose,
                                        Lignin,
                                        Carbohyrates_rel_abund,
                                        arom15_rel_abund,
                                        arom16_or_COO_rel_abund,
                                        
                                        Aliphatic_rel_abund,
                                        C,N,CN,d15N,BD,Mineral))


ef <- envfit(res.pca ~ ., data=quanti.sup,perm=999)
#ef <- envfit(res.pca ~ ., data=env.dat,perm=999)
ef

quanti.coord <- cor(quanti.sup, res.pca$x)
quanti.cos2 <- quanti.coord^2




PCA_M<- fviz_pca_ind(res.pca,
                       axes=c(1,2),
                       label="none",pointsize=2,
                       habillage = dataset5$Treatment,
                       #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                       addEllipses = TRUE,
                       ellipse.type = "convex",
                       invisible = "quali")+
  # geom_text_repel(label=dataset5$Site,
  #                 nudge_x = 0,
  #                 nudge_y = -1,
  #                 size= 3.0,
  #                 box.padding = unit(0.1, "lines"), force = 2,
  #                 segment.colour = NA,
  #                 max.overlaps = 20)+
  ggtitle("")+
  # xlab("PC1 (79.8%)") +
  # ylab("PC2 (11.7%)")+
  theme_bw()+
  theme(legend.position = "right",
        legend.title=element_blank())+
  xlab("PC1 (63.2% of variance)")+
  ylab("PC2 (18.8% of variance)")

PCA_M <- fviz_add(PCA_M,quanti.coord*0.08,color = "blue",geom = "arrow",linetype = "solid",
         labelsize = 2.5,repel=TRUE)

PCA_M
###################################################################################

dataset5 <- subset(dataset4, Depth=="45-50 cm")
env.dat <- dataset5[,7:41] 
colnames(env.dat)[which(names(env.dat) == "carb")] <- "Carbohyrates_rel_abund"
colnames(env.dat)[which(names(env.dat) == "arom15")] <- "arom15_rel_abund"
colnames(env.dat)[which(names(env.dat) == "arom16")] <- "arom16_or_COO_rel_abund"
colnames(env.dat)[which(names(env.dat) == "acids")] <- "Acids_rel_abund"
colnames(env.dat)[which(names(env.dat) == "aliph28")] <- "Aliph28_PI"
colnames(env.dat)[which(names(env.dat) == "aliph29")] <- "Aliph29_PI"
colnames(env.dat)[which(names(env.dat) == "sum_arom")] <- "Sum_Aromatics_rel_abund"
colnames(env.dat)[which(names(env.dat) == "sum_aliph")] <- "Aliphatic_rel_abund"


res.pca <- prcomp(dataset5$spc2,center = TRUE,scale = FALSE)
quanti.sup <- subset(env.dat, select =c(Cellulose,
                                        Lignin,
                                        Carbohyrates_rel_abund,
                                        arom15_rel_abund,
                                        arom16_or_COO_rel_abund,
                                        Acids_rel_abund,
                                        Aliphatic_rel_abund,
                                        C,N,CN,d13C,d15N,BD,Mineral))

quanti.sup <- subset(env.dat, select =c(Cellulose,
                                        Lignin,
                                        Carbohyrates_rel_abund,
                                        arom15_rel_abund,
                                        arom16_or_COO_rel_abund,
                                        Acids_rel_abund,
                                        
                                        C,N,d13C,d15N,BD,Mineral))

ef <- envfit(res.pca ~ ., data=quanti.sup,perm=999)
#ef <- envfit(res.pca ~ ., data=env.dat,perm=999)
ef

quanti.coord <- cor(quanti.sup, res.pca$x)
quanti.cos2 <- quanti.coord^2



PCA_B<- fviz_pca_ind(res.pca,
                     axes=c(1,2),
                     label="none",pointsize=2,
                     habillage = dataset5$Treatment,
                     #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                     addEllipses = TRUE,
                     ellipse.type = "convex",
                     invisible = "quali")+
  # geom_text_repel(label=dataset5$Site,
  #                 nudge_x = 0,
  #                 nudge_y = -1,
  #                 size= 3.0,
  #                 box.padding = unit(0.1, "lines"), force = 2,
  #                 segment.colour = NA,
  #                 max.overlaps = 20)+
  ggtitle("")+
  # xlab("PC1 (79.8%)") +
  # ylab("PC2 (11.7%)")+
  theme_bw()+
  theme(legend.position = "right",
        legend.title=element_blank())+
  xlab("PC1 (63.2% of variance)")+
  ylab("PC2 (12.0% of variance)")

PCA_B <- fviz_add(PCA_B,quanti.coord*0.08,color = "blue",geom = "arrow",linetype = "solid",
                  labelsize = 2.5,repel=TRUE)

PCA_B



####################################################################################
####################################################################################
ggarrange(PCA_T,PCA_M,PCA_B,
          labels = c(" A) 0-5 cm", "B) 15-20 cm ", "C) 45-50 cm"),
          ncol = 3,
          common.legend = TRUE)

####################################################################################
#####################################################################################






res.pca <- PCA(dataset5$spc2,
               quanti.sup = env,graph=FALSE)

print(res.pca)
plot.PCA(dataset5$spc2,axes=c(1,2),choix = c("ind","var"))
plot.PCA(dataset5$spc2, axes = c(1,2), choix=c("ind", "var"))

plot(res.pca, choix = "var")
fviz_pca_var(res.pca)
fviz_pca_var(res.pca, alpha.var="contrib")+
  theme_minimal()
plot(res.pca, choix = "ind")
fviz_pca_ind(res.pca)
fviz_pca_ind(res.pca, geom="text")

fviz_pca_ind(res.pca, label="none",habillage=dataset5$Treatment)

plot(res.pca, choix="var")
plot(res.pca, choix="ind",scale.unit=TRUE)
fviz_pca_ind(res.pca)
fviz_pca_ind(result1)

env <- dataset5[,6:11]


?fviz_add

PCA_Trt


?fviz_pca
?fviz_

fviz_pca_var(result1,col.var="steelblue")+
  theme_minimal()

ef <- envfit(result1 ~ C+N+CN, data = env , perm =999)
ef

ef.adj <- ef
pvals.adj <-  p.adjust(ef$vectors$pvals, method ='bonferroni')
ef.adj$vectors$pvals <- pvals.adj
ef.adj

ordiplot(result2, display = 'sites')
plot(ef)



scores (result1, display = 'sites', scaling = 0)
scores (result2, display = 'sites', scaling = 0)


scores (RDA, display = 'sites', scaling = 0)
scores(RDA, choices=c(1,2))
?scores



ordiplot (result1, display = 'sites')
ordiplot (result2 , display = 'sites')


points (PCA, pch = vegtype, col = vegtype)

View(res.pca)

summary(res.pca)
ordiplot (res.pca)
ordiplot (RDA, choices = c(1,2), type = 'n')
points (res.pca, choices = c(1,2), display = 'sites', pch = as.character (dataset5$Treatment), col = dataset5$Treatment)

ef_ell <- envfit (RDA, env, choices = c(1,2), display = 'site', na.rm = TRUE)
ordiplot (RDA, choices = c(1,2), type = 'n')
points (RDA, choices = c(1,2), display = 'sites', pch = as.character (dataset5$Treatment), col = dataset5$Treatment)
plot (ef_ell, p.max = 0.05)  # to display only significant Ell. values


View(RDA)

RDA <- rda(dataset5$spc2~C+N,data=dataset5)

?envfit

env.all<- dataset5[,6:11]

ordiplot (tbRDA, display = 'species', choices = c(3,4), type = 'n')
ordiplot (RDA, choices = c(3,4), type = 'n')
1.3697 /2.6271
3.222 / 5.639   
vltava.env <- read.delim ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/data/vltava-env.txt')
vltava.spe <- read.delim ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/data/vltava-spe.txt', row.names = 1)

ef_ell <- envfit (RDA, dataset5[,6:11], choices = c(3,4), display = 'species', na.rm = TRUE)


str(dataset5$spc2)
str(dataset5[,6:11])
points (RDA, choices = c(3,4), display = 'sites', pch = as.character (dataset5$Treatment), col = dataset5$Treatment)
plot (ef_ell, p.max = 0.05)  # to display only significant Ell. values
ordiplot (RDA)
RDA

0.21150/1.3443 
head(summary(RDA))



min(dataset5$spc2)
res.pca <- prcomp(dataset5$spc2,center = TRUE, scale = FALSE)

Dpth_T=fviz_pca_ind(res.pca,
                    axes = c(1, 2),
                    label="none",pointsize=2,
                    habillage = dataset5$Treatment,
                    palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                    addEllipses = TRUE,
                    invisible = "quali",
                    ellipse.type = "convex"
)+
  geom_text_repel(label=dataset5$Site,
                  nudge_x = 0,
                  nudge_y = -0.1,
                  size= 3.0,
                  box.padding = unit(0.1, "lines"), force = 2,
                  segment.colour = NA,
                  max.overlaps = 20)+
  ggtitle("")+
  #xlab("PC1 (86.5%)") + 
  #ylab("PC2 (7.0%)")+
  theme_bw()+
  theme(legend.position = "top",
        legend.title=element_blank())
#   geom_path(data = as.data.frame(res.pca$x[,1:2])[a,],aes(x=PC1,y = PC2,group=dataset5[a,]$Site),
#           arrow = arrow(length=unit(0.15,"cm"), ends="last", type = "closed"))
# 
# 
# a=which(dataset5$Treatment == 'Undrained'|dataset5$Treatment == 'Drained')

Dpth_T
