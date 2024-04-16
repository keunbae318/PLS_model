library(ggcorrplot)
#install.packages("ggpmisc")
library(ggpmisc)
#install.packages("GGally")
library(GGally)
######################################################################

theme<-theme(axis.line = element_line(colour = "black"),
             axis.text=element_text(colour="black",size=12),
             axis.title=element_text(size=12),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank(),
             panel.border = element_rect(colour = "black",fill=NA),
             panel.background = element_blank(),
             legend.position = c(1, 1), 
             legend.justification = c(1,1), 
             legend.title=element_blank(),
             legend.background = element_rect(colour = NA, fill = NA),
             legend.key=element_rect(colour = NA,fill = "white"),
             legend.key.size = unit(0.020, "npc"),
             legend.text=element_text(size=12),
             plot.margin=unit(c(0.2,1,0.1,0.2),"cm"))
######################################################################
dataset5 <- dataset4
dataset5 <- subset(dataset4, Depth=="0-5 cm")
res.pca <- prcomp(dataset5$spc5,center = TRUE, scale = TRUE)

Dpth_T=fviz_pca_ind(res.pca,
                    axes = c(3, 4),
                    label="none",pointsize=2,
                    habillage = dataset5$Treatment,
                    palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                    addEllipses = TRUE,
                    invisible = "quali",
                    ellipse.type = "convex"
)+
  geom_text_repel(label=dataset5$Site,
                  nudge_x = 0,
                  nudge_y = -1,
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
# geom_path(data = as.data.frame(res.pca$x[,1:2])[a,],aes(x=PC1,y = PC2,group=dataset5[a,]$Site),
#         arrow = arrow(length=unit(0.15,"cm"), ends="last", type = "closed"))


#a=which(dataset5$Treatment == 'Undrained'|dataset5$Treatment == 'Drained')

Dpth_T$data

###
dataset5 <- subset(dataset4, Depth=="15-20 cm")
res.pca <- prcomp(dataset5$spc5,center = TRUE, scale = FALSE)


Correlations=cbind(res.pca$x[,1:3],dataset5[,6:11])
colnames(Correlations)[1] <- "PC1"
colnames(Correlations)[2] <- "PC2"
colnames(Correlations)[3] <- "PC3"

?ggpairs
A <- ggpairs(Correlations)
A
A <- ggpairs(Correlations)

Scores <- cbind(res.pca$x[,1:3],dataset5[,1:11])

write.csv(Correlations,"C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/MIR/Data/PC.csv")
corr <- round(cor(Correlations),1)
p.mat <- cor_pmat(Correlations)

######################################

Cor_T <-  ggcorrplot(corr,
                     hc.order = FALSE,
                     type="upper",
                     #outline.color = "White",
                     #method = "circle",
                     lab=TRUE,
                     insig = "blank",
                     sig.level = 0.05,
                     p.mat = p.mat)

Cor_T

# linear_fig_1=ggplot(Correlations,aes(x=PC1,y=BD))+
#   theme+
#   geom_point()+
#   geom_smooth(method = "lm")+
#   stat_poly_eq(formula = y~x,
#                aes(label=paste(..eq.label..,..rr.label..,..p.value.label..,sep="*','~")),
#                parse=TRUE,
#                label.x.npc = 0.15,label.y.npc = 0.75,size=5.5)+
#   ylab(expression("Bulk Density ("*g~cm^-1*")"))+
#   xlab(expression("PC1"))
# 
# linear_fig_1



dataset5 <- subset(dataset4, Depth=="15-20 cm")
res.pca <- prcomp(dataset5$spc5,center = TRUE, scale = TRUE)

Dpth_M <- fviz_pca_ind(res.pca,
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
                  nudge_y = -1,
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
# geom_path(data = as.data.frame(res.pca$x[,1:2])[a,],aes(x=PC1,y = PC2,group=dataset5[a,]$Site),
#         arrow = arrow(length=unit(0.15,"cm"), ends="last", type = "closed"))


#a=which(dataset5$Treatment == 'Undrained'|dataset5$Treatment == 'Drained')

Dpth_M


###
Correlations=cbind(Dpth_M$data[,3:4],dataset5[,6:11])
colnames(Correlations)[1] <- "PC1"
colnames(Correlations)[2] <- "PC2"


corr <- round(cor(Correlations),1)
p.mat <- cor_pmat(Correlations)



######################################

Cor_M <-  ggcorrplot(corr,
                     hc.order = FALSE,
                     type="upper",
                     #outline.color = "White",
                     #method = "circle",
                     lab=TRUE,
                     insig = "blank",
                     sig.level = 0.05,
                     p.mat = p.mat)

Cor_M


#############################################################################################
dataset5 <- subset(dataset4, Depth=="45-50 cm")
res.pca <- prcomp(dataset5$spc5,center = TRUE, scale = TRUE)

Dpth_B <- fviz_pca_ind(res.pca,
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
                  nudge_y = -1,
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
# geom_path(data = as.data.frame(res.pca$x[,1:2])[a,],aes(x=PC1,y = PC2,group=dataset5[a,]$Site),
#         arrow = arrow(length=unit(0.15,"cm"), ends="last", type = "closed"))


#a=which(dataset5$Treatment == 'Undrained'|dataset5$Treatment == 'Drained')

Dpth_B

###
Correlations=cbind(Dpth_B$data[,3:4],dataset5[,6:11])
colnames(Correlations)[1] <- "PC1"
colnames(Correlations)[2] <- "PC2"
write.csv(Correlations,"C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/MIR/Data/PC.csv")

corr <- round(cor(Correlations),1)
p.mat <- cor_pmat(Correlations)

######################################

Cor_B <-  ggcorrplot(corr,
                     hc.order = FALSE,
                     type="upper",
                     #outline.color = "White",
                     #method = "circle",
                     lab=TRUE,
                     insig = "blank",
                     sig.level = 0.05,
                     p.mat = p.mat)

Cor_B


####################################################################################
####################################################################################
ggarrange(Cor_T,Cor_M,Cor_B,
          labels = c(" A) 0-5 cm", "B) 15-20 cm ", "C) 45-50 cm"),
          ncol = 3,
          common.legend = TRUE)

####################################################################################
#####################################################################################

#Detail


Correlations=cbind(Dpth_B$data[,3:4],dataset5[,6:11])
colnames(Correlations)[1] <- "PC1"
colnames(Correlations)[2] <- "PC2"
ggpairs(Correlations)

A <- ggpairs(Correlations)
A <- ggpairs(Correlations)
