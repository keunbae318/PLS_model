# Exemple taken from the vignette of the chemometric package
library(chemometrics)

library(mvoutlier)
dataset5 <- dataset4
Test <- dataset4
Test <- subset(dataset5, Depth=="0-5 cm")
Test <- subset(dataset5, Depth=="15-20 cm")
Test <- subset(dataTestset4, Depth=="45-50 cm")
Test <- dataset5


res2 <- pcout(x=Test$spc5,makeplot = T)
outliers <-as.character(Test$ID[res2$x.dist1>9])
outliers

Test<- subset(Test,!ID %in% (outliers))
Test<- subset(Test,!ID %in% c("Zwarte Beek_R_T",
                             "Drentse Aa_U_T",
                             "Drentse Aa_D_T",
                             "Drentse Aa_R_T",
                             "Kiel_D_T",
                             "Suwalszczyzna_D_T"))


Test <- subset(Test, Depth=="0-5 cm")
Test <- subset(Test, Depth=="15-20 cm")
Test <- subset(dataTestset4, Depth=="45-50 cm")

#Test<- subset(Test,Test$Site == "Drentse Aa"))


res.pca <- prcomp(Test$spc3,center = TRUE, scale = FALSE)

Dpth_T=fviz_pca_ind(res.pca,
                    axes = c(1, 2),
                    label="none",pointsize=2,
                    habillage = Test$Treatment,
                    palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                    addEllipses = TRUE,
                    invisible = "quali",
                    ellipse.type = "convex"
)+
  geom_text_repel(label=Test$Site,
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

Dpth_T=fviz_pca_ind(res.pca,
                    axes = c(1, 2),
                    label="none",pointsize=2,
                    habillage = Test$Site,
                    #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                    addEllipses = TRUE,
                    invisible = "quali",
                    ellipse.type = "convex"
)+
  geom_text_repel(label=Test$Site,
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



Scores <- cbind(res.pca$x[,1:3],Test[,1:11])

Scores2 <-  Scores %>%  
  gather(PC1:PC3,
         key = parameter,
         value = Scores)

Scores$Site=reorder(Scores$Site,Scores$PC1,mean)
ggline(Scores2,"Site","Scores",
       #add = c("mean_se"),
       group = "Treatment",
       shape = "Treatment",
       linetype = "Treatment",
       color = "Treatment",
       palette = c("#00AFBB", "#E7B800","red"),
       facet.by = c("parameter","Depth"),scale="free")+
  theme(axis.text.x = element_text(angle = 45,hjust=1))





#######################################################################
########################################################################
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

Loading1
