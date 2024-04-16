dataset5 <- subset(dataset4, Treatment=="Undrained")
dataset5 <- subset(dataset4, Treatment=="Rewetted")
dataset5 <- subset(dataset4, Treatment=="Drained")

#####################
#Statstical test
#####################
res.pca <- prcomp(dataset5$spc5,center = TRUE,scale = FALSE)
#Treatment#
adonis2(dataset5$spc5~Treatment,dataset5, 
        method = "euclidean",
        strata = dataset5$Site)
pairwise.adonis2(dataset5$spc5 ~ Treatment, 
                 data = dataset5, 
                 strata = 'Site',
                 method = "euclidean")
#Depth#
adonis2(dataset5$spc5~Depth,dataset5,
        method = "euclidean",
        strata = dataset5$Site)
pairwise.adonis2(dataset5$spc5 ~ Depth, 
                 data = dataset5, 
                 strata = 'Site',
                 method = "euclidean")

##############################################################################
##############################################################################
#All dataset of Depth n Treatment
dataset5 <- dataset4
dataset5 <- subset(dataset4, Treatment=="Undrained")
dataset5 <- subset(dataset4, Treatment=="Rewetted")
dataset5 <- subset(dataset4, Treatment=="Drained")
res.pca <- prcomp(dataset5$spc5,center = TRUE, scale = FALSE)




#Depth

PCA_Dpth<- fviz_pca_ind(res.pca,
                        label="none",pointsize=2,
                        habillage = dataset5$Depth,
                        palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                        addEllipses = FALSE,
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
        legend.title=element_blank())
PCA_Dpth

#Site
PCA_Site<- fviz_pca_ind(res.pca,
                        label="none",pointsize=2,
                        habillage = dataset5$Site,
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
        legend.title=element_blank())
PCA_Site


####################################################################################
####################################################################################
ggarrange(PCA_Trt,PCA_Dpth,
          labels = c("A)Regime","B) Depth "),
          ncol = 1,widths = c(0.5, 1),
          common.legend = FALSE)


####################################################################################
####################################################################################

dataset5 <- subset(dataset4, Depth=="0-5 cm")
res.pca <- prcomp(dataset5$spc5,center = TRUE, scale = FALSE)

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

###
Correlations=cbind(Dpth_B$data[,3:4],dataset5[,6:11])



#################################################################################################

dataset5 <- subset(dataset4, Depth=="15-20 cm")
#dataset5 <- subset(dataset5,!Site %in% ("Gutzkow"))
res.pca <- prcomp(dataset5$spc5,center = TRUE, scale = FALSE)

Dpth_M=fviz_pca_ind(res.pca,
                    axes = c(1, 2),
                    label="none",pointsize=2,
                    habillage = dataset5$Treatment,
                    palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                    addEllipses = TRUE,
                    ellipse.type = "convex",
                    invisible = "quali"
)+
  geom_text_repel(label=dataset5$Site,
                  nudge_x = 0,
                  nudge_y = -0.1,
                  size= 3.0,
                  box.padding = unit(0.1, "lines"), force = 2,
                  segment.colour = NA,
                  max.overlaps = 20)+
  ggtitle("")+
  # xlab("PC1 (79.7%)") + 
  # ylab("PC2 (10.9%)")+
  theme_bw()+
  theme(legend.position = "top",
        legend.title=element_blank())

Dpth_M


# Contributions of variables to PC1
fviz_contrib(res.pca, choice = "var", axes = 1:2, top = 100)


fviz_pca_biplot(res.pca,repel = TRUE,select.var = list(contrib=100))


Dpth_M_biplot=fviz_pca_biplot(res.pca,select.var = list(contrib=50),
                              axes = c(1, 2),
                              label="none",pointsize=2,
                              habillage = dataset5$Treatment,
                              palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                              addEllipses = TRUE,
                              ellipse.type = "convex",
                              invisible = "quali"
)+
  geom_text_repel(label=dataset5$Site,
                  nudge_x = 0,
                  nudge_y = -0.1,
                  size= 3.0,
                  box.padding = unit(0.1, "lines"), force = 2,
                  segment.colour = NA,
                  max.overlaps = 20)+
  ggtitle("")+
  # xlab("PC1 (79.7%)") + 
  # ylab("PC2 (10.9%)")+
  theme_bw()+
  theme(legend.position = "top",
        legend.title=element_blank())


Dpth_M_biplot
####################################################################################

dataset5 <- subset(dataset4, Depth=="45-50 cm")
#dataset5 <- subset(dataset5,!Site %in% c("Arlon","Binnenveld","Gutzkow","Zwarte Beek"))
res.pca <- prcomp(dataset5$spc5,center = TRUE, scale = FALSE)

Dpth_B=fviz_pca_ind(res.pca,
                    label="none",pointsize=2,
                    habillage = dataset5$Treatment,
                    palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                    addEllipses = TRUE,
                    ellipse.type = "convex",
                    invisible = "quali"
)+
  geom_text_repel(label=dataset5$Site,
                  nudge_x = 0,
                  nudge_y = -0.05,
                  size= 3.0,
                  box.padding = unit(0.1, "lines"), force = 2,
                  segment.colour = NA,
                  max.overlaps = 20)+
  ggtitle("")+
  # xlab("PC1 (79.8%)") +
  # ylab("PC2 (11.7%)")+
  theme_bw()+
  theme(legend.position = "top",
        legend.title=element_blank())
Dpth_B

####################################################################################
####################################################################################
ggarrange(Dpth_T,Dpth_M,Dpth_B,
          labels = c(" A) 0-5 cm", "B) 15-20 cm ", "C) 45-50 cm"),
          ncol = 3,
          common.legend = TRUE)

####################################################################################
#####################################################################################

#Loading graph#########
dataset5 <- dataset4

#####################################
###PC_loading_graph_0-5 cm
#####################################
dataset5 <- subset(dataset4, Depth=="0-5 cm")
res.pca <- prcomp(dataset5$spc5,center = TRUE, scale = FALSE)


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
#####################################
###PC_loading_graph_15-20 cm
#####################################
dataset5 <- subset(dataset4, Depth=="15-20 cm")
res.pca <- prcomp(dataset5$spc5,center = TRUE, scale = FALSE)


loadings_data <- as_tibble(res.pca$rotation) %>% 
  mutate(INDEX = row.names(res.pca$rotation)) %>% 
  relocate(INDEX)
loadings_data$INDEX <- as.numeric(loadings_data$INDEX)


loading.spc <- as_tibble(loadings_data[,1:3]) %>% 
  gather(PC1:PC2,
         key = PC_axis,
         value = Loadings)

loading.spc$INDEX <- as.numeric(loading.spc$INDEX)

Loading2=ggplot(loading.spc,aes(x=INDEX,y=Loadings,group=PC_axis,colour=PC_axis))+ 
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

Loading2
####################################################################################
####################################################################################
ggarrange(Loading1,Loading2,
          labels = c("   A) 0-5 cm", " B) 15-20 cm "),
          ncol = 1,
          common.legend = TRUE)

###Loading graph for explanatory variables

fviz_eig(res.pca, addlabels = TRUE)

######################################################################################

#####################################
dataset5 <- subset(dataset4, Depth=="45-50 cm")
res.pca <- prcomp(dataset5$spc5,center = TRUE, scale = FALSE)


loadings_data <- as_tibble(res.pca$rotation) %>% 
  mutate(INDEX = row.names(res.pca$rotation)) %>% 
  relocate(INDEX)
loadings_data$INDEX <- as.numeric(loadings_data$INDEX)


loading.spc <- as_tibble(loadings_data[,1:4]) %>% 
  gather(PC1:PC3,
         key = PC_axis,
         value = Loadings)

loading.spc$INDEX <- as.numeric(loading.spc$INDEX)


Loading3=ggplot(loading.spc,aes(x=INDEX,y=Loadings,group=PC_axis,colour=PC_axis))+ 
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
  scale_x_reverse(limit=c(4000,500),expand=c(0,0))+
  grids(linetype = "dashed")+
  #scale_color_manual(values=c("#0073C2FF", "#EFC000FF"))+
  stat_peaks(span = 21, geom = "point", colour = "black")+
  stat_valleys(span = 21, geom = "point", colour = "black")+
  stat_peaks(span = 21, geom = "text", colour = "black", vjust = 0.2, hjust = -0.2,angle=90,
             label.fmt = "%3.0f"
  )+
  stat_valleys(span = 21, geom = "text", colour = "black", vjust = 0.2, hjust = -0.2,angle=270,
               label.fmt = "%3.0f")+
  ylim(-0.28,0.24)
Loading3
####################################################################################




