dataset5 <- dataset4
dataset5 <- subset(dataset4, Depth=="0-5 cm")
dataset5 <- subset(dataset4, Depth=="15-20 cm")
dataset5 <- subset(dataset4, Depth=="45-50 cm")

res.pca <- prcomp(dataset5$spc2,center = TRUE, scale = FALSE)

Scores <- cbind(res.pca$x[,1:4],dataset5[,1:11])

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

ggline(Scores2,"Site","Scores",
       add = c("mean_se"),
       group = "Treatment",
       shape = "Treatment",
       linetype = "Treatment",
       color = "Treatment",
       palette = c("#00AFBB", "#E7B800","red"),
       facet.by = c("parameter"),scale="free")+
  theme(axis.text.x = element_text(angle = 45,hjust=1))


ggline(Scores2,"Site","Scores",
       add = c("mean_se"),
       group = "Depth",
       shape = "Depth",
       linetype = "Depth",
       color = "Depth",
       palette = c("#00AFBB", "#E7B800","red"),
       facet.by = c("parameter"),scale="free")+
  theme(axis.text.x = element_text(angle = 45,hjust=1))

  


###
Correlations=cbind(Dpth_M$data[,3:4],dataset5[,6:11])
colnames(Correlations)[1] <- "PC1"
colnames(Correlations)[2] <- "PC2"

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
