##Package#
library(ggcorrplot)
library(ggpmisc)
library(GGally)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(scales)
library(ggpubr)
library(tidyr)
#install.packages("vctrs")

library(grid)
#install.packages("ggpmisc")
#install.packages("glue")
## function for axis


########################################################################
str(Variables.dat)
Correlations <- subset(Variables.dat, select= -c(ID,Country,Site,Treatment,TsR))
Correlations <- subset(Variables.dat, 
                       select = c(carb,arom15,arom16,acids,aliph28,aliph28, CN,Lignin,sum_arom,HI_arom,aromatic_est
                                  ,Longitude, Latitude))
Correlations <- subset(Variables.dat, 
                       select = c(carb,arom15,arom16,acids,aliph28,aliph29, CN, d13C, d15N, sum_arom, sum_aliph
                                  ,Longitude, Latitude,Temperature, Precipitation))

###Simply regression to find correlation
cor <- ggpairs(Correlations)
#cor
#############################################################################
#aggregate#


Variables2.dat <-Variables.dat %>% 
  group_by(Depth,Treatment) %>% 
  summarize(carb=mean(carb),Cellulose=mean(Cellulose),
            arom15=mean(arom15),Lignin=mean(Lignin),
            arom16=mean(arom16),arom1516=mean(arom15+arom16)) %>% 
  mutate(Treatment=as.factor(Treatment))


Variables2.dat$Treatment=factor(Variables2.dat$Treatment,levels = c("D","R","U"),
                        labels = c("Drained","Rewetted","Undrained"))

##################################################################################
#linear regression for carbohydrates and celluse and arom15 and lignin
##################################################################################
A <-ggplot(Variables2.dat, aes(x = carb, y = Cellulose)) +
  geom_smooth(method="lm")+
  geom_point(aes(color = Treatment))+
  stat_poly_eq(formula = y ~ poly(x, 1, raw = TRUE), 
               aes(label = paste(..eq.label.., ..adj.rr.label.., ..p.value.label.., sep = "*`,`~")), 
               parse = TRUE,
               label.x.npc = "left",
               vstep = 0.05
               )+ # sets vertical spacing
  theme_bw()+
  #scale_x_continuous(label=scientific_10)+
  scale_x_log10(breaks = trans_breaks('log10', function(x) 10^x),
                labels = trans_format('log10', math_format(10^.x)))+
  xlab(expression("MIR"[carb]*"/"*"area"))+
  ylab("Cellulose (%)")+
  theme(legend.position = "top",
        legend.title=element_blank())
A
######################################################


B <-ggplot(Variables2.dat, aes(x = arom15, y = Lignin)) +
  geom_smooth(method="lm")+
  geom_point(aes(color = Treatment))+
  stat_poly_eq(formula = y ~ poly(x, 1, raw = TRUE), 
               aes(label = paste(..eq.label.., ..adj.rr.label.., ..p.value.label.., sep = "*`,`~")), 
               parse = TRUE,
               label.x.npc = "left",
               vstep = 0.05)+ # sets vertical spacing
  theme_bw()+
  #scale_x_continuous(label=scientific_10)+
  scale_x_log10(breaks = trans_breaks('log10', function(x) 10^x),
                labels = trans_format('log10', math_format(10^.x)))+
  xlab(expression("MIR"[arom15]*"/"*"area"))+
  ylab("Lignin (%)")+
  theme(legend.position = "top",
        legend.title=element_blank())

B  
##arrangemment of figures 
ggarrange(A,B,labels = c("A","B"),
          common.legend = TRUE,legend = "top",ncol=2)

##################################################################################
##################################################################################


str(Variables.dat)

Quan_variable.dat = subset(Variables.dat,select = c(Site,Treatment,Depth,carb,arom15,arom16,sum_aliph,acids))

Quan_variable2.dat=Quan_variable.dat %>% 
  gather(carb:acids, key=Compounds, value = value)

Quan_variable2.dat$Depth=factor(Quan_variable2.dat$Depth,levels = c("T","M","B"),
                      labels = c("0-5 cm","15-20 cm","45-50 cm"))

Quan_variable2.dat$Treatment=factor(Quan_variable2.dat$Treatment,levels = c("D","R","U"),
                          labels = c("Drained","Rewetted","Undrained"))

# Quan_variable2.dat$Compounds=factor(Quan_variable2.dat$Compounds,levels = c("carb","sum_arom","sum_aliph","acids"),
#                     labels = c("Carbohydrates","Aromatics","Aliphatic","Acids"))

Quan_variable2.dat$Compounds=factor(Quan_variable2.dat$Compounds,levels = c("carb","arom15","arom16",
                                                                            "acids","sum_aliph"),
                                    labels = c("Carbohydrates","Aromatics(arom15)",expression("Aromatics(arom16)~or~COO^{-1}"),"Acids(COOH)","Aliphatics"))



ggboxplot(Quan_variable2.dat, x="Depth", "value", shape="Treatment",color="Treatment",
          #add = c("jitter"),
          )+
  xlab("")+
  ylab(expression("MIR"[intensity]*"/"*"area"))+
  theme_bw()+
  theme(legend.position = "top",
        legend.title=element_blank())+
  facet_wrap(~Compounds, labeller = label_parsed,scales="free")

########################################################################################
##d15N vs carbohydrates## aggregate
str(Variables.dat)

Variables2.dat <-Variables.dat %>% 
  group_by(Depth,Treatment,Site) %>% 
  summarize(carb=mean(carb),d15N=mean(d15N)) %>% 
  mutate(Treatment=as.factor(Treatment))


Variables2.dat$Treatment=factor(Variables2.dat$Treatment,levels = c("D","R","U"),
                                labels = c("Drained","Rewetted","Undrained"))

##################################################################################
#linear regression for carbohydrates and celluse and carbohydrates and lignin
##################################################################################
plot <-ggplot(Variables2.dat, aes(x = d15N, y = carb)) +
  geom_smooth(method="lm")+
  geom_point(aes(color = Treatment))+
  stat_poly_eq(formula = y ~ poly(x, 1, raw = TRUE), 
               aes(label = paste(..eq.label.., ..adj.rr.label.., ..p.value.label.., sep = "*`,`~")), 
               parse = TRUE,
               label.x.npc = "left",
               vstep = 0.05
  )+ # sets vertical spacing
  theme_bw()+
  #scale_x_continuous(label=scientific_10)+
  scale_x_log10(breaks = trans_breaks('log10', function(x) 10^x),
                labels = trans_format('log10', math_format(10^.x)))+
  #xlab(expression("MIR"[carb]*"/"*"area"))+
  #ylab("Cellulose (%)")+
  theme(legend.position = "top",
        legend.title=element_blank())
plot









res.pca <- prcomp(Quan_variable.dat[,4:7],center = TRUE,scale = TRUE)

fviz_pca_biplot(res.pca,
                axes=c(1,2),
                label="none",pointsize=2,
                habillage = Quan_variable.dat$Treatment,
                #palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                addEllipses = TRUE,
                ellipse.type = "convex",
                invisible = "quali")



ggpar(p, palette = "jco")

library(tidyr)
library(dplyr)
library(plyr)
detach(package:plyr)

str(Variables.dat)
New.test <- Variables.dat %>% 
  group_by(Site,Treatment) %>% 
  summarise(carb.ave = mean(carb),
            aliph28.ave = mean(aliph28),
            aliph29.ave = mean(aliph29),
            arom15.ave = mean(arom15),
            arom16.ave = mean(arom16),
            temp.avg = mean(Temperature),
            Precip.avg = mean(Precipitation),
            Longitude.avg = mean(Longitude))

New.test

A <- ggpairs(New.test)
A
  ggplot(New.test, aes(x = temp.avg, y = aliph28.ave)) +
  geom_smooth(method="lm")+
  geom_point()+
  stat_poly_eq(formula = y ~ poly(x, 1, raw = TRUE), 
               aes(label = paste(..eq.label.., ..adj.rr.label.., ..p.value.label.., sep = "*`,`~")), 
               parse = TRUE,
               label.x.npc = "left",
               vstep = 0.05)+ # sets vertical spacing
  theme_bw()+
  scale_x_continuous(label=scientific_10)
