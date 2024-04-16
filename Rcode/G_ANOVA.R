###2019 September#######
#rm(list = ls())
getwd()
EA.dat=read.csv("C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/MIR/Data/EA_data.csv")

##Explore data#######################
str(EA.dat)


###Orgonize dataframe###################
EA.dat$Site=as.factor(EA.dat$Site)
str(EA.dat)
##library file#######################################
library(nlme)
library(lsmeans)
library(multcomp)
library(geoR)
library(car)
library(ggplot2)
library(ggpubr)
library(tidyr)
########################################################################
New.dat=EA.dat %>% 
  gather(C:BD,key = parameter,value = y)
########################################################################
##Lable########
########################################################################


New.dat$Depth=factor(New.dat$Depth,levels = c("T","M","B"),
                    labels = c("0-5 cm","15-20 cm","45-50 cm"))

New.dat$Treatment=factor(New.dat$Treatment,levels = c("U","R","D"),
                        labels = c("Undrained","Rewetted","Drained"))


str(EA.dat)

########################################################################
#Figure1
#as POM difference in 2018 BR 
Dataset.dat= subset(New.dat, 
                      #New.dat$Site=="Breton"&
                      #New.dat$Date=="2018 Sep"&
                      #New.dat$Crop %in% c("Spring-grain","Perennial-grain")&
                      #New.dat$Depth=="0-5 cm"&
                      New.dat$parameter %in% c("C","N","CN","d13C","d15N","BD"))



ggboxplot(Dataset.dat,"Treatment","y",
          scales ="free",add="jitter",
          ncol=1,fill="Treatment",palette = c(rgb(0.1,0.1,0.7,0.8),"grey","grey", rgb(0.8,0.1,0.3,0.6), rgb(0.8,0.1,0.3,0.6)))+
  ylab(expression("POM-C Concentration ("*g~C~kg^-1~soil*")"))+
  theme(legend.position = "none")+
  facet_grid(Depth~parameter,scales = "free")


#+
  annotate(geom="text",x=1, y=10, label="a",color="black",fontface=3,size=5)+
  annotate(geom="text",x=2, y=11, label="ab",color="black",fontface=3,size=5)+
  annotate(geom="text",x=3, y=11, label="ab",color="black",fontface=3,size=5)+
  annotate(geom="text",x=4, y=12, label="b",color="black",fontface=3,size=5)+
  annotate(geom="text",x=5, y=12, label="b",color="black",fontface=3,size=5)
#text(x= 1, y= 3, labels= "some text")
text(x= 2, y= 3, labels= "more \n text")
text(x= 3, y= 3, labels= "red text", col= "red")
