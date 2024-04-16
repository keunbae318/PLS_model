############################################################
#####Elemental analysis for graph
############################################################

getwd()
EA.dat=read.csv("C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/Elemental_Analysis/Rdata/EA_data.csv")

##Explore data#######################
str(EA.dat)


###Orgonize dataframe###################
EA.dat$Site=as.factor(EA.dat$Site)

EA.dat$Depth=factor(EA.dat$Depth,levels = c("Top","Middle","Bottom"),
                    labels = c("0-5 cm","15-20 cm","45-50 cm"))

EA.dat$Treatment=factor(EA.dat$Treatment,levels = c("Undrained","Rewetted","Drained"),
                        labels = c("Undrained","Rewetted","Drained"))

levels(EA.dat$Treatment)
table(EA.dat$Treatment)
##library file#######################################
library(nlme)
#install.packages("lsmeans")
library(lsmeans)
#install.packages("geoR")
library(multcomp)
library(geoR)
library(car)
library(dplyr)
library(tidyr)
library(ggpubr)
######################################################

str(EA.dat)

format.dat <- EA.dat %>% 
  gather(C:BD,key = parameter, value = y)
#############################################################
##Lable###

format.dat$parameter=factor(format.dat$parameter,levels = c("C","N","CN","d13C","d15N","BD"),
                            labels = c(expression("C(%)"),
                                       expression("N(%)"),
                                       expression("C:N"),
                                       expression("d13C"),
                                       expression("d15N"),
                                       expression("BD")))

ggbarplot(format.dat,x="Treatment",y="y",color="Treatment",facet.by = c("parameter","Depth"),
          add="mean_se",scales="free")#+
  # facet(scales = "free_y")#+
  # facet_grid(scale="free",labeller = label_parsed)


merged_data <- EA.dat
merged_data$C <- merged_data$C*10

custom_labels <- labeller(parameter = c("C" = "Peat C concentration (mg/g)", 
                                        "Shannon" = "Shannon diversity index (H)"))

ggboxplot(merged_data, "Treatment", "C",fill = "Treatment", facet.by = "Depth")+
  scale_fill_manual(values=c("#00AFBB","#E7B800","#FC4E07"))+
  xlab(paste0(""))+
  ylab(paste0("Peat C concentration (mg/g)"))+
  theme(legend.position = "top",
        legend.title=element_blank())
