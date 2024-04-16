getwd()
EA.dat=read.csv("C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/MIR/Data/EA_data.csv")

##Explore data#######################
str(EA.dat)


###Orgonize dataframe###################
EA.dat$Site=as.factor(EA.dat$Site)

EA.dat$Depth=factor(EA.dat$Depth,levels = c("T","M","B"),
                      labels = c("0-5 cm","15-20 cm","45-50 cm"))

EA.dat$Treatment=factor(EA.dat$Treatment,levels = c("U","R","D"),
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
###############################################################
###############################################################
#View(Dataset.dat)
##Data##
######################POM C##########################                 
ANOVA.dat <- EA.dat
#######
N <- N
########
ANOVA.lme=lme(BD~Treatment*Depth,random =~1|Site,na.omit(ANOVA.dat))
#####################################################
##Assumption##
plot(ANOVA.lme) #homoscedasticity
qqnorm(ANOVA.lme,~resid(.,type="p")|Site,abline = c(0,1)) #normality of resid
qqnorm(ANOVA.lme,~ranef(.,standard=T),abline=c(0,1)) # ind. of random effects
plot(ANOVA.lme,Site~resid(.)) #normality of ranef
shapiro.test(resid(ANOVA.lme))
shapiro.test(ranef(ANOVA.lme)[,1])
########################################################
summary(ANOVA.lme)
anova(ANOVA.lme)

#If not interaction
########################################################
ANOVA.lme=lme(BD~Treatment+Depth,random =~1|Site,na.omit(ANOVA.dat))

summary(ANOVA.lme)
anova(ANOVA.lme)

cld(lsmeans(ANOVA.lme,~Treatment))
cld(lsmeans(ANOVA.lme,~Depth))

###If assumptions are violated, the data need to be transformd with Boxcoxfit##
boxcoxfit(na.omit(ANOVA.dat$BD))
ANOVA.bcx=bcPower(ANOVA.dat$BD,lambda =0.03288139  )

ANOVA.dat$ANOVA.bcx=ANOVA.bcx
ANOVA.lme=lme(ANOVA.bcx~Treatment*Depth,random =~1|Site,na.omit(ANOVA.dat))
#Output#######################
#############################
summary(ANOVA.lme)
anova(ANOVA.lme)
#if interaction not significant###
ANOVA.lme=lme(ANOVA.bcx~Treatment+Depth,random =~1|Site,na.omit(ANOVA.dat))


summary(ANOVA.lme)
anova(ANOVA.lme)

cld(lsmeans(ANOVA.lme,~Treatment))
cld(lsmeans(ANOVA.lme,~Depth))
#########################################################
Table1=cld(lsmeans(ANOVA.lme,~Treatment))
Table2=cld(lsmeans(ANOVA.lme,~Depth))

Table_C=as.data.frame(Table1)
Table_F=as.data.frame(Table2)

Table_C <- Table_C[c("1", "2", "3"), ]
Table_F <- Table_F[c("1", "2", "3"), ]

names(Table_C)[names(Table_C)=="Treatment"] <-"Treatment"
names(Table_F)[names(Table_F)=="Depth"] <-"Treatment"

Table=rbind(Table_C,Table_F)
Table$C_value=10^(log10(Table$lsmean*0.03288139 +1)/0.03288139)
write.csv(Table,"R:/Guillermo's Lab/POM/Research/POM/POM Data/Data statical analysis/Sequence_POM/Re_analysis_table/Table.csv")
Table

