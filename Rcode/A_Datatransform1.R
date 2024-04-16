#MIR spectroscopy

# The workspace is cleared
rm(list=ls()) # The workspace is cleared
graphics.off()

#Packages
library(prospectr)
library(dplyr)
library(readr)
library(purrr)
library(readr)
library(tidyverse)
library(tibble)

#source("fun/misc.R") # source the file
#Workload
getwd()
setwd("C:/Users/u0152843/OneDrive - KU Leuven/My Drive/KULeuven_PhD/C1_RE_PEAT/Mid_Infrared/Data/Rawdata/")


################################################################################
#Data combination process with formatting.
################################################################################
##Plate1###
data_join1 <- list.files(path = "1/exp/", # Identify all CSV files
                        pattern = "*.csv", full.names = TRUE) %>% 
  lapply(read_csv) %>%                              # Store all files in list
  reduce(full_join, by = "Created as New Dataset")                      # Full-join data sets into one data set 


b<-t(data_join1[-1,])
colnames(b) <-round(as.numeric(b[1, ]))  
c<-b[-1,]
Plate_1 <- rowid_to_column(as.data.frame(c), var = "ID")
Plate_1[,1]<-rownames(c)
Plate_1$Plate  <- "P1"
Plate_1$New_ID  <-paste(Plate_1$Plate,"_",as.numeric(gsub("([0-9]+).*$", "\\1", Plate_1$ID)),sep="")
Plate_1$Sample_NO <-as.numeric(gsub("([0-9]+).*$", "\\1", Plate_1$ID))
Plate_1<-Plate_1[order(Plate_1$Sample_NO),]
rownames(Plate_1)=Plate_1$Sample_NO

##Plate2###
data_join2 <- list.files(path = "2/exp/", # Identify all CSV files
                         pattern = "*.csv", full.names = TRUE) %>% 
  lapply(read_csv) %>%                              # Store all files in list
  reduce(full_join, by = "Created as New Dataset")                      # Full-join data sets into one data set 


b<-t(data_join2[-1,])
colnames(b) <-round(as.numeric(b[1, ]))  
c<-b[-1,]
Plate_2 <- rowid_to_column(as.data.frame(c), var = "ID")
Plate_2[,1]<-rownames(c)
Plate_2$Plate  <- "P2"
Plate_2$New_ID  <-paste(Plate_2$Plate,"_",as.numeric(gsub("([0-9]+).*$", "\\1", Plate_2$ID)),sep="")
Plate_2$Sample_NO <-as.numeric(gsub("([0-9]+).*$", "\\1", Plate_2$ID))
Plate_2<-Plate_2[order(Plate_2$Sample_NO),]
rownames(Plate_2)=Plate_2$Sample_NO


##Plate3###
data_join3 <- list.files(path = "3/exp/", # Identify all CSV files
                         pattern = "*.csv", full.names = TRUE) %>% 
  lapply(read_csv) %>%                              # Store all files in list
  reduce(full_join, by = "Created as New Dataset")                      # Full-join data sets into one data set 


b<-t(data_join3[-1,])
colnames(b) <-round(as.numeric(b[1, ]))  
c<-b[-1,]
Plate_3 <- rowid_to_column(as.data.frame(c), var = "ID")
Plate_3[,1]<-rownames(c)
Plate_3$Plate  <- "P3"
Plate_3$New_ID  <-paste(Plate_3$Plate,"_",as.numeric(gsub("([0-9]+).*$", "\\1", Plate_3$ID)),sep="")
Plate_3$Sample_NO <-as.numeric(gsub("([0-9]+).*$", "\\1", Plate_3$ID))
Plate_3<-Plate_3[order(Plate_3$Sample_NO),]
rownames(Plate_3)=Plate_3$Sample_NO


##Plate4###
data_join4 <- list.files(path = "4/exp/", # Identify all CSV files
                         pattern = "*.csv", full.names = TRUE) %>% 
  lapply(read_csv) %>%                              # Store all files in list
  reduce(full_join, by = "Created as New Dataset")                      # Full-join data sets into one data set 


b<-t(data_join4[-1,])
colnames(b) <-round(as.numeric(b[1, ]))  
c<-b[-1,]
Plate_4 <- rowid_to_column(as.data.frame(c), var = "ID")
Plate_4[,1]<-rownames(c)
Plate_4$Plate  <- "P4"
Plate_4$New_ID  <-paste(Plate_4$Plate,"_",as.numeric(gsub("([0-9]+).*$", "\\1", Plate_4$ID)),sep="")
Plate_4$Sample_NO <-as.numeric(gsub("([0-9]+).*$", "\\1", Plate_4$ID))
Plate_4<-Plate_4[order(Plate_4$Sample_NO),]
rownames(Plate_4)=Plate_4$Sample_NO


##Plate5###
data_join5 <- list.files(path = "5/exp/", # Identify all CSV files
                         pattern = "*.csv", full.names = TRUE) %>% 
  lapply(read_csv) %>%                              # Store all files in list
  reduce(full_join, by = "Created as New Dataset")                      # Full-join data sets into one data set 


b<-t(data_join5[-1,])
colnames(b) <-round(as.numeric(b[1, ]))  
c<-b[-1,]
Plate_5 <- rowid_to_column(as.data.frame(c), var = "ID")
Plate_5[,1]<-rownames(c)
Plate_5$Plate  <- "P5"
Plate_5$New_ID  <-paste(Plate_5$Plate,"_",as.numeric(gsub("([0-9]+).*$", "\\1", Plate_5$ID)),sep="")
Plate_5$Sample_NO <-as.numeric(gsub("([0-9]+).*$", "\\1", Plate_5$ID))
Plate_5<-Plate_5[order(Plate_5$Sample_NO),]
rownames(Plate_5)=Plate_5$Sample_NO

##################################################################################################
#Data rbind with subsetting##
dataset=rbind(Plate_1,Plate_2,Plate_3,Plate_4,Plate_5)
dataset=dataset[1:468,]

######################################
#Construction for data frame
######################################

x=1:4
y=c("Top","Middle","Bottom")
z=c("Arlon","Zwarte Beek","Binnenveld","Drentse Aa","Gutzkow","Peene mouth","Kiel","Recknitz","Biebrza",
    "Rospuda","Suwalszczyzna","Mazury","Anglesey")
v=c("Undrained","Drained","Rewetted")

dataset$Rep=rep(x)
dataset$Rep=as.character(dataset$Rep)
dataset$Depth=rep(y, each=4)
dataset$Site=rep(z, each=36)
dataset$Treatment=rep(v, each=12)

dataset$Country<-ifelse(dataset$Site %in% c("Arlon","Zwarte Beek"),"Belgium",
                        ifelse(dataset$Site %in% c("Binnenveld","Drentse Aa"),"Netherlands",
                        ifelse(dataset$Site %in% c("Gutzkow","Peene mouth","Kiel","Recknitz"),"Germany",
                        ifelse(dataset$Site %in% c("Biebrza","Rospuda","Suwalszczyzna","Mazury"),"Poland","UK"))))


dataset<- subset(dataset, select = -c(ID,Plate,New_ID,Sample_NO))
dataset$ID<-paste(dataset$Site,"_",substr(dataset$Treatment,1,1),"_",substr(dataset$Depth,1,1),"_",dataset$Rep,sep="")

dataset<-dataset %>% 
  select(ID,Country,Site,Treatment,Depth,Rep,everything())

str(dataset)
dataset$spc <- dataset[, grep("^[1-9]+",colnames(dataset))]
dataset <- dataset[, !grepl("^[1-9]+",colnames(dataset))]

#identify all character columns
chars <- sapply(dataset$spc, is.character)
#convert all character columns to numeric
dataset$spc[ , chars] <- as.data.frame(apply(dataset$spc[ , chars], 2, as.numeric))


rownames(dataset) <- 1:468

# #Final dataset
# 
# dataset$spc<- -log10(dataset$spc)
 #dataset$spc<- 2-log10(dataset$spc*100)
########################################################################
#Data conversion from reflection to absorbance after KBr corrections
########################################################################


data_KBr <- list.files(path = "6/exp/", # Identify all CSV files
                       pattern = "*.csv", full.names = TRUE) %>%
  lapply(read_csv) %>%                              # Store all files in list
  reduce(full_join, by = "Created as New Dataset")                      # Full-join data sets into one data set

data_KBr


b<-t(data_KBr[-1,])
colnames(b) <-round(as.numeric(b[1, ]))
c<-b[-1,]
Plate_6 <- rowid_to_column(as.data.frame(c), var = "ID")
Plate_6[,1]<-rownames(c)
Plate_6$Plate  <- "P6"
Plate_6$New_ID  <-paste(Plate_6$Plate,"_",as.numeric(gsub("([0-9]+).*$", "\\1", Plate_6$ID)),sep="")
Plate_6$Sample_NO <-as.numeric(gsub("([0-9]+).*$", "\\1", Plate_6$ID))
Plate_6<-Plate_6[order(Plate_6$Sample_NO),]
rownames(Plate_6)=Plate_6$Sample_NO

Plate_6$spc <- Plate_6[, grep("^[1-9]+",colnames(Plate_6))]
Plate_6 <- Plate_6[, !grepl("^[1-9]+",colnames(Plate_6))]

#identify all character columns
chars <- sapply(Plate_6$spc, is.character)
#convert all character columns to numeric
Plate_6$spc[ , chars] <- as.data.frame(apply(Plate_6$spc[ , chars], 2, as.numeric))

Plate_6<-Plate_6[,! names(Plate_6)%in% c("New_ID","Plate","ID","Sample_NO")]

KBr_spc.dat <-colMeans(Plate_6)
KBr_spc.dat <-as.data.frame(t(KBr_spc.dat))

# ################################################################################
# ###FINALL_dataset######################################################
# ###########################################################################
#
dataset$spc<- 2-log10(sweep(dataset$spc, 2, unlist(KBr_spc.dat[1, ]), `/`)*100)
#dataset$spc<- sweep(dataset$spc, 2, unlist(KBr_spc.dat[1, ]), `/`)*100

#save(dataset, file = "MIRs.RData")


#######################################################################################################
