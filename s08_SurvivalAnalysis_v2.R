#Survival Analysis
setwd("~/Desktop/PatelLab/NewAnalysis/")
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)

#load in data sets: Randomforest classifications and Patient data
load("~/Desktop/PatelLab/031023_Inhouse_Raleigh_Hongkong_GeneSymbol.rda")
HK_UCSF_BCM <- read.csv("~/Desktop/PatelLab/HK_UCSF_BCM.csv")


#determine which names aren't overlapping; may be a formatting issue
inter<-intersect(HK_UCSF_BCM$Original_name, row.names(samplesheet))

samplesheet<-samplesheet[row.names(samplesheet)%in%inter,]
samplesheet$name<-row.names(samplesheet)
HK_UCSF_BCM<-HK_UCSF_BCM[HK_UCSF_BCM$Original_name%in%inter,]


survivaldata<-merge(HK_UCSF_BCM,samplesheet, by.x = "Original_name", by.y = "name")

row.names(survivaldata)<-survivaldata$Original_name
survivaldata<-survivaldata[,c(1,23,24,40)]

colnames(survivaldata)[4]<-"TumorClass"
survivaldata<-survivaldata[,-1]
#Clean up survival data
table(survivaldata$Recurrence)
table(survivaldata$Recurrence.Free.Survival)

survivaldata$Recurrence.Free.Survival<-as.numeric(survivaldata$Recurrence.Free.Survival)
survivaldata$Recurrence.Free.Survival<-survivaldata$Recurrence.Free.Survival/12

#Survival analysis
sfit<-survfit(Surv(Recurrence.Free.Survival, Recurrence)~TumorClass, data = survivaldata)
table(survivaldata$TumorClass, survivaldata$Recurrence)

ggsurvplot(sfit)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "BCM Prediction Tumor Class Survival Plot", x = "Time (Years)", y = "Recurrence Survival Probability")



##########Survival Data split up by modules########################################################################
################################################################################################################################################
#load in patient data
HK_UCSF_BCM <- read.csv("~/Desktop/PatelLab/HK_UCSF_BCM.csv")

#load in sample sheet
load("~/Desktop/PatelLab/031023_Inhouse_Raleigh_Hongkong_GeneSymbol.rda")

#load in module data: TWO OPTIONS
#BayesPrism module data OR BayesPrism Module data
#OPTION 1: BayesPrism Module Data
Moduledeconv<-bayesprism_combined_deconvolution <- read.csv("~/Desktop/PatelLab/NewAnalysis/Deconvolution/bayesprism_combined_deconvolution.csv", row.names=1)
colnames(Moduledeconv)<-gsub("^C","cluster", colnames(Moduledeconv))

#Option 2: Cibersort Module data
Moduledeconv<- CIBERSORTx_Results <- read.csv("~/Desktop/PatelLab/NewAnalysis/Deconvolution/CIBERSORTx_BCMsc_Results.csv", row.names=1)
Moduledeconv<-Moduledeconv[,1:8]



#transform module ratios into z-scale
Moduledeconv <- scale(Moduledeconv)


#combine module ratios with patient vital data and tumorclass data
Moduledeconv<-as.data.frame(Moduledeconv)
Moduledeconv$name<-row.names(Moduledeconv)

samplesheet$name<-row.names(samplesheet)
samplesheet<-samplesheet[,-1]

HK_UCSF_BCM<-HK_UCSF_BCM[,c(23:24,4) ]

part1<-merge(Moduledeconv, samplesheet, by.x = "name", by.y = "name")
survivaldata<-merge(part1, HK_UCSF_BCM, by.x = "name", by.y = "Original_name")

#only 613 overlapping bulk samples. Possibly lack of total overlap due to name formatting

#assign long label names to row.names
row.names(survivaldata)<-survivaldata$name
survivaldata<-survivaldata[,c(-1
                            )]

#Clean up survival data
survivaldata$Recurrence<-as.numeric(survivaldata$Recurrence)
survivaldata$Recurrence.Free.Survival<-survivaldata$Recurrence.Free.Survival/12


# set Z-scale cut-offs for high and low expression
survivaldata$cluster1 <- ifelse(survivaldata$cluster1 > 0, 'High', 'Low')
survivaldata$cluster2 <- ifelse(survivaldata$cluster2 > 0, 'High', 'Low')
survivaldata$cluster3 <- ifelse(survivaldata$cluster3 > 0, 'High', 'Low')
survivaldata$cluster6 <- ifelse(survivaldata$cluster6 > 0, 'High', 'Low')
survivaldata$cluster5 <- ifelse(survivaldata$cluster5 > 0, 'High', 'Low')
survivaldata$cluster7 <- ifelse(survivaldata$cluster7 > 0, 'High', 'Low')
survivaldata$cluster8 <- ifelse(survivaldata$cluster8 > 0, 'High', 'Low')
survivaldata$cluster9 <- ifelse(survivaldata$cluster9 > 0, 'High', 'Low')
survA<-survivaldata[survivaldata$class == "A",]
survB<-survivaldata[survivaldata$class == "B",]
survC<-survivaldata[survivaldata$class == "C",]

####C1 survival plots

#Class A
sfit<-survfit(Surv(Recurrence.Free.Survival, Recurrence)~cluster1, data = survA)

ggsurvplot(sfit, pval = T, risk.table = T)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "BP Module C1 Survival Plot Class A", x = "Time (Years)", y = "Recurrence Survival")
#Class B
sfit<-survfit(Surv(Recurrence.Free.Survival, Recurrence)~cluster1, data = survB)

ggsurvplot(sfit, pval = T, risk.table = T)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "BP Module C1 Survival Plot Class B", x = "Time (Years)", y = "Recurrence Survival")

#Class C
sfit<-survfit(Surv(Recurrence.Free.Survival, Recurrence)~cluster1, data = survC)

ggsurvplot(sfit, pval = T, risk.table = T)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "BP Module C1 Survival Plot Class C", x = "Time (Years)", y = "Recurrence Survival")

####C2 survival plots
#Class A
sfit<-survfit(Surv(Recurrence.Free.Survival, Recurrence)~cluster2, data = survA)

ggsurvplot(sfit, pval = T, risk.table = T)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "BP Module C2 Survival Plot Class A", x = "Time (Years)", y = "Recurrence Survival")

#Class B
sfit<-survfit(Surv(Recurrence.Free.Survival, Recurrence)~cluster2, data = survB)

ggsurvplot(sfit, pval = T, risk.table = T)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "BP Module C2 Survival Plot Class B", x = "Time (Years)", y = "Recurrence Survival")

#Class C
sfit<-survfit(Surv(Recurrence.Free.Survival, Recurrence)~cluster2, data = survC)

ggsurvplot(sfit, pval = T, risk.table = T)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "BP Module C2 Survival Plot Class C", x = "Time (Years)", y = "Recurrence Survival")

###########C3 survival plot

#Class A
sfit<-survfit(Surv(Recurrence.Free.Survival, Recurrence)~cluster3, data = survA)

ggsurvplot(sfit, pval = T, risk.table = T)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "BP Module C3 Survival Plot Class A", x = "Time (Years)", y = "Recurrence Survival")

#Class B
sfit<-survfit(Surv(Recurrence.Free.Survival, Recurrence)~cluster3, data = survB)

ggsurvplot(sfit, pval = T, risk.table = T)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "BP Module C3 Survival Plot Class B", x = "Time (Years)", y = "Recurrence Survival")

#Class C
sfit<-survfit(Surv(Recurrence.Free.Survival, Recurrence)~cluster3, data = survC)

ggsurvplot(sfit, pval = T, risk.table = T)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "BP Module C3 Survival Plot Class C", x = "Time (Years)", y = "Recurrence Survival")
#######C5 survival plot

#Class A
sfit<-survfit(Surv(Recurrence.Free.Survival, Recurrence)~cluster5, data = survA)

ggsurvplot(sfit, pval = T, risk.table = T)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "BP Module C5 Survival Plot Class A", x = "Time (Years)", y = "Recurrence Survival")

#Class B
sfit<-survfit(Surv(Recurrence.Free.Survival, Recurrence)~cluster5, data = survB)

ggsurvplot(sfit, pval = T, risk.table = T)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "BP Module C5 Survival Plot Class B", x = "Time (Years)", y = "Recurrence Survival")

#Class C
sfit<-survfit(Surv(Recurrence.Free.Survival, Recurrence)~cluster5, data = survC)

ggsurvplot(sfit, pval = T, risk.table = T)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "BP Module C5 Survival Plot Class C", x = "Time (Years)", y = "Recurrence Survival")
############C6 survival plot

#Class A
sfit<-survfit(Surv(Recurrence.Free.Survival, Recurrence)~cluster6, data = survA)

ggsurvplot(sfit, pval = T, risk.table = T)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "BP Module C6 Survival Plot Class A", x = "Time (Years)", y = "Recurrence Survival")

#Class B
sfit<-survfit(Surv(Recurrence.Free.Survival, Recurrence)~cluster6, data = survB)

ggsurvplot(sfit, pval = T, risk.table = T)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "BP Module C6 Survival Plot Class B", x = "Time (Years)", y = "Recurrence Survival")

#Class C
sfit<-survfit(Surv(Recurrence.Free.Survival, Recurrence)~cluster6, data = survC)

ggsurvplot(sfit, pval = T, risk.table = T)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "BP Module C6 Survival Plot Class C", x = "Time (Years)", y = "Recurrence Survival")

############C7 survival plot

#Class A
sfit<-survfit(Surv(Recurrence.Free.Survival, Recurrence)~cluster7, data = survA)

ggsurvplot(sfit, pval = T, risk.table = T)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "BP Module C7 Survival Plot Class A", x = "Time (Years)", y = "Recurrence Survival")

#Class B
sfit<-survfit(Surv(Recurrence.Free.Survival, Recurrence)~cluster7, data = survB)

ggsurvplot(sfit, pval = T, risk.table = T)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "BP Module C7 Survival Plot Class B", x = "Time (Years)", y = "Recurrence Survival")

#Class C
sfit<-survfit(Surv(Recurrence.Free.Survival, Recurrence)~cluster7, data = survC)

ggsurvplot(sfit, pval = T, risk.table = T)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "BP Module C7 Survival Plot Class C", x = "Time (Years)", y = "Recurrence Survival")


############C8 survival plot

#Class A
sfit<-survfit(Surv(Recurrence.Free.Survival, Recurrence)~cluster8, data = survA)

ggsurvplot(sfit, pval = T, risk.table = T)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "BP Module C8 Survival Plot Class A", x = "Time (Years)", y = "Recurrence Survival")

#Class B
sfit<-survfit(Surv(Recurrence.Free.Survival, Recurrence)~cluster8, data = survB)

ggsurvplot(sfit, pval = T, risk.table = T)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "BP Module C8 Survival Plot Class B", x = "Time (Years)", y = "Recurrence Survival")

#Class C
sfit<-survfit(Surv(Recurrence.Free.Survival, Recurrence)~cluster8, data = survC)

ggsurvplot(sfit, pval = T, risk.table = T)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "BP Module C8 Survival Plot Class C", x = "Time (Years)", y = "Recurrence Survival")

############C9 survival plot

#Class A
sfit<-survfit(Surv(Recurrence.Free.Survival, Recurrence)~cluster9, data = survA)

ggsurvplot(sfit, pval = T, risk.table = T)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "BP Module C9 Survival Plot Class A", x = "Time (Years)", y = "Recurrence Survival")

#Class B
sfit<-survfit(Surv(Recurrence.Free.Survival, Recurrence)~cluster9, data = survB)

ggsurvplot(sfit, pval = T, risk.table = T)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "BP Module C9 Survival Plot Class B", x = "Time (Years)", y = "Recurrence Survival")

#Class C
sfit<-survfit(Surv(Recurrence.Free.Survival, Recurrence)~cluster9, data = survC)

ggsurvplot(sfit, pval = T, risk.table = T)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "BP Module C9 Survival Plot Class C", x = "Time (Years)", y = "Recurrence Survival")

