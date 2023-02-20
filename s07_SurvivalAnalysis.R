#Survival Analysis
library(dplyr)
library(survival)
library(survminer)

#load in data sets: Randomforest classifications and Patient data
SeqTumorClinicalData_12_30_22_anon <- read.csv("~/Desktop/PatelLab/SeqTumorClinicalData_12_30_22 anon.csv")
randomforestclass <- read.csv("~/Desktop/PatelLab/Analysis_Results/randomforestclass.csv", row.names=1)

#both data sets have differently formatted names. Need to merge the two data sets.
#load in bulk sample data which has both formatted names
load("~/Desktop/PatelLab/09192022_GE_RAW_Normalized_N330PrimarySamples_GeneSymbol.rda")

Names<-data.frame("Short" = samples$Sample.Name, "Long" = row.names(samples))
colnames(randomforestclass)<-"RandForest_Prediction"
randomforestclass$rowname<-row.names(randomforestclass)

formatrandforest<-merge(randomforestclass, Names, by.x = "rowname", by.y = "Long")

survivaldata<-merge(formatrandforest, SeqTumorClinicalData_12_30_22_anon, by.x = "Short", by.y = "Sample.Name")
row.names(survivaldata)<-survivaldata$rowname
colnames(survivaldata)[which(colnames(survivaldata)=="RandForest_Prediction")]<-"TumorClass"
survivaldata<-survivaldata[,-2]
#Clean up survival data
table(survivaldata$Vital.Status)
survivaldata<-survivaldata[-which(survivaldata$Vital.Status=="Lost to Follow-up"),]
survivaldata<-survivaldata[-which(survivaldata$Vital.Status==""),]
survivaldata$Vital.Status<-gsub("alive", "Alive", survivaldata$Vital.Status)
survivaldata$Vital.Status<-gsub("ALive", "Alive", survivaldata$Vital.Status)
survivaldata$Vital.Status<-gsub("Alive", 0, survivaldata$Vital.Status)
survivaldata$Vital.Status<-gsub("Dead", 1, survivaldata$Vital.Status)
survivaldata$Vital.Status<-as.numeric(survivaldata$Vital.Status)

survivaldata$Recurrence.Free.Time[survivaldata$Recurrence.Free.Time=="#NUM!"]<-NA
survivaldata$Recurrence.Free.Time[survivaldata$Recurrence.Free.Time==""]<-NA
survivaldata$Recurrence.Free.Time<-as.numeric(survivaldata$Recurrence.Free.Time)
survivaldata$Recurrence.Free.Time<-survivaldata$Recurrence.Free.Time/12

survivaldata<-survivaldata[-which(survivaldata$Recurrence==""),]
survivaldata<-survivaldata[-which(survivaldata$Recurrence=="n/a"),]
survivaldata$Recurrence<-gsub("N", 0, survivaldata$Recurrence)
survivaldata$Recurrence<-gsub("Y", 1, survivaldata$Recurrence)
survivaldata$Recurrence<-as.numeric(survivaldata$Recurrence)

#Survival analysis
sfit<-survfit(Surv(Recurrence.Free.Time, Recurrence)~TumorClass, data = survivaldata)
table(survivaldata$TumorClass, survivaldata$Vital.Status)
ggsurvplot(sfit)+
  labs(title = "Random Forest Prediction Survival Plot", x = "Time (Years)", y = "Recurrence Survival Probability")

summary(sfit)

######################################################################################################################
###############################################################################################################################
####################################################################################################################################
#repeat above steps but with BCM predictions
#load in data sets: BCM classifications and Patient data
SeqTumorClinicalData_12_30_22_anon <- read.csv("~/Desktop/PatelLab/SeqTumorClinicalData_12_30_22 anon.csv")
InhouseMeningioma_CLIPPR_results <- read.csv("~/Desktop/PatelLab/FinalFigureAnalysis/InhouseMeningioma_CLIPPR_results.csv", row.names=1)

#fix InhouseMengingioma_CLIPPR_results
load("~/Desktop/PatelLab/09192022_GE_RAW_Normalized_N330PrimarySamples_GeneSymbol.rda")


#load in newpredsamples
newpred<-InhouseMeningioma_CLIPPR_results[,c(1,7)]
head(newpred)
intersecting<-intersect(samples$Sample.Name, row.names(newpred))
index<-""
for(x in 1:length(intersecting)){
  index<-c(index,which(row.names(newpred)==intersecting[x]))
}
index<-index[-1]
index<-as.numeric(index)

newpred<-newpred[-index,]


index<-""
for(x in 1:length(intersecting)){
  index<-c(index,which(samples$Sample.Name==intersecting[x]))
}
index<-index[-1]
index<-as.numeric(index)

samples<-samples[-index,]

newpred<-newpred[c(-3,-4),]
newpred$samplename<-samples[1:2,1]
excludedsamples<-newpred
##two samples whose names needed manually fixing. The format of the original names was switched. 
row.names(InhouseMeningioma_CLIPPR_results)[1]<-excludedsamples[1,3]
row.names(InhouseMeningioma_CLIPPR_results)[2]<-excludedsamples[2,3]

#merge data sets
InhouseMeningioma_CLIPPR_results$Name<-row.names(InhouseMeningioma_CLIPPR_results)

survivaldata<-merge(InhouseMeningioma_CLIPPR_results, SeqTumorClinicalData_12_30_22_anon, by.x = "Name", by.y = "Sample.Name")

row.names(survivaldata)<-survivaldata$Name
colnames(survivaldata)[which(colnames(survivaldata)=="predictions_meta_class")]<-"TumorClass"
survivaldata<-survivaldata[,-1]
#Clean up survival data
table(survivaldata$Vital.Status)
survivaldata<-survivaldata[-which(survivaldata$Vital.Status=="Lost to Follow-up"),]
survivaldata<-survivaldata[-which(survivaldata$Vital.Status==""),]
survivaldata$Vital.Status<-gsub("alive", "Alive", survivaldata$Vital.Status)
survivaldata$Vital.Status<-gsub("ALive", "Alive", survivaldata$Vital.Status)
survivaldata$Vital.Status<-gsub("Alive", 0, survivaldata$Vital.Status)
survivaldata$Vital.Status<-gsub("Dead", 1, survivaldata$Vital.Status)
survivaldata$Vital.Status<-as.numeric(survivaldata$Vital.Status)


survivaldata$Recurrence.Free.Time[survivaldata$Recurrence.Free.Time=="#NUM!"]<-NA
survivaldata$Recurrence.Free.Time[survivaldata$Recurrence.Free.Time==""]<-NA
survivaldata$Recurrence.Free.Time<-as.numeric(survivaldata$Recurrence.Free.Time)
survivaldata$Recurrence.Free.Time<-survivaldata$Recurrence.Free.Time/12

table(survivaldata$Recurrence)
survivaldata<-survivaldata[-which(survivaldata$Recurrence==""),]
survivaldata<-survivaldata[-which(survivaldata$Recurrence=="n/a"),]
survivaldata$Recurrence<-gsub("N", 0, survivaldata$Recurrence)
survivaldata$Recurrence<-gsub("Y", 1, survivaldata$Recurrence)
survivaldata$Recurrence<-as.numeric(survivaldata$Recurrence)


#Survival analysis
sfit<-survfit(Surv(Recurrence.Free.Time, Recurrence)~TumorClass, data = survivaldata)
table(survivaldata$TumorClass, survivaldata$Vital.Status)

ggsurvplot(sfit)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "BCM Prediction Tumor Class Survival Plot", x = "Time (Years)", y = "Recurrence Survival Probability")
setwd("~/Desktop/PatelLab/FinalFigureAnalysis/")



##########Survival Data split up by modules########################################################################
################################################################################################################################################
#load in patient data
SeqTumorClinicalData_12_30_22_anon <- read.csv("~/Desktop/PatelLab/SeqTumorClinicalData_12_30_22 anon.csv")

#load in module data: TWO OPTIONS
#Cibersort module data OR BayesPrism Module data
#OPTION 1: BayesPrism Module Data
Moduledeconv<-read.csv("~/Desktop/PatelLab/Analysis_Results/Deconvolution/bayesprism_deconvolution.csv")
colnames(Moduledeconv)[1]<-"Mixture"
Moduledeconv<-Moduledeconv[,1:6]
#Option 2: Cibersort Module data
Moduledeconv<-read.csv("~/Desktop/PatelLab/Analysis_Results/Deconvolution/Cibersort_Module_Deconvolution.csv")
Moduledeconv<-Moduledeconv[,1:6]


#load in newpredsample object with BCM projections
InhouseMeningioma_CLIPPR_results <- read.csv("~/Desktop/PatelLab/FinalFigureAnalysis/InhouseMeningioma_CLIPPR_results.csv", row.names=1)


#need to switch format of sample names in patient data to match sample names in module ratios
#load in bulk sample data which has both formatted names
load("~/Desktop/PatelLab/09192022_GE_RAW_Normalized_N330PrimarySamples_GeneSymbol.rda")

#fix InhouseMengingioma_CLIPPR_results
#load in newpredsamples
newpred<-InhouseMeningioma_CLIPPR_results[,c(1,7)]
head(newpred)
intersecting<-intersect(samples$Sample.Name, row.names(newpred))
index<-""
for(x in 1:length(intersecting)){
  index<-c(index,which(row.names(newpred)==intersecting[x]))
}
index<-index[-1]
index<-as.numeric(index)

newpred<-newpred[-index,]


index<-""
for(x in 1:length(intersecting)){
  index<-c(index,which(samples$Sample.Name==intersecting[x]))
}
index<-index[-1]
index<-as.numeric(index)

samples<-samples[-index,]

newpred<-newpred[c(-3,-4),]
newpred$samplename<-samples[1:2,1]
excludedsamples<-newpred
##two samples whose names needed manually fixing. The format of the original names was switched. 
row.names(InhouseMeningioma_CLIPPR_results)[1]<-excludedsamples[1,3]
row.names(InhouseMeningioma_CLIPPR_results)[2]<-excludedsamples[2,3]

#need to switch format of sample names in patient data to match sample names in module ratios
#combine inhousemenigioma, Moduledeconv, and SeqTumorclinicaldata

#reload in bulk sample data
load("~/Desktop/PatelLab/09192022_GE_RAW_Normalized_N330PrimarySamples_GeneSymbol.rda")

#add both format names to the seqtumorclinicaldata obj
Names<-data.frame("Short" = samples$Sample.Name, "Long" = row.names(samples))

formatpatient<-merge(Names, SeqTumorClinicalData_12_30_22_anon, by.x = "Short", by.y = "Sample.Name")
which(colnames(formatpatient)=="Recurrence")
formatpatient<-formatpatient[,c(1,2,25,27,31)]

#add both format names to the Inhousemeningioma obj
InhouseMeningioma_CLIPPR_results$Sample.Name<-row.names(InhouseMeningioma_CLIPPR_results)
formatinhouse<-merge(Names, InhouseMeningioma_CLIPPR_results, by.x = "Short", by.y = "Sample.Name" )
formatinhouse<-formatinhouse[, c(1,2,which(colnames(formatinhouse)=="predictions_meta_class"))]
colnames(formatinhouse)<-c("Short", "Long", "TumorClass")

#transform module ratios into z-scale
Moduledeconv[2:6] <- t(scale(t(Moduledeconv[2:6])))

#rename modules C1: turq, C2: blue, C3: green, C4: yellow, C5: brown
newmodulenames<-c("Mixture","C1", "C2", "C3", "C4", "C5")
colnames(Moduledeconv)<-newmodulenames

#combine module ratios with patient vital data and tumorclass data
part1<-merge(formatpatient, Moduledeconv, by.x = "Long", by.y = "Mixture")
survivaldata<-merge(part1, formatinhouse, by.x = "Long", by.y = "Long")

#assign long label names to row.names
row.names(survivaldata)<-survivaldata$Long
survivaldata<-survivaldata[,c(-1,-2, -11)]

#Clean up survival data
table(survivaldata$Vital.Status)
survivaldata<-survivaldata[-which(survivaldata$Vital.Status=="Lost to Follow-up"),]
survivaldata<-survivaldata[-which(survivaldata$Vital.Status==""),]
survivaldata$Vital.Status<-gsub("alive", "Alive", survivaldata$Vital.Status)
survivaldata$Vital.Status<-gsub("ALive", "Alive", survivaldata$Vital.Status)
survivaldata$Vital.Status<-gsub("Alive", 0, survivaldata$Vital.Status)
survivaldata$Vital.Status<-gsub("Dead", 1, survivaldata$Vital.Status)
survivaldata$Vital.Status<-as.numeric(survivaldata$Vital.Status)


survivaldata$Recurrence.Free.Time[survivaldata$Recurrence.Free.Time=="#NUM!"]<-NA
survivaldata$Recurrence.Free.Time[survivaldata$Recurrence.Free.Time==""]<-NA
survivaldata$Recurrence.Free.Time<-as.numeric(survivaldata$Recurrence.Free.Time)
survivaldata$Recurrence.Free.Time<-survivaldata$Recurrence.Free.Time/12

survivaldata<-survivaldata[-which(survivaldata$Recurrence==""),]
survivaldata<-survivaldata[-which(survivaldata$Recurrence=="n/a"),]
survivaldata$Recurrence<-gsub("N", 0, survivaldata$Recurrence)
survivaldata$Recurrence<-gsub("Y", 1, survivaldata$Recurrence)
survivaldata$Recurrence<-as.numeric(survivaldata$Recurrence)


# set Z-scale cut-offs for high and low expression
survivaldata$C1 <- ifelse(survivaldata$C1 > 0, 'High', 'Low')
survivaldata$C2 <- ifelse(survivaldata$C2 > 0, 'High', 'Low')
survivaldata$C3 <- ifelse(survivaldata$C3 > 0, 'High', 'Low')
survivaldata$C4 <- ifelse(survivaldata$C4 > 0, 'High', 'Low')
survivaldata$C5 <- ifelse(survivaldata$C5 > 0, 'High', 'Low')
#C1 survival plot

sfit<-survfit(Surv(Recurrence.Free.Time, Recurrence)~C1+TumorClass, data = survivaldata)

ggsurvplot(sfit, pval = T)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "Cibersort Module C1 Survival Plot", x = "Time (Years)", y = "Recurrence Survival Probability")

#C2 survival plot

sfit<-survfit(Surv(Recurrence.Free.Time, Recurrence)~C2+TumorClass, data = survivaldata)

ggsurvplot(sfit, pval = T)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "Cibersort Module C2 Survival Plot", x = "Time (Years)", y = "Recurrence Survival Probability")
#C3 survival plot

sfit<-survfit(Surv(Recurrence.Free.Time, Recurrence)~C3+TumorClass, data = survivaldata)

ggsurvplot(sfit, pval = T)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "Cibersort Module C3 Survival Plot", x = "Time (Years)", y = "Recurrence Survival Probability")
#C4 survival plot

sfit<-survfit(Surv(Recurrence.Free.Time, Recurrence)~C4+TumorClass, data = survivaldata)

ggsurvplot(sfit, pval = T)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "Cibersort Module C4 Survival Plot", x = "Time (Years)", y = "Recurrence Survival Probability")
#C5 survival plot

sfit<-survfit(Surv(Recurrence.Free.Time, Recurrence)~C5+TumorClass, data = survivaldata)

ggsurvplot(sfit, pval = T)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "Cibersort Module C5 Survival Plot", x = "Time (Years)", y = "Recurrence Survival Probability")



