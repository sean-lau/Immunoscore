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


#Survival analysis
sfit<-survfit(Surv(Recurrence.Free.Time, Vital.Status)~TumorClass, data = survivaldata)
table(survivaldata$TumorClass, survivaldata$Vital.Status)
ggsurvplot(sfit)+
  labs(title = "Random Forest Prediction Survival Plot", x = "Time (Months)")

summary(sfit)

########
#################
######################
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

#Survival analysis
sfit<-survfit(Surv(Recurrence.Free.Time, Vital.Status)~TumorClass, data = survivaldata)
table(survivaldata$TumorClass, survivaldata$Vital.Status)

r <- ggsurvplot(sfit)+ ###For ggsurvplot to work, make sure the vital status variable is not factored i.e. don't use as.factor() or unfactor the variable
  labs(title = "BCM Prediction Tumor Class Survival Plot", x = "Time (Months)")
setwd("~/Desktop/PatelLab/FinalFigureAnalysis/")
ggsave('BCMprediction_tumorclass_survivalplot.pdf', plot = r, width=15, height=10)

?ggsurvplot
