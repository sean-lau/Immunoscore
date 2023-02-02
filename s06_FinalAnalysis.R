#Create final figures
#This is final analysis 

##use new prediction annotation for tumor clouse from "InhouseMeningioma_CLIPPR..."
#need to merge new predictions with old data set. Not all names in either set were congruent. Steps
#below were used to determine the samples that needed to be manually fixed. 

#load InhouseMeningioma_Clippr_results and bulk sample data
InhouseMeningioma_CLIPPR_results <- read.csv("InhouseMeningioma_CLIPPR_results.csv")
load("~/Desktop/PatelLab/09192022_GE_RAW_Normalized_N330PrimarySamples_GeneSymbol.rda")


newpred<-InhouseMeningioma_CLIPPR_results[,c(1,8)]
head(newpred)
intersecting<-intersect(samples$Sample.Name, newpred$ID)
index<-""
for(x in 1:length(intersecting)){
  index<-c(index,which(newpred$ID==intersecting[x]))
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
##two samples whose names needed manually fixing. Two samples did not overlap with the samples from the bulk data
InhouseMeningioma_CLIPPR_results[1,1]<-excludedsamples[1,3]
InhouseMeningioma_CLIPPR_results[2,1]<-excludedsamples[2,3]

##After above steps, InhouseMeningioma... should have all overlapping genes
##reload sample data
load("~/Desktop/PatelLab/09192022_GE_RAW_Normalized_N330PrimarySamples_GeneSymbol.rda")
#continue analysis
newpred<-InhouseMeningioma_CLIPPR_results[,c(1,2,8)]
samples$rowname<-row.names(samples)
newpredsamples<-merge(samples,newpred, by.x = "Sample.Name", by.y = "ID")
row.names(newpredsamples)<-newpredsamples$rowname
table(newpredsamples$Sample.Name)

###EXTRA CHECKPOINT (NOT necessary for main code)
#Samples with repeat names 16-01-031 (2),18-01-010 (2), 19-01-094 (2), 20-01-027 (3)
#subset newpredsamples to identify if repeats are an issue
subset<-newpredsamples[newpredsamples$Sample.Name=="20-01-027",]
subset$predictions_meta_class
#not an issue. The ID's are the same for different bulk samples because they are likely from the same patient

##pheatmap analysis
Moduledeconv<-read.csv("~/Desktop/PatelLab/Analysis_Results/Deconvolution/Cibersort_Module_deconvolution.csv")
Moduledeconv<-Moduledeconv[,1:6]
ID<-row.names(newpredsamples)
index<-""
for(x in 1:length(ID)){
  index<-c(index, which(Moduledeconv[,1]==ID[x]))
}
index<-index[-1]
index<-as.numeric(index)
Moduledeconv<-Moduledeconv[index,]

row.names(Moduledeconv)<-Moduledeconv[,1]
Moduledeconv<-Moduledeconv[,-1]
Moduledeconv<-t(Moduledeconv)

#rename modules C1: turq, C2: blue, C3: green, C4: yellow, C5: brown
newmodulenames<-c("C1", "C2", "C3", "C4", "C5")
row.names(Moduledeconv)<-newmodulenames

#new prediction annotations
sampleclass<-as.data.frame(newpredsamples$predictions_meta_class)
row.names(sampleclass)<-row.names(newpredsamples)
colnames(sampleclass)<-"TumorClass"
shortbulknames<-newpredsamples[,1]
names(shortbulknames)<-row.names(newpredsamples)

#recover free time annotations
#load in recover time data
SeqTumorClinicalData_12_30_22_anon <- read.csv("~/Desktop/PatelLab/SeqTumorClinicalData_12_30_22 anon.csv")
recover<-SeqTumorClinicalData_12_30_22_anon[,c(1,26)]
recover$Recurrence.Free.Time[recover$Recurrence.Free.Time=="#NUM!"]<-NA
recover$Recurrence.Free.Time[recover$Recurrence.Free.Time==0]<-NA
Names<-as.data.frame(newpredsamples[,1])
Names$mainlab<-row.names(newpredsamples)
##note: there are two formats of name labels in data sets so must synchronize data based on one format
names(Names)<-c("datelab", "mainlab"
)
recover<-merge(recover,Names, by.x = "Sample.Name", by.y = "datelab")
row.names(recover)<-recover$mainlab
#subset Moduledeconv to only those in the recover freetime annotations
Moduledeconv<-read.csv("~/Desktop/PatelLab/Analysis_Results/Deconvolution/Cibersort_Module_deconvolution.csv")
Moduledeconv<-Moduledeconv[,1:6]
ID<-row.names(recover)
index<-""
for(x in 1:length(ID)){
  index<-c(index, which(Moduledeconv[,1]==ID[x]))
}
index<-index[-1]
index<-as.numeric(index)
Moduledeconv<-Moduledeconv[index,]

row.names(Moduledeconv)<-Moduledeconv[,1]
Moduledeconv<-Moduledeconv[,-1]
Moduledeconv<-t(Moduledeconv)

#rename modules C1: turq, C2: blue, C3: green, C4: yellow, C5: brown
newmodulenames<-c("C1", "C2", "C3", "C4", "C5")
row.names(Moduledeconv)<-newmodulenames

##clean up recover data frame
recover<-recover[,c(-1,-3)]
recover<-as.data.frame(recover)
row.names(recover)<-colnames(Moduledeconv)
colnames(recover)<-"Recover Time"
recover$`Recover Time`<-as.numeric(recover$`Recover Time`)
class(recover$`Recover Time`)


########add tumor class and randomforest predictions to recover
recover$names<-row.names(recover)
sampleclass$ID<-row.names(sampleclass)
recovertumorclass<-merge(recover, sampleclass, by.x = "names", by.y = "ID")

#load in data frame with random forest class predictions
randomforestclass <- read.csv("~/Desktop/PatelLab/Analysis_Results/randomforestclass.csv", row.names=1)

colnames(randomforestclass)<-"Randomforest_ClassPred"
randomforestclass$names<-row.names(randomforestclass)
recovertumorclass<-merge(recovertumorclass, randomforestclass, by.x = "names", by.y = "names")
row.names(recovertumorclass)<-recovertumorclass$names
recovertumorclass<-recovertumorclass[,-1]

#Manually change the bulk sample names to short names
recover<-SeqTumorClinicalData_12_30_22_anon[,c(1,26)]
recover$Recurrence.Free.Time[recover$Recurrence.Free.Time=="#NUM!"]<-NA
recover$Recurrence.Free.Time[recover$Recurrence.Free.Time==0]<-NA
Names<-as.data.frame(newpredsamples[,1])
Names$mainlab<-row.names(newpredsamples)
recover$mainlab<-row.names(recover)
shortbulknames<-recover$Sample.Name
names(shortbulknames)<-recover$mainlab

#create pheatmap
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
breaks <- seq(-2, 2, length = 16)
color <- colorRampPalette(rev(brewer.pal(11,'RdBu')))(length(breaks))
annoCol<-list(TumorClass=c(A='#28A860', B ='#0042ED', C='#FE4733'),Randomforest_ClassPred=c(A='#28A860', B ='#0042ED', C='#FE4733') )

r<-pheatmap(Moduledeconv, annotation_col = recovertumorclass, scale = "row", color = color, breaks = breaks, fontsize = 8, 
            labels_col = shortbulknames, annotation_colors = annoCol)
setwd("~/Desktop/PatelLab/FinalFigureAnalysis/")
ggsave('heatmap_modulebulkdata_deconvolution_withrecover.pdf', plot = r, width=30, height=5)
