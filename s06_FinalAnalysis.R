#Create final figures
#This is final analysis 

##use new prediction annotation for tumor class from "InhouseMeningioma_CLIPPR..."
#need to merge new predictions with old data set. Not all names in either set were congruent. Steps
#below were used to determine the samples that needed to be manually fixed. 

#load InhouseMeningioma_Clippr_results and bulk sample data
InhouseMeningioma_CLIPPR_results <- read.csv("~/Desktop/PatelLab/FinalFigureAnalysis/InhouseMeningioma_CLIPPR_results.csv", row.names=1)
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
##After above steps, InhouseMeningioma... should have all overlapping samples. Steps below require the modified version of InhouseMeningioma
##reload sample data
load("~/Desktop/PatelLab/09192022_GE_RAW_Normalized_N330PrimarySamples_GeneSymbol.rda")
#continue analysis
newpred<-InhouseMeningioma_CLIPPR_results[,c(1,7)]
samples$rowname<-row.names(samples)
newpred$ID<-row.names(newpred)
newpredsamples<-merge(samples,newpred, by.x = "Sample.Name", by.y = "ID")
row.names(newpredsamples)<-newpredsamples$rowname

###EXTRA CHECKPOINT (NOT necessary for main code)
#Samples with repeat names 16-01-031 (2),18-01-010 (2), 19-01-094 (2), 20-01-027 (3)
#subset newpredsamples to identify if repeats are an issue
subset<-newpredsamples[newpredsamples$Sample.Name=="20-01-027",]
subset$predictions_meta_class
#not an issue. The ID's are the same for different bulk samples because they are likely from the same patient

##pheatmap analysis
##Cibersort Deconvolution
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
newmodulenames<-c("C1_Cibersort", "C2_Cibersort", "C3_Cibersort", "C4_Cibersort", "C5_Cibersort")
row.names(Moduledeconv)<-newmodulenames

#load in bayesprism deconvolution
bayesprism_deconvolution <- read.csv("~/Desktop/PatelLab/Analysis_Results/Deconvolution/bayesprism_deconvolution.csv", row.names=1)
bpmodulenames<-c("C1_BP", "C2_BP", "C3_BP", "C4_BP", "C5_BP")
colnames(bayesprism_deconvolution)<-bpmodulenames

#merge bayes prism and cibersort deconvolution
Moduledeconv<-t(Moduledeconv)
Moduledeconv<-as.data.frame(Moduledeconv)
Moduledeconv$rownames<-row.names(Moduledeconv)

bayesprism_deconvolution$rownames<-row.names(bayesprism_deconvolution)

Moduledeconv<-merge(Moduledeconv, bayesprism_deconvolution, by.x = "rownames", by.y = "rownames")

row.names(Moduledeconv)<-Moduledeconv$rownames
Moduledeconv<-Moduledeconv[,-1]

Moduledeconv<-t(Moduledeconv)
  
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
Moduledeconvindex<-read.csv("~/Desktop/PatelLab/Analysis_Results/Deconvolution/Cibersort_Module_deconvolution.csv")
Moduledeconvindex<-Moduledeconvindex[,1:6]
ID<-row.names(recover)
index<-""
for(x in 1:length(ID)){
  index<-c(index, which(Moduledeconvindex[,1]==ID[x]))
}
index<-index[-1]
index<-as.numeric(index)
Moduledeconv<-Moduledeconv[,index]


##clean up recover data frame
recover<-recover[,c(-1,-3)]
recover<-as.data.frame(recover)
row.names(recover)<-colnames(Moduledeconv)
colnames(recover)<-"Recurrence Free Time"
recover$`Recurrence Free Time`<-as.numeric(recover$`Recurrence Free Time`)
class(recover$`Recurrence Free Time`)


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

?pheatmap(0)
r<-pheatmap(Moduledeconv, annotation_col = recovertumorclass, scale = "row", color = color, breaks = breaks, fontsize = 8, 
            labels_col = shortbulknames, annotation_colors = annoCol, cluster_rows = F)
setwd("~/Desktop/PatelLab/FinalFigureAnalysis/")
ggsave('heatmap_modulebulkdata_deconvolution_withrecover.pdf', plot = r, width=30, height=5)

###Create geom boxplot (ggplot: facet_wrap with C1_Ciber, C1_BP, ... ###################################################
#ALSO x-axis with my Tumor classes, and y axis the ratios?#############################################################
#add tumorclass 
Moduledeconv<-t(Moduledeconv)
neworder<-(sort(row.names(Moduledeconv)))
Moduledeconv<-Moduledeconv[neworder,]

neworder<-(sort(row.names(recovertumorclass)))
recovertumorclass<-recovertumorclass[neworder,]
Moduledeconv<-as.data.frame(Moduledeconv)

Moduledeconv$tumorclass<-as.factor(recovertumorclass$TumorClass)

#create new dataframe with three columns only: ratios, tumorclass, and module/software
for(x in 1:10){
  if(x ==1){
  boxplotdf<-data.frame("Ratios" = Moduledeconv[,x], "TumorClass" = Moduledeconv[,11], 
                        "Module_Program" = rep(colnames(Moduledeconv)[x], nrow(Moduledeconv)))
  }else{ 
    new<-data.frame("Ratios" = Moduledeconv[,x], "TumorClass" = Moduledeconv[,11], 
                    "Module_Program" = rep(colnames(Moduledeconv)[x], nrow(Moduledeconv)))
    boxplotdf<-rbind(boxplotdf, new)
  }
}


#create boxplot
library(ggplot2)
# Basic box plot
b<-ggplot(boxplotdf, aes(x=TumorClass, y=Ratios, fill = TumorClass)) + 
  geom_boxplot()+ facet_wrap(~ Module_Program) + labs(x = "Tumor Class", y = "Deconvolution Ratios")+
  scale_fill_manual(values=c('#28A860','#0042ED','#FE4733'))+
  theme(legend.position = "none")

setwd("~/Desktop/PatelLab/FinalFigureAnalysis/")
ggsave('boxplot_deconvolutionmodulesbyclass.pdf', plot = b, width=10, height=8)
