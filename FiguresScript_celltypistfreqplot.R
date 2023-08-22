#Create Cell Typist Plots, Module Score Boxplots, + PCA Tumor Type plot

#Load in immune cell clusters
setwd("~/Desktop/PatelLab/NewAnalysis")
load("072523_allLocations_immunecellsonly_modulescores.rda")
immune<-seuratObj

###Frequency bar plots
#Class A
classA<-subset(immune, subset = TumorType == "A")

types<-levels(factor(classA$celltypist))
freq<-1:length(levels(factor(classA$celltypist)))
for (ident in 1:length(levels(factor(classA$celltypist)))) {
  freq[ident]<-sum(classA$celltypist==types[ident])/length(classA$celltypist)
}


ClassAMajVot<-data.frame(types,freq,"TumorType" = rep("A",length(types)))
#class B
classB<-subset(immune, subset = TumorType == "B")

types<-levels(factor(classB$celltypist))
freq<-1:length(levels(factor(classB$celltypist)))
for (ident in 1:length(levels(factor(classB$celltypist)))) {
  freq[ident]<-sum(classB$celltypist==types[ident])/length(classB$celltypist)
}


ClassBMajVot<-data.frame(types,freq,"TumorType" = rep("B",length(types)))
#Class C
classC<-subset(immune, subset = TumorType == "C")

types<-levels(factor(classC$celltypist))
freq<-1:length(levels(factor(classC$celltypist)))
for (ident in 1:length(levels(factor(classC$celltypist)))) {
  freq[ident]<-sum(classC$celltypist==types[ident])/length(classC$celltypist)
}


ClassCMajVot<-data.frame(types,freq, "TumorType" = rep("C",length(types)))


#bind data frames
MajVotingFreq<-rbind(ClassAMajVot,ClassBMajVot, ClassCMajVot)
#create bar plot
library(tidyverse)
ggplot(data=MajVotingFreq, aes(x=freq, y=types, fill = TumorType))+
geom_bar(stat="identity", width=0.5, position=position_dodge())+
theme(axis.text.x = element_text(size = 5))+
  xlab("Relative Frequency")+
  ylab("Cell Typist Annotations")+
  scale_fill_manual(values=c('#28A860','#0042ED', "#FE4733"))

######create boxplot with module scores and Cell typist annotations'
#load in gene set
PanGenes <- read.csv("~/Desktop/PatelLab/PanImmuneGeneSet/PanGenes.csv")
sets<-levels(factor(PanGenes$SetName))
sets<-sets[-which(sets=="CD103pos_CD103neg_ratio_25446897")]
sets<-sets[-which(sets=="GP11_Immune_IFN")]


#seuratObj AddModuleScore
metadata<-immune@meta.data
ModuleScoreIndex<-which(colnames(metadata)=="GeneSet1"):114

##Split gene sets into three groups in order for cleaner visualization
subset1<-1:36
subset2<-37:72
subset3<-73:107
sets<-sets[subset3]
ModuleScoreIndex<-ModuleScoreIndex[subset3]



##create boxplot with module scores and Cell typist annotations'
#Class A
classA<-subset(immune, subset = TumorType == "A")
metadata<-classA@meta.data
for (geneset in 1:length(sets)) {
  if(geneset ==1){
    score<-metadata[,ModuleScoreIndex[geneset]]
  } else {
    score<-c(score,metadata[,ModuleScoreIndex[geneset]])
  }
}

for(geneset in 1:length(sets)){
  if(geneset == 1){
    setname<-rep(sets[geneset],nrow(metadata))
  } else {
    setname<-c(setname,rep(sets[geneset],nrow(metadata)))
  }
}

ClassAScore<-data.frame(setname,score,"TumorType" = rep("A",length(score)))

#Class B
classB<-subset(immune, subset = TumorType == "B")
metadata<-classB@meta.data
for (geneset in 1:length(sets)) {
  if(geneset ==1){
    score<-metadata[,ModuleScoreIndex[geneset]]
  } else {
    score<-c(score,metadata[,ModuleScoreIndex[geneset]])
  }
}

for(geneset in 1:length(sets)){
  if(geneset == 1){
    setname<-rep(sets[geneset],nrow(metadata))
  } else {
    setname<-c(setname,rep(sets[geneset],nrow(metadata)))
  }
}

ClassBScore<-data.frame(setname,score,"TumorType" = rep("B",length(score)))

#Class C
classC<-subset(immune, subset = TumorType == "C")
metadata<-classC@meta.data
for (geneset in 1:length(sets)) {
  if(geneset ==1){
    score<-metadata[,ModuleScoreIndex[geneset]]
  } else {
    score<-c(score,metadata[,ModuleScoreIndex[geneset]])
  }
}

for(geneset in 1:length(sets)){
  if(geneset == 1){
    setname<-rep(sets[geneset],nrow(metadata))
  } else {
    setname<-c(setname,rep(sets[geneset],nrow(metadata)))
  }
}

ClassCScore<-data.frame(setname,score,"TumorType" = rep("C",length(score)))

##bind rows together
combineddf<-rbind(ClassAScore,ClassBScore,ClassCScore)
#make ggplot
library(tidyverse)
ggplot(data=combineddf, aes(x=score, y=setname, fill = TumorType))+
  geom_boxplot(outlier.size=0.5)+
  theme(axis.text.x = element_text(size = 8))+
  xlab("Module Score")+
  ylab("Pan Gene Set")+
  theme(axis.text = element_text(size = 8))+xlim(-1.5,4.5)

##################PCA analysis
#load in gene set
PanGenes <- read.csv("~/Desktop/PatelLab/PanImmuneGeneSet/PanGenes.csv")
sets<-levels(factor(PanGenes$SetName))
sets<-sets[-which(sets=="CD103pos_CD103neg_ratio_25446897")]
sets<-sets[-which(sets=="GP11_Immune_IFN")]


#seuratObj AddModuleScore
metadata<-immune@meta.data
ModuleScoreIndex<-which(colnames(metadata)=="GeneSet1"):114
mdata<-immune@meta.data
immunepca<-prcomp(mdata[,ModuleScoreIndex], center = TRUE, scale. = TRUE)
library(devtools)
library(ggfortify)
autoplot(immunepca, data = mdata, colour = 'TumorType', )
summary(immunepca)


