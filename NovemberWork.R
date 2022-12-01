##Novemeber Workspace 
#remove tumor clusters 1, 3, 10, 15, 24, 8, 12, 16, 19

clusters<-levels(sub$seurat_clusters)
clusters<-clusters[c(-2,-4,-11,-16,-25, -9, -13, -17, -20)]
immune<- subset(sub, cells=colnames(sub)[sub$seurat_clusters %in% clusters])
immune$TumorType
###Frequency bar plots
#Class A
classA<-subset(immune, subset = TumorType == "ClassA")

types<-levels(factor(classA$CellTypistAnnotation))
freq<-1:length(levels(factor(classA$CellTypistAnnotation)))
for (ident in 1:length(levels(factor(classA$CellTypistAnnotation)))) {
  freq[ident]<-sum(classA$CellTypistAnnotation==types[ident])/length(classA$CellTypistAnnotation)
}


ClassAMajVot<-data.frame(types,freq,"TumorType" = rep("A",length(types)))
#class B
classB<-subset(immune, subset = TumorType == "ClassB")

types<-levels(factor(classB$CellTypistAnnotation))
freq<-1:length(levels(factor(classB$CellTypistAnnotation)))
for (ident in 1:length(levels(factor(classB$CellTypistAnnotation)))) {
  freq[ident]<-sum(classB$CellTypistAnnotation==types[ident])/length(classB$CellTypistAnnotation)
}


ClassBMajVot<-data.frame(types,freq,"TumorType" = rep("B",length(types)))
#Class C
classC<-subset(immune, subset = TumorType == "ClassC")

types<-levels(factor(classC$CellTypistAnnotation))
freq<-1:length(levels(factor(classC$CellTypistAnnotation)))
for (ident in 1:length(levels(factor(classC$CellTypistAnnotation)))) {
  freq[ident]<-sum(classC$CellTypistAnnotation==types[ident])/length(classC$CellTypistAnnotation)
}


ClassCMajVot<-data.frame(types,freq, "TumorType" = rep("C",length(types)))


#bind data frames
MajVotingFreq<-rbind(ClassAMajVot,ClassBMajVot, ClassCMajVot)
#create bar plot
install.packages("tidyverse")
library(tidyverse)
ggplot(data=MajVotingFreq, aes(x=freq, y=types, fill = TumorType))+
geom_bar(stat="identity", width=0.5, position=position_dodge())+
theme(axis.text.x = element_text(size = 5))+
  xlab("Frequency")+
  ylab("Cell Typist Annotations")

####################################AddModuleScore
#Call in PanImmune gene panel
PanGenes <- read.csv("~/Desktop/PatelLab/PanImmuneGeneSet/PanGenes.csv")
#AddModuleScore

sets<-levels(factor(PanGenes$SetName))

#check which gene sets have 0 intersecting genes with Seuratobj
genecountperset<-1:length(sets)
for(x in 1:length(sets)){
  SetofGenes<-PanGenes[PanGenes$SetName==sets[x],2]
  intersection<-intersect(row.names(seuratObj),SetofGenes)
  genecountperset[x]<-length(intersection)
}


##gene sets with no matching genes "CD103pos_CD103neg_ratio_25446897" "GP11_Immune_IFN"   
sets<-sets[-which(sets=="CD103pos_CD103neg_ratio_25446897")]
sets<-sets[-which(sets=="GP11_Immune_IFN")]

###Create data frame with average module scores 
#List with all gene sets
for(set in 1:length(sets)){
  if(set==1){
    SetofGenes<-list(PanGenes[PanGenes$SetName==sets[set],2])
  } else {
    newset<-list(PanGenes[PanGenes$SetName==sets[set],2])
    SetofGenes<-append(SetofGenes, newset)  }
}

#seuratObj AddModuleScore
immune<-Seurat::AddModuleScore(immune, features = SetofGenes, name = "GeneSet")
metadata<-immune@meta.data
ModuleScoreIndex<-which(colnames(metadata)=="GeneSet1"):122

##create boxplot with module scores and Cell typist annotations'
#Class A
classA<-subset(immune, subset = TumorType == "ClassA")
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
classB<-subset(immune, subset = TumorType == "ClassB")
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
classC<-subset(immune, subset = TumorType == "ClassC")
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
#make ggplot
ggplot(data=ClassCScore, aes(x=score, y=setname, fill = TumorType))+
  geom_boxplot(outlier.size=0.5)+
  theme(axis.text.x = element_text(size = 5))+
  xlab("Module Score")+
  ylab("Pan Gene Set")+
  theme(axis.text = element_text(size = 5))+xlim(-1.5,4.5)


##################PCA analysis
mdata<-immune@meta.data
row.names(mdata)
immunepca<-prcomp(mdata[,ModuleScoreIndex], center = TRUE, scale. = TRUE)
install.packages("devtools")
library(devtools)
install.packages("ggfortify")
library(ggfortify)
autoplot(immunepca, data = mdata, colour = 'TumorType', )
summary(immunepca)


