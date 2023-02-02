##REMOVE Tumor/Non-immunecell clusters
#remove tumor clusters 1, 3, 10, 15, 24, 8, 12, 16, 19
#load in updated seurat Obj
load("~/Desktop/PatelLab/102022_SCELL_UCSF_BCM_WITH_newannotations_addmodulescore_seuratObj.rda")
sub<- subset(seuratObj, cells=colnames(seuratObj)[seuratObj$orig.ident %in% c("BCM_Front", "BCM_Post", "MSC3", "MSC5", "MSC1", "MSC4")])

sub$TumorType <- rep("", length(sub$seurat_clusters))

sub$TumorType[sub$orig.ident %in% c("BCM_Front", "BCM_Post", "MSC3")] <- "ClassC"

sub$TumorType[sub$orig.ident %in% c("MSC5")] <- "ClassB"

sub$TumorType[sub$orig.ident %in% c("MSC1", "MSC4")] <- "ClassA"

#subset only immune tumors
clusters<-levels(sub$seurat_clusters)
clusters<-clusters[c(-2,-4,-11,-16,-25, -9, -13, -17, -20)]
immune<- subset(sub, cells=colnames(sub)[sub$seurat_clusters %in% clusters])
save(immune, file = "~/Desktop/PatelLab/notumorseuratObj.Rdata")

#Load in nontumor cell clusters
load("~/Desktop/PatelLab/notumorseuratObj.Rdata")

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

######create boxplot with module scores and Cell typist annotations'
#load in gene set
PanGenes <- read.csv("~/Desktop/PatelLab/PanImmuneGeneSet/PanGenes.csv")
sets<-levels(factor(PanGenes$SetName))
sets<-sets[-which(sets=="CD103pos_CD103neg_ratio_25446897")]
sets<-sets[-which(sets=="GP11_Immune_IFN")]


#seuratObj AddModuleScore
metadata<-immune@meta.data
ModuleScoreIndex<-which(colnames(metadata)=="GeneSet1"):122

##Split gene sets into three groups in order for cleaner visualization
sets<-sets[1:36]
ModuleScoreIndex<-ModuleScoreIndex[1:36]

subset1<-sets[1:36]
subset2<-sets[37:72]
subset3<-sets[73:107]

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
mdata<-immune@meta.data
row.names(mdata)
immunepca<-prcomp(mdata[,ModuleScoreIndex], center = TRUE, scale. = TRUE)
install.packages("devtools")
library(devtools)
install.packages("ggfortify")
library(ggfortify)
autoplot(immunepca, data = mdata, colour = 'TumorType', )
summary(immunepca)


