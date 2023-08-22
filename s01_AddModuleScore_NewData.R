#ADD MODULE SCORE WORKSPACE
library(Seurat)

setwd("~/Desktop/PatelLab/NewAnalysis")

load("~/Desktop/PatelLab/060223_Heidelberg_SCELL_seuratObj.rda")


######filter out non-immune cells
table(seuratObj$seurat_clusters)
#Immune clusters are "13", "3", "18", "16". 18/3/16 are macrophages/microglias. Cluster 13 is T-cells.
sub<- subset(seuratObj, cells=colnames(seuratObj)[seuratObj$seurat_clusters %in% c("13", "3", "18", "16")])

#save as rda file
save(sub, file = "072523_Heidellberg_SCELL_seuratObj_immunecellsonly.rda")

#Call in PanImmune gene panel
PanGenes <- read.csv("~/Desktop/PatelLab/PanImmuneGeneSet/PanGenes.csv")
#AddModuleScore

sets<-levels(factor(PanGenes$SetName))

#check which gene sets have 0 intersecting genes with Seuratobj
genecountperset<-1:length(sets)
for(x in 1:length(sets)){
  SetofGenes<-PanGenes[PanGenes$SetName==sets[x],2]
  intersection<-intersect(row.names(sub),SetofGenes)
  genecountperset[x]<-length(intersection)
}

sets[genecountperset==0]

write.csv(sets, file = "genesets_PanImmunegenes.csv")
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

#Load in seuratObj
#load("~/Desktop/PatelLab/091622_SCELL_UCSF_BCM_seuratObj_ann.rda")

#seuratObj AddModuleScore
seuratObj<-Seurat::AddModuleScore(sub, features = SetofGenes, name = "GeneSet")
SingleCellTypeSummary <- read.csv("~/Desktop/PatelLab/NewAnalysis/SingleCellTypeSummary.csv", row.names=1)

SingleCellTypeSummary$Type

TumorType<- rep("", length(seuratObj$seurat_clusters))
ident<-seuratObj$orig.ident

for(x in 1:length(ident)){
  index<-row.names(SingleCellTypeSummary)%in%ident[x]
  TumorType[x]<-SingleCellTypeSummary[index,1]
}

seuratObj$TumorType<-TumorType


#add cell typist annotations
load("072523_allLocations_immunecellsonly_modulescores.rda")
predicted_labels <- read.csv("~/Desktop/PatelLab/NewAnalysis/predicted_labels.csv", row.names=1)

seuratObj$celltypist<-predicted_labels$majority_voting


#remove cells with NA class
seuratObj<-subset(seuratObj, cells=colnames(seuratObj)[!is.na(seuratObj$TumorType)])


save(seuratObj, file = "072523_allLocations_immunecellsonly_modulescores.rda")






#Old ks test. Can skip this step
#ClassAModuleScores <- read.csv("~/Desktop/PatelLab/Analysis_Results/ClassAModuleScores.csv", row.names=1)
#ClassBModuleScores <- read.csv("~/Desktop/PatelLab/Analysis_Results/ClassBModuleScores.csv", row.names=1)
#ClassCModuleScores <- read.csv("~/Desktop/PatelLab/Analysis_Results/ClassCModuleScores.csv", row.names=1)
#kstestAB<-ks.test(ClassAModuleScores[,1],ClassBModuleScores[,1])
###Compare gene sets among different classes
for(x in 1:length(sets)){
  if(x == 1){
  kstestAB<-ks.test(ClassAModuleScores[,x],ClassBModuleScores[,x])
  kstestAC<-ks.test(ClassAModuleScores[,x],ClassCModuleScores[,x])
  kstestBC<-ks.test(ClassBModuleScores[,x],ClassCModuleScores[,x])
  kstestanalysis<-data.frame(kstestAB$p.value, kstestAC$p.value, kstestBC$p.value)
  } else {
    kstestAB<-ks.test(ClassAModuleScores[,x],ClassBModuleScores[,x])
    kstestAC<-ks.test(ClassAModuleScores[,x],ClassCModuleScores[,x])
    kstestBC<-ks.test(ClassBModuleScores[,x],ClassCModuleScores[,x])
    newkstestanalysis<-data.frame(kstestAB$p.value, kstestAC$p.value, kstestBC$p.value)
    kstestanalysis<-rbind(kstestanalysis,newkstestanalysis)
  }
}

row.names(kstestanalysis)<-sets
row.names(kstestDvalue)<-sets
colnames(kstestanalysis)<-c("A vs B", "A vs C","B vs C")

interestinggenesets<-which(kstestanalysis$`A vs B`!=0|kstestanalysis$`A vs C`!=0|kstestanalysis$`B vs C`!=0)
interestgenesetkstest<-kstestanalysis[interestinggenesets,]

write.csv(kstestanalysis, "~/Desktop/PatelLab/kstestanalysis_pvalues.csv")
write.csv(interestgenesetkstest, "~/Desktop/PatelLab/interestsets_kstest_pvalues.csv")


