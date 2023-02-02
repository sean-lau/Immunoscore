#ADD MODULE SCORE WORKSPACE
#AddModuleScore
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

sets[genecountperset==0]

write.csv(sets, "~/Desktop/PatelLab/genesets_PanImmunegenes.csv")
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
load("~/Desktop/PatelLab/091622_SCELL_UCSF_BCM_seuratObj_ann.rda")

#seuratObj AddModuleScore
seuratObj<-Seurat::AddModuleScore(seuratObj, features = SetofGenes, name = "GeneSet")

##new Seurat Obj with new annotations and module scores from this analysis saved 
#as 102022_SCELL_BCM_WITH_newannotation...

#Class A
#create data frame Genesets vs Cell barcodes
ClassA <- subset(sub, cells=colnames(sub)[sub$TumorType %in% "ClassA"])

ClassA<-Seurat::AddModuleScore(ClassA, features = SetofGenes)

MetaDataA<-ClassA@meta.data
ClassAScores<-MetaDataA[,which(colnames(MetaDataA)=="Cluster1"):121]
colnames(ClassAScores)<-sets
write.csv(ClassAScores,"~/Desktop/PatelLab/ClassAModuleScores.csv")
#Class B
ClassB <- subset(sub, cells=colnames(sub)[sub$TumorType %in% "ClassB"])

ClassB<-Seurat::AddModuleScore(ClassB, features = SetofGenes)

MetaDataB<-ClassB@meta.data
ClassBScores<-MetaDataB[,which(colnames(MetaDataB)=="Cluster1"):121]
colnames(ClassBScores)<-sets

write.csv(ClassBScores,"~/Desktop/PatelLab/ClassBModuleScores.csv")

#Class C
ClassC <- subset(sub, cells=colnames(sub)[sub$TumorType %in% "ClassC"])

ClassC<-Seurat::AddModuleScore(ClassC, features = SetofGenes)

MetaDataC<-ClassC@meta.data
ClassCScores<-MetaDataC[,which(colnames(MetaDataC)=="Cluster1"):121]
colnames(ClassCScores)<-sets

write.csv(ClassCScores,"~/Desktop/PatelLab/ClassCModuleScores.csv")


ClassAModuleScores <- read.csv("~/Desktop/PatelLab/Analysis_Results/ClassAModuleScores.csv", row.names=1)
ClassBModuleScores <- read.csv("~/Desktop/PatelLab/Analysis_Results/ClassBModuleScores.csv", row.names=1)
ClassCModuleScores <- read.csv("~/Desktop/PatelLab/Analysis_Results/ClassCModuleScores.csv", row.names=1)
kstestAB<-ks.test(ClassAModuleScores[,1],ClassBModuleScores[,1])
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

##Create pheatmap of module scores from all genesets with WGCNA modules as annotations
load("~/Desktop/PatelLab/notumorseuratObj.Rdata")
#load in WGCNA results
load("~/Desktop/PatelLab/Analysis_Results/WGCNA/wgcna_results.rda")

#visualizisation
moduleLabels = net$colors


moduleColors = labels2colors(net$colors)
names(moduleColors)<-names(moduleLabels)

#Call in PanImmune gene panel
PanGenes <- read.csv("~/Desktop/PatelLab/PanImmuneGeneSet/PanGenes.csv")
sets<-levels(factor(PanGenes$SetName))

##gene sets with no matching genes "CD103pos_CD103neg_ratio_25446897" "GP11_Immune_IFN"   
sets<-sets[-which(sets=="CD103pos_CD103neg_ratio_25446897")]
sets<-sets[-which(sets=="GP11_Immune_IFN")]


metadata<-immune@meta.data
ModuleScoreIndex<-which(colnames(metadata)=="GeneSet1"):122
#Matrix of data
metadf<-metadata[,ModuleScoreIndex]
colnames(metadf)<-sets

#correlational matrix of scores
Figure1A<-cor(metadf)
#add module labels
moduleColors<-as.data.frame(moduleColors)
moduleColors$moduleColors<-gsub("turquoise", "C1",moduleColors$moduleColors)
moduleColors$moduleColors<-gsub("blue", "C2", moduleColors$moduleColors)
moduleColors$moduleColors<-gsub("green", "C3", moduleColors$moduleColors)
moduleColors$moduleColors<-gsub("yellow", "C4", moduleColors$moduleColors)
moduleColors$moduleColors<-gsub("brown", "C5", moduleColors$moduleColors)
moduleColors$moduleColors<-gsub("grey", "C0", moduleColors$moduleColors)


p<-pheatmap(Figure1A, fontsize = 7 , show_rownames = T, show_colnames = F, annotation_row = moduleColors)
ggsave("heatmap_cormatrix_genesignatures_Figure1A.pdf", plot = p, width = 14, height = 10)

