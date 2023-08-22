##WGNCA Analysis

##load in necessary packages and data
library(Seurat)
library(WGCNA)
options(stringsAsFactors = FALSE)

setwd("~/Desktop/PatelLab/NewAnalysis")
load("072523_allLocations_immunecellsonly_modulescores.rda")
#Call in PanImmune gene panel
PanGenes <- read.csv("~/Desktop/PatelLab/PanImmuneGeneSet/PanGenes.csv")
sets<-levels(factor(PanGenes$SetName))

##gene sets with no matching genes "CD103pos_CD103neg_ratio_25446897" "GP11_Immune_IFN"   
sets<-sets[-which(sets=="CD103pos_CD103neg_ratio_25446897")]
sets<-sets[-which(sets=="GP11_Immune_IFN")]

immune<-seuratObj
metadata<-immune@meta.data
ModuleScoreIndex<-which(colnames(metadata)=="GeneSet1"):114
#Matrix of data
metadf<-metadata[,ModuleScoreIndex]
colnames(metadf)<-sets

#load in wgcna results
load("~/Desktop/PatelLab/NewAnalysis/WGCNA/wgcna_results.rda")

#visualizisation
moduleLabels = net$colors


moduleColors = labels2colors(net$colors)
names(moduleColors)<-names(moduleLabels)

sizeGrWindow(6,6)

plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors", main = "Cluster Dendrogram",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

#Tumor Type data
samples<-rownames(metadf)
matchingrows<-match(samples,rownames(metadata))
datTrait<-metadata[matchingrows,115]

####identifying relationships among genesets within modules and modules with tumor class
#correlational analysis of modules within each tumor class
MEs = net$MEs;
geneTree = net$dendrograms[[1]]

# Define numbers of genes and samples
nGenes = ncol(metadf);
nSamples = nrow(metadf);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(metadf, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

net$colors

classeigendf<-data.frame(MEs,datTrait)



#need to create dataframe with mean expression values of each module separated by tumor class
#need visual representation of the mean scores
library(ggplot2)
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE)/sqrt(length(x[[col]])))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

##classifications: MEgreen:C5, MEbrown:C3, MEturquoise:C1, C2:MEblue, C4:MEyellow, C0:MEgrey
library(plyr)
sumdata<-data_summary(classeigendf, varname = "MEgrey", groupnames = "datTrait" )

ggplot(data = sumdata, aes(x = datTrait, y = MEgrey))+
  geom_bar(stat="identity", fill = "lightblue",
           position=position_dodge())+
  geom_errorbar(aes(ymin=MEgrey-sd*2, ymax=MEgrey+sd*2), width=.2,
                position=position_dodge(.9))+
  ggtitle("Average WGCNA Module Score by Tumor Class")+ xlab("Tumor Class")+
  ylab("C0 Module Eigen Score")



#need to ask Dr. Harmanci about which regression models might work best
#need to figure out gene signficance of each gene set in each module
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(metadf, MEs, use = "p"));
colnames(geneModuleMembership)<-c("C5", "C3", "C1", "C2", "C4", "C0")
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

library(pheatmap)

pheatmap(geneModuleMembership, fontsize = 6)
write.table(geneModuleMembership, file = '~/Desktop/PatelLab/NewAnalysis/geneModuleMembership.csv')


#ncol(MMPvalue)

#names(geneModuleMembership) = paste("MM", modNames, sep="");
#names(MMPvalue) = paste("p.MM", modNames, sep="")

#pheatmap(MMPvalue)
#rownames(geneModuleMembership)<- gsub(" ", "_",rownames(geneModuleMembership))



#create data frame with association between tumor class and module

classA<-classeigendf[classeigendf$datTrait == "A",]
Atvalue<-""
for(x in 1:6){
  t<-t.test(classA[,x])
  Atvalue<-c(Atvalue,t$statistic)
}
Atvalue<-Atvalue[-1]
Atvalue<-as.numeric(Atvalue)

classB<-classeigendf[classeigendf$datTrait == "B",]
Btvalue<-""
for(x in 1:6){
  t<-t.test(classB[,x])
  Btvalue<-c(Btvalue,t$statistic)
}
Btvalue<-Btvalue[-1]
Btvalue<-as.numeric(Btvalue)

classC<-classeigendf[classeigendf$datTrait == "C",]
Ctvalue<-""
for(x in 1:6){
  t<-t.test(classC[,x])
  Ctvalue<-c(Ctvalue,t$statistic)
}
Ctvalue<-Ctvalue[-1]
Ctvalue<-as.numeric(Ctvalue)

classmoduleassoctable<-data.frame("ClassA" = Atvalue,"ClassB"= Btvalue,"ClassC" = Ctvalue)
modcolor<-names(classeigendf)
modcolor<-modcolor[1:6]
row.names(classmoduleassoctable)<-modcolor
row.names(classmoduleassoctable)<-c("C5", "C3", "C1", "C2", "C4", "C0")
write.table(classmoduleassoctable, file = "~/Desktop/PatelLab/NewAnalysis/WGCNA/classbymodule_tvaluetable.csv")


##########
##########Create pheatmap of module scores from all genesets with WGCNA modules as annotations
######################
library(WGCNA)
load("072523_allLocations_immunecellsonly_modulescores.rda")

#load in WGCNA results
load("~/Desktop/PatelLab/NewAnalysis/WGCNA/wgcna_results.rda")

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


metadata<-seuratObj@meta.data
ModuleScoreIndex<-which(colnames(metadata)=="GeneSet1"):114
#Matrix of data
metadf<-metadata[,ModuleScoreIndex]
colnames(metadf)<-sets

#correlational matrix of scores
Figure1A<-cor(metadf)
#add module labels
##classifications: MEgreen:C5, MEbrown:C3, MEturquoise:C1, C2:MEblue, C4:MEyellow, C0:MEgrey

moduleColors<-as.data.frame(moduleColors)
moduleColors$moduleColors<-gsub("turquoise", "C1",moduleColors$moduleColors)
moduleColors$moduleColors<-gsub("blue", "C2", moduleColors$moduleColors)
moduleColors$moduleColors<-gsub("green", "C5", moduleColors$moduleColors)
moduleColors$moduleColors<-gsub("yellow", "C4", moduleColors$moduleColors)
moduleColors$moduleColors<-gsub("brown", "C3", moduleColors$moduleColors)
moduleColors$moduleColors<-gsub("grey", "C0", moduleColors$moduleColors)

library(pheatmap)
p<-pheatmap(Figure1A, fontsize = 7 , show_rownames = T, show_colnames = F, annotation_row = moduleColors)
ggsave("heatmap_cormatrix_genesignatures_Figure1A.pdf", plot = p, width = 14, height = 10)

