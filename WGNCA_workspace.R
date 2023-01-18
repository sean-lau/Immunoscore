##December Work
#Input Module Scores from each gene set into WGCNA to identify relation among the gene sets
install.packages("Seurat")
library(Seurat)
install.packages("BiocManager")
BiocManager::install("WGCNA")
library(WGCNA)
options(stringsAsFactors = FALSE)

row.names(immune)

load("~/Desktop/PatelLab/notumorseuratObj.Rdata")
#Call in PanImmune gene panel
PanGenes <- read.csv("~/Desktop/PatelLab/PanImmuneGeneSet/PanGenes.csv")
#AddModuleScore

sets<-levels(factor(PanGenes$SetName))

##gene sets with no matching genes "CD103pos_CD103neg_ratio_25446897" "GP11_Immune_IFN"   
sets<-sets[-which(sets=="CD103pos_CD103neg_ratio_25446897")]
sets<-sets[-which(sets=="GP11_Immune_IFN")]


metadata<-immune@meta.data
ModuleScoreIndex<-which(colnames(metadata)=="GeneSet1"):122
#Matrix of data
metadf<-metadata[,ModuleScoreIndex]
colnames(metadf)<-sets

##Data Cleaning
gsg <- goodSamplesGenes(metadf, verbose = 3)
gsg$allOK

##automatic network construction and module detection
#choosing a set of soft-thresholdinng powers
powers <- c(c(1:10), seq(from = 1, to = 20, by = 2))
#call newtork topology analysis function
sft <-pickSoftThreshold(metadf, powerVector = powers, verbose = 5)
###plot results
sizeGrWindow(9,5)
par(mfrow = c(1,2))
cex1 <- 0.9

#Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshld (power)", ylab = "Scale Free topology model Fit, signed R^2", type= "n", main = paste("scale independence"), ylim = c(0,1))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

abline(h=0.90,col="red")


plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
?blockwiseModules()

##Use 9 as the soft threshold power

net = blockwiseModules(metadf, maxBlockSize = 107,
                         power = 9, TOMType = "signed", minModuleSize = 3, pamStage = F,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "immunecellTOM-blockwise",
                         verbose = 3)
#visualizisation
moduleLabels = net$colors


moduleColors = labels2colors(net$colors)
sizeGrWindow(6,6)

plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors", main = "Cluster Dendrogram",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

#Tumor Type data
samples<-rownames(metadf)
matchingrows<-match(samples,rownames(metadata))
datTrait<-metadata[matchingrows,123]

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

sumdata<-data_summary(classeigendf, varname = "MEgreen", groupnames = "datTrait" )

ggplot(data = sumdata, aes(x = datTrait, y = MEgreen))+
  geom_bar(stat="identity", fill = "lightblue",
           position=position_dodge())+
  geom_errorbar(aes(ymin=MEgreen-sd*2, ymax=MEgreen+sd*2), width=.2,
                position=position_dodge(.9))+
  ggtitle("Average Module Score by Tumor Class")+ xlab("Tumor Class")+
  ylab("Green Module Eigen Score")

#need to ask Dr. Harmanci about which regression models might work best
#need to figure out gene signficance of each gene set in each module
?substring()
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(metadf, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
library(pheatmap)
pheatmap(geneModuleMembership, fontsize = 6)
ncol(MMPvalue)

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="")

pheatmap(MMPvalue)
rownames(geneModuleMembership)<- gsub(" ", "_",rownames(geneModuleMembership))

write.table(geneModuleMembership, file = '~/Desktop/PatelLab/geneModuleMembership.csv')


#create data frame with association between tumor class and module

classA<-classeigendf[classeigendf$datTrait == "ClassA",]
Atvalue<-""
for(x in 1:6){
  t<-t.test(classA[,x])
  Atvalue<-c(Atvalue,t$statistic)
}
Atvalue<-Atvalue[-1]
Atvalue<-as.numeric(Atvalue)

classB<-classeigendf[classeigendf$datTrait == "ClassB",]
Btvalue<-""
for(x in 1:6){
  t<-t.test(classB[,x])
  Btvalue<-c(Btvalue,t$statistic)
}
Btvalue<-Btvalue[-1]
Btvalue<-as.numeric(Btvalue)

classC<-classeigendf[classeigendf$datTrait == "ClassC",]
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
write.table(classmoduleassoctable, file = "~/Desktop/PatelLab/Analysis_Results/WGCNA/classbymodule_tvaluetable.csv")
