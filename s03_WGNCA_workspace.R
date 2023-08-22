##WGCNA
#Input Module Scores from each gene set into WGCNA to identify relation among the gene sets
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

##Use 5 as the soft threshold power

net = blockwiseModules(metadf, maxBlockSize = 107,
                         power = 5, TOMType = "signed", minModuleSize = 3, pamStage = F,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "immunecellTOM-blockwise",
                         verbose = 3)
save(net, file = "~/Desktop/PatelLab/NewAnalysis/WGCNA/wgcna_results.rda")
