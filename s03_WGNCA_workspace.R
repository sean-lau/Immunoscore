##WGCNA
#Input Module Scores from each gene set into WGCNA to identify relation among the gene sets
library(Seurat)
install.packages("BiocManager")
BiocManager::install("WGCNA")
library(WGCNA)
options(stringsAsFactors = FALSE)


load("~/Desktop/PatelLab/notumorseuratObj.Rdata")
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

##Use 9 as the soft threshold power

net = blockwiseModules(metadf, maxBlockSize = 107,
                         power = 9, TOMType = "signed", minModuleSize = 3, pamStage = F,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "immunecellTOM-blockwise",
                         verbose = 3)
save(net, file = "~/Desktop/PatelLab/Analysis_Results/WGCNA/wgcna_results.rda")

### Shortcut: load WGCNA results from above analysis
load("~/Desktop/PatelLab/Analysis_Results/WGCNA/wgcna_results.rda")

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

##assign single cells to module
load("~/Desktop/PatelLab/Analysis_Results/WGCNA/wgcna_results.rda")
print(moduleColors)
print(moduleLabels)
eigensig<-net$MEs
eigensig<-eigensig[-6]
MEcolor<-c("blue","yellow","turquoise","brown","green")
colnames(eigensig)<-MEcolor

##use for loop to assign each cell a module by finding the module with the highest score for each cell.

cellnum<-nrow(eigensig)
cellmodule_WGCNA<-""
for(x in 1:cellnum){
  cell<-eigensig[x,]
  if(cell[1]<0&cell[2]<0&cell[3]<0&cell[4]<0&cell[5]<0){
    cellmodule_WGCNA<-c(cellmodule_WGCNA, "N")
  } else {
  cellmodule_WGCNA<-c(cellmodule_WGCNA,colnames(cell[which.max(cell)]))
  }
}
cellmodule_WGCNA<-cellmodule_WGCNA[-1]
immune$representativemodule<-cellmodule_WGCNA
sigimmune<-subset(immune, cells=colnames(immune)[immune$representativemodule != "N"])
ncol(sigimmune)
##1999 cells with only negative eigen signature scores removed

save(sigimmune,file = "~/Desktop/PatelLab/Analysis_Results/Deconvolution/seurat_significantimmunecells.rda")

##Develop frequency plot with the number of modules per 
load("~/Desktop/PatelLab/Analysis_Results/Deconvolution/seurat_significantimmunecells.rda")
metadata<-sigimmune@meta.data
metadata<-metadata[,123:124]
table(metadata$TumorType, metadata$representativemodule)

##ggplot
metadata$representativemodule<-gsub("turquoise", "C1", metadata$representativemodule)
metadata$representativemodule<-gsub("blue", "C2", metadata$representativemodule)
metadata$representativemodule<-gsub("green", "C3", metadata$representativemodule)
metadata$representativemodule<-gsub("yellow", "C4", metadata$representativemodule)
metadata$representativemodule<-gsub("brown", "C5", metadata$representativemodule)

#create data frame with frequencies for each variation

ClassA<-metadata[metadata$TumorType=="ClassA",]
ClassAfreq<-table(ClassA$representativemodule)
ClassAfreq<-as.data.frame(ClassAfreq)
colnames(ClassAfreq)<-c("Module", "Frequency")
ClassAfreq$Frequency<-ClassAfreq$Frequency/nrow(ClassA)
ClassAfreq$Class<-rep("A", nrow(ClassAfreq))

ClassB<-metadata[metadata$TumorType=="ClassB",]
ClassBfreq<-table(ClassB$representativemodule)
ClassBfreq<-as.data.frame(ClassBfreq)
colnames(ClassBfreq)<-c("Module", "Frequency")
ClassBfreq$Frequency<-ClassBfreq$Frequency/nrow(ClassB)
ClassBfreq$Class<-rep("B", nrow(ClassBfreq))

ClassC<-metadata[metadata$TumorType=="ClassC",]
ClassCfreq<-table(ClassC$representativemodule)
ClassCfreq<-as.data.frame(ClassCfreq)
colnames(ClassCfreq)<-c("Module", "Frequency")
ClassCfreq$Frequency<-ClassCfreq$Frequency/nrow(ClassC)
ClassCfreq$Class<-rep("C", nrow(ClassCfreq))

#Combine frequencies tables

combine<-rbind(ClassAfreq, ClassBfreq, ClassCfreq)

b<-ggplot(data = combine, mapping = aes(x = Module, y = Frequency,  fill = Class))+
  geom_bar(stat = "identity", position=position_dodge())+
  scale_fill_manual(values=c('#28A860','#0042ED', "#FE4733"))+
  labs(title="Frequency of Modules by Tumor Class", 
       x="Module", y = "Relative Frequency", fill = "Tumor Class")+
  ylim(0,1)
ggsave("ModuleFrequency_barplot_bytumorclass.pdf",  plot = b, width = 10, height = 5)
