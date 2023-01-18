###Deconvolution Workspace
#How many genes are in each module
RepresentativeGenesets <- read.csv("~/Desktop/PatelLab/Analysis_Results/Deconvolution/RepresentativeGenesets.csv", header=FALSE)
colnames(RepresentativeGenesets)<-c("Gene_Sets","Module")

#Call in PanImmune gene panel
PanGenes <- read.csv("~/Desktop/PatelLab/PanImmuneGeneSet/PanGenes.csv")
genesets<-RepresentativeGenesets$Gene_Sets
sigsets<-length(RepresentativeGenesets$Gene_Sets)
for(x in 1:sigsets){
  if(x == 1){
    num<-length(PanGenes[PanGenes$SetName==genesets[x],2])
    mod<-RepresentativeGenesets[x,2]
    Genesinmodules<-data.frame(PanGenes[PanGenes$SetName==genesets[x],2],rep(genesets[x], num), rep(mod, num))
  } else {
    num<-length(PanGenes[PanGenes$SetName==genesets[x],2])
    mod<-RepresentativeGenesets[x,2]
    new<-data.frame(PanGenes[PanGenes$SetName==genesets[x],2],rep(genesets[x], num), rep(mod, num))
    Genesinmodules<-rbind(Genesinmodules,new)
  }
}
colnames(Genesinmodules)<-c("Genes","Genesets","Module")
table(Genesinmodules$Module)
length(unique(Genesinmodules$Genes))
length(unique(PanGenes$Gene))

test<-PanGenes[PanGenes$Gene=="CTSS",]
intersect(test$SetName,Genesinmodules$Genesets)

#create dataframe with unique genes note: some modules had repeats of genes
blue<-Genesinmodules[Genesinmodules$Module == "Blue",]
brown<-Genesinmodules[Genesinmodules$Module == "Brown",]
turq<-Genesinmodules[Genesinmodules$Module == "Turquoise",]
yellow<-Genesinmodules[Genesinmodules$Module == "Yellow",]
green<-Genesinmodules[Genesinmodules$Module == "Green",]
length(unique(blue$Genes))
length(unique(brown$Genes))
length(unique(turq$Genes))
length(unique(yellow$Genes))
length(unique(green$Genes))

###Bulk Sequence
GS<-row.names(data_GS)
intergenes<-intersect(GS, Genesinmodules$Genes)
#which(row.names(normalized_data_GS)==intergenes)
genename_index<-""
for(x in 1:length(intergenes)){
  genename_index<-c(genename_index,which(row.names(normalized_data_GS)==intergenes[x]))
}
genename_index<-genename_index[-1]
genename_index<-as.numeric(genename_index)

modulegenes_data<-normalized_data_GS[genename_index,]

#Add module column to modulegenes_data
genename_index<-""
for(x in 1:length(intergenes)){
  genename_index<-c(genename_index,which(Genesinmodules$Genes==intergenes[x]))
}
genename_index<-genename_index[-1]
genename_index<-as.numeric(genename_index)

intersectinggenes<-Genesinmodules[genename_index,]
write.table(table(intersectinggenes$Genes,intersectinggenes$Module), file = "~/Desktop/PatelLab/table_intersectinggenesbymodule.txt")

##Extract modules for genes that don't have overlap across modules and load in genes that have multiple modules. Combine the two data 
#frames
#Load in genes that don't overlap across modules (this was manually done using excel)
nonmultimodule_genes <- read.csv("Desktop/PatelLab/Analysis_Results/Deconvolution/nonmultimodule_genes.csv")
nonmultigenes<-nonmultimodule_genes$Genes

for(x in 1:length(nonmultigenes)){
  if(x == 1){
    nonmulti<-intersectinggenes[intersectinggenes$Genes==nonmultigenes[x],]
  }else{
  sset<-intersectinggenes[intersectinggenes$Genes==nonmultigenes[x],]
  nonmulti<-rbind(nonmulti,sset)
  }
}

multimodule_genes <- read.csv("Desktop/PatelLab/Analysis_Results/Deconvolution/multimodule_genes.csv")
nonmulti<-nonmulti[!duplicated(nonmulti$Genes),]
nonmulti<-nonmulti[,-2]

multi<-data.frame(multimodule_genes,rep("Multiple",length(multimodule_genes)))
colnames(multi)<-c("Genes", "Module")
modules<-rbind(nonmulti,multi)
row.names(modules)<-modules$Genes
modules<-modules[-1]


#Extract class data for samples
bulkdataclass<-data.frame(samples$class)
row.names(bulkdataclass)<-row.names(samples)
colnames(bulkdataclass)<-"TumorClass"

#Heat Map
library(pheatmap)
install.packages("RColorBrewer")
library(RColorBrewer)
breaks <- seq(-2, 2, length = 16)
color <- colorRampPalette(rev(brewer.pal(11,'RdBu')))(length(breaks))
pheatmap(modulegenes_data, annotation_col = bulkdataclass, annotation_row = modules, show_rownames = F, show_colnames = F, scale = "row",
         color= color,breaks = breaks)

#load in single cell data after WGCNA analysis
load("~/Desktop/PatelLab/Analysis_Results/Deconvolution/seurat_significantimmunecells.rda")

#Load in single cell data matrix
immunecellscounts <- as.matrix(Seurat::GetAssayData(sigimmune, slot = "counts"))

#cell type labels vector
celltype<-sigimmune$TumorType

#cell state label vector
cellstate<-sigimmune$representativemodule

table(cbind.data.frame(cellstate, celltype))

#save bulk data to input into Cibersort
write.table(data_GS, file = "~/Desktop/PatelLab/Analysis_Results/Deconvolution/bulkdata_unnormalized.txt", quote = F, sep = "\t")

##Subset the single cell data by the genes in the modules
subsetgenes<-unique(Genesinmodules$Genes)
immunecellscounts <- as.matrix(Seurat::GetAssayData(sigimmune, slot = "counts"))
#create index of genes to subset matrix
geneindex <- ""
matrixgenes<-row.names(immunecellscounts)
length<-length(subsetgenes)
for(x in 1:length){
  geneindex<-c(geneindex,which(matrixgenes==subsetgenes[x]))
}
geneindex<-geneindex[-1]
geneindex<-as.numeric(geneindex)
subset.matrix <- immunecellscounts[geneindex, ] # Pull the raw expression matrix from the original Seurat object containing only the genes of interest
colnames(subset.matrix)<-sigimmune$representativemodule
#write table to input into cibersort
write.table(subset.matrix, file = "~/Desktop/PatelLab/Analysis_Results/Deconvolution/subset.matrix.Cibersort.txt", quote = F, sep = "\t")

#create a pheatmap with both LM22 and Module ratios
Cibserort_deconvolution_LM22_bulk <- read.csv("Desktop/PatelLab/Analysis_Results/Deconvolution/Cibserort_deconvolution_LM22_bulk.csv")
Cibersort_Module_deconvolution <- read.csv("Desktop/PatelLab/Analysis_Results/Deconvolution/Cibersort_Module_deconvolution.csv")

#combine the two dataframes
Cibersort_Module_deconvolution<-Cibersort_Module_deconvolution[,c(-7,-8,-9)]
Cibserort_deconvolution_LM22_bulk<-Cibserort_deconvolution_LM22_bulk[,c(-24,-25,-26)]
deconv<-data.frame(Cibersort_Module_deconvolution,Cibserort_deconvolution_LM22_bulk)
deconv<-deconv[,-7]
row.names(deconv)<-deconv$Mixture
deconv<-deconv[,-1]
deconv<-t(deconv)
deconv<-as.data.frame(deconv)

#vector with tumor class
sampleclass<-as.data.frame(samples$class)
row.names(sampleclass)<-row.names(samples)
colnames(sampleclass)<-"TumorClass"

#create pheatmap
library(pheatmap)
library(RColorBrewer)
pheatmap(deconv, annotation_col = sampleclass, show_colnames = F, cluster_rows = F, scale = "row",
         color= color,breaks = breaks)
ncol(deconv)
length(sampleclass)

##correlational matrix between LM22 and Module
Moduledeconv<-Cibersort_Module_deconvolution[,-1]
row.names(Moduledeconv)<-Cibersort_Module_deconvolution$Mixture
LM22deconv<-Cibserort_deconvolution_LM22_bulk[-1]
row.names(LM22deconv)<-Cibserort_deconvolution_LM22_bulk$Mixture
cormatrix<-cor(Moduledeconv,LM22deconv)
cormatrix<-cormatrix[,-7]
pheatmap(cormatrix)
write.table(cormatrix, file = '~/Desktop/PatelLab/Analysis_Results/Deconvolution/cormatrix.txt', quote = F, sep = "\t")
