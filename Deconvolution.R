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
GS<-row.names(normalized_data_GS)
intergenes<-intersect(GS, Genesinmodules$Genes)
which(row.names(normalized_data_GS)==intergenes)
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

##Extract modules for genes that don't have overlap across modules and load in genes that have multiple modules. Combine the two data frames
nonmultigenes<-nonmultimodule_genes$Genes

for(x in 1:length(nonmultigenes)){
  if(x == 1){
    nonmulti<-intersectinggenes[intersectinggenes$Genes==nonmultigenes[x],]
  }else{
  sset<-intersectinggenes[intersectinggenes$Genes==nonmultigenes[x],]
  nonmulti<-rbind(nonmulti,sset)
  }
}

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
pheatmap(modulegenes_data, annotation_col = bulkdataclass, annotation_row = modules, show_rownames = F, show_colnames = F, scale = "row")

