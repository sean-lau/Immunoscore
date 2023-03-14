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

###################Bulk Sequence############################################################
####################################################################################################
##load in BCM bulk sample data
load("~/Desktop/PatelLab/09192022_GE_RAW_Normalized_N330PrimarySamples_GeneSymbol.rda")
#hongkong data
load("~/Desktop/PatelLab/031023_Inhouse_Raleigh_Hongkong_GeneSymbol.rda")
data_GS<-data.frame(data_GS, data)
GS<-row.names(data_GS)

intergenes<-intersect(GS, Genesinmodules$Genes)
#which(row.names(normalized_data_GS)==intergenes)
genename_index<-""
for(x in 1:length(intergenes)){
  genename_index<-c(genename_index,which(row.names(data_GS)==intergenes[x]))
}
genename_index<-genename_index[-1]
genename_index<-as.numeric(genename_index)

modulegenes_data<-data_GS[genename_index,]

#Add module column to modulegenes_data
genename_index<-""
for(x in 1:length(intergenes)){
  genename_index<-c(genename_index,which(Genesinmodules$Genes==intergenes[x]))
}
genename_index<-genename_index[-1]
genename_index<-as.numeric(genename_index)

intersectinggenes<-Genesinmodules[genename_index,]

##Extract modules for genes that don't have overlap across modules and load in genes that have multiple modules. Combine the two data 
#frames
#Load in genes that don't overlap across modules (this was manually done using excel)
nonmultimodule_genes <- read.csv("~/Desktop/PatelLab/Analysis_Results/Deconvolution/nonmultimodule_genes.csv", sep="")
nonmultigenes<-nonmultimodule_genes$Genes

for(x in 1:length(nonmultigenes)){
  if(x == 1){
    nonmulti<-intersectinggenes[intersectinggenes$Genes==nonmultigenes[x],]
  }else{
    sset<-intersectinggenes[intersectinggenes$Genes==nonmultigenes[x],]
    nonmulti<-rbind(nonmulti,sset)
  }
}

multimodule_genes <- read.csv("~/Desktop/PatelLab/Analysis_Results/Deconvolution/multimodule_genes.csv", sep="")
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

#hongkong class data
hkclass<-samplesheet$class
hkclass<-as.data.frame(hkclass)
names(hkclass)<-"TumorClass"
row.names(hkclass)<-row.names(samplesheet)

#combine class data from both data sets
bulkdataclass<-rbind(bulkdataclass, hkclass)



#load in single cell data after subtype clustering
load("~/Desktop/PatelLab/Analysis_Results/Deconvolution/subtypecluster_cells.rda")

#Load in single cell data matrix
immunecellscounts <- as.matrix(Seurat::GetAssayData(subclustcells, slot = "counts"))


#cell type labels vector
#load in cluster data
celltype<-subclustcells$subtypecluster

#cell state label vector
cellstate<-subclustcells$TumorType


##Subset the single cell data by the genes in the modules. Not all genes in single cell data were 
#included in modules by WGCNA.
subsetgenes<-unique(Genesinmodules$Genes)
immunecellscounts <- as.matrix(Seurat::GetAssayData(subclustcells, slot = "counts"))
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


###Bayes Prism##########################################################################################
########################################################################################################
library(BayesPrism)

##need to split up cell states within the cell types

for(x in 1:length(cellstate)){
  cellstate[x]<-paste0(celltype[x],"_",cellstate[x])
}

table(cbind.data.frame(cellstate, celltype))


#transpose bulk and sc data
data_GS<-t(data_GS)
subset.matrix<-t(subset.matrix)


#identify whether cell states and types are of quality
plot.cor.phi(input=subset.matrix,
             input.labels=cellstate,
             title="cell state correlation",
             #specify pdf.prefix if need to output to pdf
             #pdf.prefix="gbm.cor.cs", 
             cexRow=0.2, cexCol=0.2,
             margins=c(2,2))

plot.cor.phi(input=subset.matrix, 
             input.labels=celltype, 
             title="cell type correlation",
             #specify pdf.prefix if need to output to pdf
             #pdf.prefix="gbm.cor.ct",
             cexRow=0.5, cexCol=0.5,
)

#filter outlier genes
sc.stat <- plot.scRNA.outlier(
  input=subset.matrix, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=celltype,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE #return the data used for plotting. 
  #pdf.prefix="gbm.sc.stat" specify pdf.prefix if need to output to pdf
)

bk.stat <- plot.bulk.outlier(
  bulk.input=data_GS,#make sure the colnames are gene symbol or ENSMEBL ID 
  sc.input=subset.matrix, #make sure the colnames are gene symbol or ENSMEBL ID 
  cell.type.labels=celltype,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE
  #pdf.prefix="gbm.bk.stat" specify pdf.prefix if need to output to pdf
)

sc.dat.filtered <- cleanup.genes (input=subset.matrix,
                                  input.type="count.matrix",
                                  species="hs", 
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1") ,
                                  exp.cells=5)

###number of genes filtered in each category: 
#Rb      Mrp other_Rb     chrM   MALAT1 
#0        0        0        0        0 
#A total of  0  genes from Rb Mrp other_Rb chrM MALAT1  have been excluded 
#A total of  9  gene expressed in fewer than  5  cells have been excluded 

#checking for concordance of expression for different types of genes
plot.bulk.vs.sc (sc.input = sc.dat.filtered,
                 bulk.input = data_GS
                 #pdf.prefix="gbm.bk.vs.sc" specify pdf.prefix if need to output to pdf
) ##may show an error due to "no complete element pairs" but this error shouldn't hinder the rest of the analysis

#mainly only see concordance among protein coding, so filter for protein coding genes
sc.dat.filtered.pc <-  select.gene.type (sc.dat.filtered,
                                         gene.type = "protein_coding")

#filter for significant genes to decrease cpu load
diff.exp.stat <- get.exp.stat(sc.dat=sc.dat.filtered.pc[,colSums(sc.dat.filtered.pc>0)>3],# filter genes to reduce memory use
                              cell.type.labels=celltype,
                              cell.state.labels=cellstate,
                              psuedo.count=0.1, #a numeric value used for log2 transformation. =0.1 for 10x data, =10 for smart-seq. Default=0.1.
                              cell.count.cutoff=50, # a numeric value to exclude cell state with number of cells fewer than this value for t test. Default=50.
                              n.cores=1 #number of threads
)

sc.dat.filtered.pc.sig <- select.marker(sc.dat=sc.dat.filtered.pc,
                                        stat=diff.exp.stat,
                                        pval.max=0.01,
                                        lfc.min=0.1)

#create bayes prism object
myPrism <- new.prism(
  reference=sc.dat.filtered.pc.sig, 
  mixture=data_GS,
  input.type="count.matrix", 
  cell.type.labels = celltype, 
  cell.state.labels = cellstate,
  key= NULL,
  outlier.cut=0.01,
  outlier.fraction=0.1,
)


#run bayes prism
bp.res <- run.prism(prism = myPrism, n.cores=50)



