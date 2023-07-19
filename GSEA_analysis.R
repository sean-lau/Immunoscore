#We can consider all the genes that are specifically up and down regulated in classes A/B/C. 
#Next we convert those to mouse gene symbol. We create a list of genes (A_Up, A_down, B_Up...). 
#And then use those gene  sets in GSEA analysis to see which ones are up/down in mouse. 
#In mouse you can perform mouse A vs C group comparison and give that to FGSEA R package
#create up and down regulated gene sets with additional tumor class subcategories

#intersect(siggenes,row.names(normalized_data_GS))
SignficantGenes <- read.csv("~/Desktop/PatelLab/akash_dog_meningioma/SignficantGenes.csv")
load("~/Desktop/PatelLab/akash_dog_meningioma/human_Rdata/09192022_GE_RAW_Normalized_N330PrimarySamples_GeneSymbol.rda")
Sampleclass<-data.frame("samplename"= row.names(samples),samples$class)
class<-c("A","B","C")
geneset<-""
for(a in 1:3){
  targetsamples<-Sampleclass[Sampleclass$samples.class==class[a],]
  targetdata<-normalized_data_GS[,colnames(normalized_data_GS)%in%targetsamples$samplename]
  referencesamples<-Sampleclass[Sampleclass$samples.class!=class[a],]
  referencedata<-normalized_data_GS[,colnames(normalized_data_GS)%in%referencesamples$samplename]
  siggenes<-SignficantGenes[SignficantGenes$Class==class[a],1]
 targetgenemeans<-apply(targetdata[row.names(targetdata)%in%siggenes,],1, function(x) {
   mean(x)})
 referencegenemeans<-apply(referencedata[row.names(referencedata)%in%siggenes,],1, function(x) {
   mean(x)})
 for(b in 1:length(targetgenemeans)){
   if(targetgenemeans[b]>referencegenemeans[b]){
     geneset<-c(geneset, paste0("UP_",class[a]))
   }else{
     geneset<-c(geneset, paste0("DOWN_", class[a]))
   }
 }
}

geneset<-geneset[-1]
table(geneset)

 #for(x in 1:length(siggenes)){
  #  if(mean(as.numeric(targetdata[which(row.names(targetdata)==siggenes[1]),]))>
  #     mean(referencedata[which(row.names(referencedata)==siggenes[x]),])){
   #   geneset<-c(geneset, print0("UP_",class[a]))
  #  }else (
  #    geneset<-c(geneset, print0("DOWN_", class[a]))
  #  )
#  }

geneset<-data.frame(SignficantGenes,geneset)


library(fgsea)
library(DESeq2)
library(limma)

#gsea analysis
#load in clusters
filtered_data_kmeans <- read.csv("~/Desktop/PatelLab/akash_dog_meningioma/filtered_data_kmeans.csv", 
                                 row.names=1)

dogclass<-filtered_data_kmeans[1:62,]

#load in dog data

setwd("~/Desktop/PatelLab/akash_dog_meningioma/")

#path <- "/Volumes/harmanci/baylor/akash/Canine_Meningioma/STAR_BAM/STAR_BAM"

data2 <- read.table("./data/counts.txt", header=T)
data <- data2[, -c(1:6)]
rownames(data) <- data2[,1]
colnames(data) <- gsub("_R1_001.fastq.gz_Aligned.sortedByCoord.out.bam", "", colnames(data))

geneConv <- read.delim("~/Desktop/PatelLab/akash_dog_meningioma/Canine_Human_Gene_Conversion.txt")
common <- intersect(rownames(data) , geneConv[,1] )
geneConv <- geneConv[match(common, geneConv[,1]),2]
data_mouse <- data[match(common, rownames(data)), ]
rownames(data_mouse) <- geneConv
data<-data_mouse
library(DESeq2)
library(limma)
coldata <- as.data.frame(rep(TRUE, each= length(colnames(data))))
rownames(coldata)<- colnames(data)
colnames(coldata)<- c("class")

coldata$class<-factor(dogclass$class)

dds<-DESeqDataSetFromMatrix(countData=(data),colData=coldata, ~ class)
dds <- DESeq(dds)
resultsNames(dds)

#load in dog data
res <- results(dds, contrast=c("class","dogA","dogC"), pAdjustMethod="fdr")
res$geneSymbol <- rownames(res)
res <- data.frame(res)

library(magrittr)
library(dplyr)
library(tibble)
res2 <- res %>%
  dplyr::select(geneSymbol, stat) %>%
  na.omit() %>%
  distinct() %>%
  group_by(geneSymbol) %>%
  summarize(stat=mean(stat))
ranks <- deframe(res2)


#create gene sets list
geneSets <- list()

# Loop over unique gene sets
uniqueGeneSets <- unique(geneset$geneset)
for (x in uniqueGeneSets) {
  # Subset the dataframe for the current gene set
  geneSet <- geneset[geneset$geneset == x,1]
  
  # Store the gene set as a named vector in the list
  geneSets[[x]] <- geneSet
}

# Convert the list to a named list
geneSets <- setNames(geneSets, uniqueGeneSets)
fgseaRes <- fgseaMultilevel(pathways = geneSets,
                            stats    = ranks)



#genesets is the human A_Up, A_down, B_Up... genesets

#Analysis of leading up_A genes between dogsA and dogsB
LeadingUP_A<-fgseaRes$leadingEdge[[2]]
rld <- vst(dds, blind=T)
normalized_data <- assay(rld)
colnames(normalized_data) <- colnames(data)

dogdata<-normalized_data[row.names(normalized_data)%in%LeadingUP_A,]
library(pheatmap)
library(RColorBrewer)
classdata<-data.frame("dog class" = dogclass$class)
row.names(classdata)<-row.names(dogclass)
breaks <- seq(-2, 2, length = 16)
color <- colorRampPalette(rev(brewer.pal(11,'RdBu')))(length(breaks))


pheatmap(dogdata, annotation_col = classdata,
         scale = "row", color = color, breaks = breaks, fontsize = 6)



#Analysis of leading up_C
LeadingUP_C<-fgseaRes$leadingEdge[[4]]

dogdata<-normalized_data[row.names(normalized_data)%in%LeadingUP_C,]
classdata<-data.frame("dog class" =dogclass$class)
row.names(classdata)<-row.names(dogclass)
breaks <- seq(-2, 2, length = 16)
color <- colorRampPalette(rev(brewer.pal(11,'RdBu')))(length(breaks))


pheatmap(dogdata, annotation_col = classdata,
         scale = "row", color = color, breaks = breaks, fontsize = 6)

#write fgsea results
fgseaRes<-data.frame(fgseaRes[,-8])
write.csv(fgseaRes, file = "fgsea_results.csv")
