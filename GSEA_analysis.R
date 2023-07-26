#We can consider all the genes that are specifically up and down regulated in classes A/B/C. 
#Next we convert those to mouse gene symbol. We create a list of genes (A_Up, A_down, B_Up...). 
#And then use those gene  sets in GSEA analysis to see which ones are up/down in mouse. 
#In mouse you can perform mouse A vs C group comparison and give that to FGSEA R package
#create up and down regulated gene sets with additional tumor class subcategories

#log2foldchange shows upregulation and downregulation; upregulation > 0 and downregulation < 0 
#intersect(siggenes,row.names(normalized_data_GS))
SignficantGenes <-read.csv("~/Desktop/PatelLab/ClassSpecificGenes_InhouseRaleighHK_031723.csv", row.names=1)
SignficantGenes<-SignficantGenes[,c(2,5,8)]
SignficantGenes<-SignficantGenes[!is.na(SignficantGenes$AvsOthers_log2FoldChange),]
SignficantGenes$GENESYMBOL<-row.names(SignficantGenes)
row.names(SignficantGenes)<-NULL
upA<-SignficantGenes[SignficantGenes$AvsOthers_log2FoldChange>1.5,]
upA$Class<-rep("up_A", nrow(upA))

downA<-SignficantGenes[SignficantGenes$AvsOthers_log2FoldChange<=-1.5,]
downA$Class<-rep("down_A", nrow(downA))


upB<-SignficantGenes[SignficantGenes$BvsOthers_log2FoldChange>=1.5,]
upB$Class<-rep("up_B", nrow(upB))

downB<-SignficantGenes[SignficantGenes$BvsOthers_log2FoldChange<=-1.5,]
downB$Class<-rep("down_B", nrow(downB))


upC<-SignficantGenes[SignficantGenes$CvsOthers_log2FoldChange>=1.5,]
upC$Class<-rep("up_C", nrow(upC))

downC<-SignficantGenes[SignficantGenes$CvsOthers_log2FoldChange<=-1.5,]
downC$Class<-rep("down_C", nrow(downC))

126+393+191+212+100+120
SignficantGenes<-rbind(upA,downA, upB,downB, upC,downC)


SignficantGenes<-SignficantGenes[!duplicated(SignficantGenes$GENESYMBOL)&!duplicated(SignficantGenes$GENESYMBOL, fromLast = TRUE),]

SignficantGenes<-SignficantGenes[,-c(1:3)]
table(SignficantGenes$Class)


#check to make sure there are unique genes
#duplicated(row.names(SigA), row.names(SiC))
#sum(duplicated(row.names(SigB), row.names(SigC))
#)

#Sigcom<-rbind(SigA, SigB)
#length(unique(row.names(Sigcom)))


load("~/Desktop/PatelLab/akash_dog_meningioma/human_Rdata/09192022_GE_RAW_Normalized_N330PrimarySamples_GeneSymbol.rda")
#Sampleclass<-data.frame("samplename"= row.names(samples),samples$class)
#class<-c("A","B","C")
#geneset<-""
#for(a in 1:3){
#  targetsamples<-Sampleclass[Sampleclass$samples.class==class[a],]
#  targetdata<-normalized_data_GS[,colnames(normalized_data_GS)%in%targetsamples$samplename]
#  referencesamples<-Sampleclass[Sampleclass$samples.class!=class[a],]
#  referencedata<-normalized_data_GS[,colnames(normalized_data_GS)%in%referencesamples$samplename]
#  siggenes<-row.names(SignficantGenes[SignficantGenes$Class==class[a],])
#targetgenemeans<-apply(targetdata[row.names(targetdata)%in%siggenes,],1, function(x) {
#   mean(x)})
# referencegenemeans<-apply(referencedata[row.names(referencedata)%in%siggenes,],1, function(x) {
#   mean(x)})
# for(b in 1:length(targetgenemeans)){
#   if(targetgenemeans[b]>referencegenemeans[b]){
#     geneset<-c(geneset, paste0("UP_",class[a]))
#   }else{
#     geneset<-c(geneset, paste0("DOWN_", class[a]))
#   }
# }
#}

#geneset<-geneset[-1]
write.table(table(SignficantGenes$Class), file = "geneset.csv")

geneset<-SignficantGenes
colnames(geneset)<-c("GENESYMBOL","geneset")

 #for(x in 1:length(siggenes)){
  #  if(mean(as.numeric(targetdata[which(row.names(targetdata)==siggenes[1]),]))>
  #     mean(referencedata[which(row.names(referencedata)==siggenes[x]),])){
   #   geneset<-c(geneset, print0("UP_",class[a]))
  #  }else (
  #    geneset<-c(geneset, print0("DOWN_", class[a]))
  #  )
#  }

#geneset<-data.frame(SignficantGenes,geneset)


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


test<-geneset[geneset$GENESYMBOL%in%names(ranks),]
table(test$geneset)
names(ranks)
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

geneSets
fgseaRes <- fgseaMultilevel(pathways = geneSets,
                            stats    = ranks)


#genesets is the human A_Up, A_down, B_Up... genesets

#Analysis of leading up_A genes between dogsA and dogsB
LeadingUP_A<-fgseaRes$leadingEdge[[4]]
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



#Analysis of leading down_B
Leadingdown_B<-fgseaRes$leadingEdge[[2]]

dogdata<-normalized_data[row.names(normalized_data)%in%Leadingdown_B,]
classdata<-data.frame("dog class" =dogclass$class)
row.names(classdata)<-row.names(dogclass)
breaks <- seq(-2, 2, length = 16)
color <- colorRampPalette(rev(brewer.pal(11,'RdBu')))(length(breaks))


pheatmap(dogdata, annotation_col = classdata,
         scale = "row", color = color, breaks = breaks, fontsize = 6)

#down_A leading genes
Leadingdown_A<-fgseaRes$leadingEdge[[1]]

dogdata<-normalized_data[row.names(normalized_data)%in%Leadingdown_A,]
classdata<-data.frame("dog class" =dogclass$class)
row.names(classdata)<-row.names(dogclass)
breaks <- seq(-2, 2, length = 16)
color <- colorRampPalette(rev(brewer.pal(11,'RdBu')))(length(breaks))


pheatmap(dogdata, annotation_col = classdata,
         scale = "row", color = color, breaks = breaks, fontsize = 6)

#down_C leading genes

Leadingdown_C<-fgseaRes$leadingEdge[[3]]

dogdata<-normalized_data[row.names(normalized_data)%in%Leadingdown_C,]
classdata<-data.frame("dog class" =dogclass$class)
row.names(classdata)<-row.names(dogclass)
breaks <- seq(-2, 2, length = 16)
color <- colorRampPalette(rev(brewer.pal(11,'RdBu')))(length(breaks))


pheatmap(dogdata, annotation_col = classdata,
         scale = "row", color = color, breaks = breaks, fontsize = 6)

#up_b leading genes
LeadingUp_B<-fgseaRes$leadingEdge[[5]]

dogdata<-normalized_data[row.names(normalized_data)%in%LeadingUp_B,]
classdata<-data.frame("dog class" =dogclass$class)
row.names(classdata)<-row.names(dogclass)
breaks <- seq(-2, 2, length = 16)
color <- colorRampPalette(rev(brewer.pal(11,'RdBu')))(length(breaks))


pheatmap(dogdata, annotation_col = classdata,
         scale = "row", color = color, breaks = breaks, fontsize = 6)
#up_C leading genes
LeadingUp_C<-fgseaRes$leadingEdge[[6]]

dogdata<-normalized_data[row.names(normalized_data)%in%LeadingUp_C,]
classdata<-data.frame("dog class" =dogclass$class)
row.names(classdata)<-row.names(dogclass)
breaks <- seq(-2, 2, length = 16)
color <- colorRampPalette(rev(brewer.pal(11,'RdBu')))(length(breaks))


pheatmap(dogdata, annotation_col = classdata,
         scale = "row", color = color, breaks = breaks, fontsize = 6)

#write fgsea results
#fgseaRes<-data.frame(fgseaRes[,-8])
#write.csv(fgseaRes, file = "fgsea_results.csv")
