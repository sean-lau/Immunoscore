##December Work
install.packages("BiocManager")
BiocManager::install("WGCNA")
library(WGCNA)
options(stringsAsFactors = FALSE)

load("~/Desktop/PatelLab/notumorseuratObj.Rdata")

#Matrix of data
immunecellscounts <- as.matrix(Seurat::GetAssayData(immune, slot = "counts"))

##Data Cleaning
gsg <- goodSamplesGenes(immunecellscounts, verbose = 3)
gsg$allOK

sampleTree <- hclust(dist(immunecellscounts), method = 'average')
sizeGrWindow(12,9)

par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, sub = "", xlab = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

##automatic network construction and module detection
#choosing a set of soft-thresholdinng powers
powers <- c(c(1:10), seq(from = 1, to = 20, by = 2))
#call newtork topology analysis function
sft <-pickSoftThreshold(immunecellscounts, powerVector = powers, verbose = 5)
###plot results
sizeGrWindow(9,5)
par(mfrow = c(1,2))
cex1 <- 0.9
?sizeGrWindow

?par()
#Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshld (power)", ylab = "Scale Free topology model Fit, signed R^2", type= "n", main = paste("scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

abline(h=0.90,col="red")
view(sft$fitIndices)

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


##Use 19 as the soft threshold power

